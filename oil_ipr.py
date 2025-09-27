#!/usr/bin/env python3
"""
Oil Well Inflow Performance Relationship (IPR) Calculator

Supports:
- Darcy linear IPR using productivity index J (from input or computed from reservoir properties)
- Vogel IPR (solution-gas drive) with q_max estimated from a test point or from J
- Piecewise undersaturated method: linear above P_b and Vogel below P_b

Units: Field units are assumed throughout unless otherwise noted.
- Pressure: psi
- Permeability k: md
- Thickness h: ft
- Viscosity mu: cP
- Formation volume factor B: rb/STB
- Radii r_e, r_w: ft
- Rates: STB/day

Formula references:
- Darcy productivity index (steady-state radial):
  J = 0.00708 * k * h / [ mu * B * ( ln(r_e / r_w) - 0.75 + s ) ]

- Darcy IPR: q = J * (P_r - P_wf)

- Vogel (1978) saturated oil IPR:
  q / q_max = 1 - 0.2 * (P_wf / P_r) - 0.8 * (P_wf / P_r)^2

- q_max from one test (q_t at P_wf,t):
  q_max = q_t / [ 1 - 0.2 * (P_wf,t / P_r) - 0.8 * (P_wf,t / P_r)^2 ]

- q_max from J near P_r (matching slope at small drawdown):
  dq/dP_wf|_{P_wf=P_r} = -J = -1.8 * q_max / P_r  =>  q_max = J * P_r / 1.8

- Piecewise undersaturated (P_r > P_b):
  For P_wf >= P_b: q = J * (P_r - P_wf)
  Let q_b = J * (P_r - P_b)
  Compute q_max via Vogel at P_wf = P_b with reservoir pressure P_r:
    q_max = q_b / [ 1 - 0.2 * (P_b / P_r) - 0.8 * (P_b / P_r)^2 ]
  For P_wf <= P_b: use Vogel with that q_max and P_r

Author: CLI utility intended for quick engineering calculations.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple


# ------------------------- Core Calculations -------------------------


def compute_productivity_index(
    permeability_md: float,
    thickness_ft: float,
    viscosity_cp: float,
    formation_volume_factor_rb_per_stb: float,
    drainage_radius_ft: float,
    wellbore_radius_ft: float,
    skin: float = 0.0,
) -> float:
    """Compute productivity index J [STB/day/psi] using steady-state radial Darcy flow.

    Field units formula:
    J = 0.00708 * k * h / ( mu * B * ( ln(re/rw) - 0.75 + s ) )

    Raises ValueError for invalid inputs.
    """
    if permeability_md <= 0 or thickness_ft <= 0:
        raise ValueError("Permeability and thickness must be positive")
    if viscosity_cp <= 0 or formation_volume_factor_rb_per_stb <= 0:
        raise ValueError("Viscosity and formation volume factor must be positive")
    if drainage_radius_ft <= 0 or wellbore_radius_ft <= 0:
        raise ValueError("Radii must be positive")
    if drainage_radius_ft <= wellbore_radius_ft:
        raise ValueError("Drainage radius must exceed wellbore radius")

    log_term = math.log(drainage_radius_ft / wellbore_radius_ft)
    denominator = (
        viscosity_cp
        * formation_volume_factor_rb_per_stb
        * (log_term - 0.75 + skin)
    )
    if denominator <= 0:
        raise ValueError(
            "Invalid denominator in J calculation. Check mu, B, radii, and skin values."
        )

    J = 0.00708 * permeability_md * thickness_ft / denominator
    return J


def darcy_rate(reservoir_pressure_psi: float, pwf_psi: float, productivity_index: float) -> float:
    """Compute oil rate [STB/day] using Darcy linear IPR: q = J * (Pr - Pwf)."""
    drawdown = reservoir_pressure_psi - pwf_psi
    return max(0.0, productivity_index * max(0.0, drawdown))


def vogel_fraction(pwf_over_pr: float) -> float:
    """Return the Vogel normalized fraction f = 1 - 0.2*y - 0.8*y^2, where y = Pwf/Pr.

    The returned value is clamped at >= 0.
    """
    f = 1.0 - 0.2 * pwf_over_pr - 0.8 * pwf_over_pr * pwf_over_pr
    return max(0.0, f)


def vogel_qmax_from_test(
    reservoir_pressure_psi: float, q_test_stb_per_day: float, pwf_test_psi: float
) -> float:
    """Estimate q_max for Vogel from a single test point at known Pwf.

    q_max = q_test / f(Pwf/Pr)
    """
    if reservoir_pressure_psi <= 0:
        raise ValueError("Reservoir pressure must be positive")
    if q_test_stb_per_day <= 0:
        raise ValueError("Test rate must be positive")
    if not (0 <= pwf_test_psi <= reservoir_pressure_psi):
        raise ValueError("Test Pwf must be within [0, Pr]")

    frac = vogel_fraction(pwf_test_psi / reservoir_pressure_psi)
    if frac <= 0:
        raise ValueError("Test point yields zero Vogel fraction; check inputs")
    return q_test_stb_per_day / frac


def vogel_qmax_from_pi(reservoir_pressure_psi: float, productivity_index: float) -> float:
    """Estimate q_max for Vogel by matching slope at Pwf = Pr: q_max = J * Pr / 1.8."""
    if reservoir_pressure_psi <= 0 or productivity_index <= 0:
        raise ValueError("Reservoir pressure and J must be positive")
    return productivity_index * reservoir_pressure_psi / 1.8


def vogel_rate_qmax(
    reservoir_pressure_psi: float, pwf_psi: float, qmax_stb_per_day: float
) -> float:
    """Compute oil rate [STB/day] via Vogel with known q_max.

    q = q_max * (1 - 0.2 * (Pwf/Pr) - 0.8 * (Pwf/Pr)^2)
    """
    if reservoir_pressure_psi <= 0 or qmax_stb_per_day < 0:
        raise ValueError("Reservoir pressure must be positive and q_max non-negative")
    if pwf_psi < 0:
        raise ValueError("Pwf must be non-negative")
    y = min(1.0, max(0.0, pwf_psi / reservoir_pressure_psi))
    return max(0.0, qmax_stb_per_day * vogel_fraction(y))


def piecewise_rate(
    reservoir_pressure_psi: float,
    bubble_point_psi: float,
    pwf_psi: float,
    productivity_index: float,
) -> Tuple[float, float]:
    """Compute rate using piecewise (undersaturated) method.

    Returns (q_stb_per_day, q_max_stb_per_day).
    """
    if bubble_point_psi <= 0:
        raise ValueError("Bubble point must be positive")
    if reservoir_pressure_psi <= bubble_point_psi:
        raise ValueError("Piecewise requires Pr > Pb")
    if productivity_index <= 0:
        raise ValueError("Productivity index J must be positive")

    q_b = darcy_rate(reservoir_pressure_psi, bubble_point_psi, productivity_index)
    # Compute q_max using Vogel at Pwf = Pb with reservoir pressure Pr
    denom = vogel_fraction(bubble_point_psi / reservoir_pressure_psi)
    if denom <= 0:
        raise ValueError("Invalid Vogel denominator at Pb; check pressures")
    q_max = q_b / denom

    if pwf_psi >= bubble_point_psi:
        q = darcy_rate(reservoir_pressure_psi, pwf_psi, productivity_index)
        return (q, q_max)
    else:
        q = vogel_rate_qmax(reservoir_pressure_psi, pwf_psi, q_max)
        return (q, q_max)


# ------------------------- Curve Generators -------------------------


def generate_pressure_range(pmin: float, pmax: float, pstep: float) -> List[float]:
    if pstep <= 0:
        raise ValueError("Pressure step must be positive")
    if pmin > pmax:
        raise ValueError("pmin cannot exceed pmax")
    values: List[float] = []
    p = pmin
    # Include pmax, ensure finite iterations
    max_iters = 1000000
    iters = 0
    while p <= pmax + 1e-9:
        values.append(round(p, 6))
        p += pstep
        iters += 1
        if iters > max_iters:
            raise RuntimeError("Pressure range iteration overflow; check step size")
    if values[-1] != pmax:
        values.append(pmax)
    return values


def curve_darcy(
    reservoir_pressure_psi: float,
    productivity_index: float,
    pmin: float,
    pmax: float,
    pstep: float,
) -> List[Tuple[float, float]]:
    pwfs = generate_pressure_range(pmin, pmax, pstep)
    return [
        (pwf, darcy_rate(reservoir_pressure_psi, pwf, productivity_index)) for pwf in pwfs
    ]


def curve_vogel(
    reservoir_pressure_psi: float,
    qmax_stb_per_day: float,
    pmin: float,
    pmax: float,
    pstep: float,
) -> List[Tuple[float, float]]:
    pwfs = generate_pressure_range(pmin, pmax, pstep)
    return [
        (pwf, vogel_rate_qmax(reservoir_pressure_psi, pwf, qmax_stb_per_day))
        for pwf in pwfs
    ]


def curve_piecewise(
    reservoir_pressure_psi: float,
    bubble_point_psi: float,
    productivity_index: float,
    pmin: float,
    pmax: float,
    pstep: float,
) -> Tuple[List[Tuple[float, float]], float]:
    pwfs = generate_pressure_range(pmin, pmax, pstep)
    rows: List[Tuple[float, float]] = []
    q_max_final: float = 0.0
    for pwf in pwfs:
        q, q_max = piecewise_rate(
            reservoir_pressure_psi, bubble_point_psi, pwf, productivity_index
        )
        rows.append((pwf, q))
        q_max_final = q_max
    return rows, q_max_final


# ------------------------- CSV Writer -------------------------


def write_csv(path: str, rows: Sequence[Tuple[float, float]]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pwf_psi", "q_stb_per_day"])
        for pwf, q in rows:
            writer.writerow([f"{pwf:.6f}", f"{q:.6f}"])


# ------------------------- CLI -------------------------


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Oil well IPR calculator (Darcy, Vogel, Piecewise)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument(
        "--method",
        choices=["auto", "darcy", "vogel", "piecewise"],
        default="auto",
        help="Calculation method",
    )
    p.add_argument("--pr", type=float, required=True, help="Average reservoir pressure, psi")
    p.add_argument("--pb", type=float, default=None, help="Bubble point pressure, psi")

    # Productivity index (direct) or properties to compute J
    p.add_argument("--J", type=float, default=None, help="Productivity index, STB/day/psi")
    p.add_argument("--k", type=float, help="Permeability, md")
    p.add_argument("--h", type=float, help="Net pay thickness, ft")
    p.add_argument("--mu", type=float, help="Oil viscosity, cP")
    p.add_argument("--B", type=float, help="Oil formation volume factor, rb/STB")
    p.add_argument("--re", type=float, help="Drainage radius, ft")
    p.add_argument("--rw", type=float, help="Wellbore radius, ft")
    p.add_argument("--s", type=float, default=0.0, help="Skin factor")

    # Vogel calibration options
    p.add_argument("--q-test", type=float, dest="q_test", help="Test rate, STB/day")
    p.add_argument("--pwf-test", type=float, dest="pwf_test", help="Test Pwf, psi")

    # Output controls
    p.add_argument("--pwf", type=float, default=None, help="Single-point Pwf, psi")
    p.add_argument("--pmin", type=float, default=0.0, help="Curve min Pwf, psi")
    p.add_argument(
        "--pmax",
        type=float,
        default=None,
        help="Curve max Pwf, psi (defaults to Pr)",
    )
    p.add_argument("--pstep", type=float, default=50.0, help="Curve Pwf step, psi")
    p.add_argument("--csv", type=str, default=None, help="Write CSV to path")
    p.add_argument("--aof", action="store_true", help="Print AOF estimate if available")

    return p.parse_args(argv)


def resolve_productivity_index(args: argparse.Namespace) -> float | None:
    if args.J is not None:
        if args.J <= 0:
            raise ValueError("J must be positive")
        return args.J

    # If any of the properties are provided, require all necessary ones
    prop_fields = [args.k, args.h, args.mu, args.B, args.re, args.rw]
    if any(v is not None for v in prop_fields):
        missing = {
            name
            for name, v in zip(
                ["k", "h", "mu", "B", "re", "rw"], prop_fields
            )
            if v is None
        }
        if missing:
            raise ValueError(
                f"Missing properties to compute J: {', '.join(sorted(missing))}"
            )
        return compute_productivity_index(
            permeability_md=float(args.k),
            thickness_ft=float(args.h),
            viscosity_cp=float(args.mu),
            formation_volume_factor_rb_per_stb=float(args.B),
            drainage_radius_ft=float(args.re),
            wellbore_radius_ft=float(args.rw),
            skin=float(args.s),
        )

    return None


def auto_select_method(args: argparse.Namespace, J: float | None) -> str:
    method = args.method
    if method != "auto":
        return method

    if args.pb is not None and args.pr > args.pb:
        return "piecewise"
    # If test is provided, prefer Vogel
    if args.q_test is not None and args.pwf_test is not None:
        return "vogel"
    # If J is provided, linear is a safe default. Vogel is possible with slope match, but be explicit only if requested.
    return "darcy"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    if args.pr <= 0:
        raise SystemExit("Reservoir pressure Pr must be positive")

    pmax = args.pmax if args.pmax is not None else args.pr
    if args.pwf is not None:
        if args.pwf < 0 or args.pwf > args.pr:
            raise SystemExit("Single-point Pwf must be within [0, Pr]")
    if args.pmin < 0 or pmax < 0 or args.pmin > pmax:
        raise SystemExit("Invalid curve Pwf bounds")

    # Resolve productivity index if needed
    J = resolve_productivity_index(args)

    # Determine method
    method = auto_select_method(args, J)

    qmax: float | None = None
    rows: List[Tuple[float, float]] = []

    if method == "darcy":
        if J is None:
            raise SystemExit(
                "Darcy method requires J or properties (k,h,mu,B,re,rw[,s])"
            )

        if args.pwf is not None:
            q = darcy_rate(args.pr, args.pwf, J)
            rows = [(args.pwf, q)]
        else:
            rows = curve_darcy(args.pr, J, args.pmin, pmax, args.pstep)

    elif method == "vogel":
        if args.q_test is not None and args.pwf_test is not None:
            qmax = vogel_qmax_from_test(args.pr, args.q_test, args.pwf_test)
        elif J is not None:
            qmax = vogel_qmax_from_pi(args.pr, J)
        else:
            raise SystemExit(
                "Vogel requires either a test point (q-test, pwf-test) or J to estimate q_max"
            )

        if args.pwf is not None:
            q = vogel_rate_qmax(args.pr, args.pwf, qmax)
            rows = [(args.pwf, q)]
        else:
            rows = curve_vogel(args.pr, qmax, args.pmin, pmax, args.pstep)

    elif method == "piecewise":
        if args.pb is None or args.pr <= args.pb:
            raise SystemExit("Piecewise requires Pb provided and Pr > Pb")
        if J is None:
            raise SystemExit(
                "Piecewise method requires J or properties (k,h,mu,B,re,rw[,s])"
            )

        if args.pwf is not None:
            q, qmax_val = piecewise_rate(args.pr, args.pb, args.pwf, J)
            qmax = qmax_val
            rows = [(args.pwf, q)]
        else:
            rows, qmax_val = curve_piecewise(
                args.pr, args.pb, J, args.pmin, pmax, args.pstep
            )
            qmax = qmax_val

    else:
        raise SystemExit(f"Unknown method: {method}")

    # Output
    if args.csv:
        write_csv(args.csv, rows)
        print(f"Wrote CSV: {args.csv} ({len(rows)} rows)")

    # Display table (brief)
    if len(rows) == 1:
        pwf, q = rows[0]
        print(f"Method: {method}")
        if J is not None:
            print(f"J = {J:.6f} STB/day/psi")
        if qmax is not None and args.aof:
            print(f"AOF = {qmax:.6f} STB/day")
        print(f"Pwf = {pwf:.2f} psi -> q = {q:.6f} STB/day")
    else:
        if J is not None:
            print(f"J = {J:.6f} STB/day/psi")
        if qmax is not None and args.aof:
            print(f"AOF = {qmax:.6f} STB/day")
        print("Pwf_psi,q_stb_per_day")
        # Print first few and last few for brevity
        preview = 5
        for i, (pwf, q) in enumerate(rows):
            if i < preview or i >= max(0, len(rows) - preview):
                print(f"{pwf:.2f},{q:.6f}")
        if len(rows) > 2 * preview:
            print(f"... ({len(rows)} rows total)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

