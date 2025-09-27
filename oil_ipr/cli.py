import argparse
import csv
import sys
from typing import Optional

from .models import (
	pi_rate,
	pi_pwf,
	vogel_qmax_from_point,
	vogel_rate_from_qmax,
	vogel_pwf_from_qmax,
	vogel_two_point_calibrate,
	fetkovich_rate,
	fetkovich_pwf,
	fetkovich_norm_rate,
	fetkovich_norm_pwf,
)


def _float(value: str) -> float:
	try:
		return float(value)
	except Exception as exc:
		raise argparse.ArgumentTypeError(str(exc)) from exc


def _build_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(prog="oil-ipr", description="Oil well inflow performance calculator")
	sub = parser.add_subparsers(dest="command", required=True)

	# rate
	rate = sub.add_parser("rate", help="Compute rate q for a given pwf")
	_rate_models = rate.add_subparsers(dest="model", required=True)

	pi = _rate_models.add_parser("pi", help="Linear PI: q = J (pr - pwf)")
	pi.add_argument("--J", type=_float, required=True, help="Productivity index")
	pi.add_argument("--pr", type=_float, required=True, help="Reservoir pressure")
	pi.add_argument("--pwf", type=_float, required=True, help="Flowing pressure")

	vogel = _rate_models.add_parser("vogel", help="Vogel with single-point calibration or qmax")
	vogel.add_argument("--pr", type=_float, required=True, help="Reservoir pressure")
	vogel.add_argument("--pwf", type=_float, required=True, help="Flowing pressure")
	vogel_group = vogel.add_mutually_exclusive_group(required=True)
	vogel_group.add_argument("--qmax", type=_float, help="Maximum rate at pwf=0")
	vogel_group.add_argument("--q-test", type=_float, help="Single test rate for calibration")
	vogel.add_argument("--pwf-test", type=_float, help="Single test pwf for calibration (required with --q-test)")

	vogel2 = _rate_models.add_parser("vogel2", help="Vogel two-point calibration then rate")
	vogel2.add_argument("--q1", type=_float, required=True)
	vogel2.add_argument("--pwf1", type=_float, required=True)
	vogel2.add_argument("--q2", type=_float, required=True)
	vogel2.add_argument("--pwf2", type=_float, required=True)
	vogel2.add_argument("--pwf", type=_float, required=True, help="Flowing pressure for prediction")
	vogel2.add_argument("--print-calibration", action="store_true", help="Print pr and qmax used")

	fetk = _rate_models.add_parser("fetkovich", help="Fetkovich power-law: q = J (pr - pwf)^n")
	fetk.add_argument("--J", type=_float, required=True)
	fetk.add_argument("--pr", type=_float, required=True)
	fetk.add_argument("--pwf", type=_float, required=True)
	fetk.add_argument("--n", type=_float, required=True)

	fetkn = _rate_models.add_parser("fetkovich-norm", help="Normalized Fetkovich: q = qmax (1 - (pwf/pr)^n)")
	fetkn.add_argument("--pr", type=_float, required=True)
	fetkn.add_argument("--pwf", type=_float, required=True)
	fetkn.add_argument("--qmax", type=_float, required=True)
	fetkn.add_argument("--n", type=_float, required=True)

	# pwf
	pwf = sub.add_parser("pwf", help="Compute pwf for a given rate")
	_pwf_models = pwf.add_subparsers(dest="model", required=True)

	pi_p = _pwf_models.add_parser("pi")
	pi_p.add_argument("--J", type=_float, required=True)
	pi_p.add_argument("--pr", type=_float, required=True)
	pi_p.add_argument("--q", type=_float, required=True)

	vogel_p = _pwf_models.add_parser("vogel")
	vogel_p.add_argument("--pr", type=_float, required=True)
	vogel_p_group = vogel_p.add_mutually_exclusive_group(required=True)
	vogel_p_group.add_argument("--qmax", type=_float)
	vogel_p_group.add_argument("--q-test", type=_float)
	vogel_p.add_argument("--pwf-test", type=_float, help="Required with --q-test")
	vogel_p.add_argument("--q", type=_float, required=True)

	fetk_p = _pwf_models.add_parser("fetkovich")
	fetk_p.add_argument("--J", type=_float, required=True)
	fetk_p.add_argument("--pr", type=_float, required=True)
	fetk_p.add_argument("--n", type=_float, required=True)
	fetk_p.add_argument("--q", type=_float, required=True)

	fetkn_p = _pwf_models.add_parser("fetkovich-norm")
	fetkn_p.add_argument("--pr", type=_float, required=True)
	afetkn_group = fetkn_p.add_mutually_exclusive_group(required=True)
	afetkn_group.add_argument("--qmax", type=_float)
	afetkn_group.add_argument("--q-test", type=_float)
	fetkn_p.add_argument("--pwf-test", type=_float, help="Required with --q-test")
	fetkn_p.add_argument("--n", type=_float, required=True)
	fetkn_p.add_argument("--q", type=_float, required=True)

	# table
	table = sub.add_parser("table", help="Generate pwf vs rate table as CSV")
	_table_models = table.add_subparsers(dest="model", required=True)

	for name in ("pi", "vogel", "fetkovich", "fetkovich-norm"):
		m = _table_models.add_parser(name)
		m.add_argument("--pwf-min", type=_float, required=True)
		m.add_argument("--pwf-max", type=_float, required=True)
		m.add_argument("--steps", type=int, default=25)
		m.add_argument("--csv", type=str, default="-", help="Output CSV path or '-' for stdout")
		if name == "pi":
			m.add_argument("--J", type=_float, required=True)
			m.add_argument("--pr", type=_float, required=True)
		elif name == "vogel":
			m.add_argument("--pr", type=_float, required=True)
			mg = m.add_mutually_exclusive_group(required=True)
			mg.add_argument("--qmax", type=_float)
			mg.add_argument("--q-test", type=_float)
			m.add_argument("--pwf-test", type=_float)
		elif name == "fetkovich":
			m.add_argument("--J", type=_float, required=True)
			m.add_argument("--pr", type=_float, required=True)
			m.add_argument("--n", type=_float, required=True)
		elif name == "fetkovich-norm":
			m.add_argument("--pr", type=_float, required=True)
			m.add_argument("--qmax", type=_float, required=True)
			m.add_argument("--n", type=_float, required=True)

	return parser


def _write_csv(rows, path: str) -> None:
	fieldnames = ["pwf", "q"]
	if path == "-":
		writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
		writer.writeheader()
		for row in rows:
			writer.writerow(row)
	else:
		with open(path, "w", newline="") as f:
			writer = csv.DictWriter(f, fieldnames=fieldnames)
			writer.writeheader()
			for row in rows:
				writer.writerow(row)


def _handle_rate(args: argparse.Namespace) -> int:
	if args.model == "pi":
		q = pi_rate(args.J, args.pr, args.pwf)
		print(f"q = {q}")
		return 0
	if args.model == "vogel":
		if args.qmax is not None:
			q = vogel_rate_from_qmax(args.pr, args.qmax, args.pwf)
		else:
			if args.pwf_test is None:
				raise SystemExit("--pwf-test is required with --q-test")
			qmax = vogel_qmax_from_point(args.pr, args.q_test, args.pwf_test)
			q = vogel_rate_from_qmax(args.pr, qmax, args.pwf)
		print(f"q = {q}")
		return 0
	if args.model == "vogel2":
		pr, qmax = vogel_two_point_calibrate(args.q1, args.pwf1, args.q2, args.pwf2)
		q = vogel_rate_from_qmax(pr, qmax, args.pwf)
		if args.print_calibration:
			print(f"pr = {pr}")
			print(f"qmax = {qmax}")
		print(f"q = {q}")
		return 0
	if args.model == "fetkovich":
		q = fetkovich_rate(args.J, args.pr, args.pwf, args.n)
		print(f"q = {q}")
		return 0
	if args.model == "fetkovich-norm":
		q = fetkovich_norm_rate(args.pr, args.qmax, args.pwf, args.n)
		print(f"q = {q}")
		return 0
	raise SystemExit("Unknown model")


def _handle_pwf(args: argparse.Namespace) -> int:
	if args.model == "pi":
		pwf = pi_pwf(args.J, args.pr, args.q)
		print(f"pwf = {pwf}")
		return 0
	if args.model == "vogel":
		if args.qmax is not None:
			pwf = vogel_pwf_from_qmax(args.pr, args.qmax, args.q)
		else:
			if args.pwf_test is None:
				raise SystemExit("--pwf-test is required with --q-test")
			qmax = vogel_qmax_from_point(args.pr, args.q_test, args.pwf_test)
			pwf = vogel_pwf_from_qmax(args.pr, qmax, args.q)
		print(f"pwf = {pwf}")
		return 0
	if args.model == "fetkovich":
		pwf = fetkovich_pwf(args.J, args.pr, args.q, args.n)
		print(f"pwf = {pwf}")
		return 0
	if args.model == "fetkovich-norm":
		if args.qmax is not None:
			pwf = fetkovich_norm_pwf(args.pr, args.qmax, args.q, args.n)
		else:
			if args.pwf_test is None:
				raise SystemExit("--pwf-test is required with --q-test")
			qmax = vogel_qmax_from_point(args.pr, args.q_test, args.pwf_test)
			pwf = fetkovich_norm_pwf(args.pr, qmax, args.q, args.n)
		print(f"pwf = {pwf}")
		return 0
	raise SystemExit("Unknown model")


def _handle_table(args: argparse.Namespace) -> int:
	pwf_min = args["pwf_min"] if isinstance(args, dict) else args.pwf_min
	pwf_max = args["pwf_max"] if isinstance(args, dict) else args.pwf_max
	steps = args["steps"] if isinstance(args, dict) else args.steps
	if steps <= 1:
		raise SystemExit("--steps must be > 1")
	rows = []
	for i in range(steps + 1):
		pwf = pwf_min + (pwf_max - pwf_min) * (i / steps)
		if args.model == "pi":
			q = pi_rate(args.J, args.pr, pwf)
		elif args.model == "vogel":
			if getattr(args, "qmax", None) is not None:
				q = vogel_rate_from_qmax(args.pr, args.qmax, pwf)
			else:
				qmax = vogel_qmax_from_point(args.pr, args.q_test, args.pwf_test)
				q = vogel_rate_from_qmax(args.pr, qmax, pwf)
		elif args.model == "fetkovich":
			q = fetkovich_rate(args.J, args.pr, pwf, args.n)
		elif args.model == "fetkovich-norm":
			q = fetkovich_norm_rate(args.pr, args.qmax, pwf, args.n)
		else:
			raise SystemExit("Unknown model")
		rows.append({"pwf": pwf, "q": q})
	_write_csv(rows, args.csv)
	return 0


def main(argv: Optional[list[str]] = None) -> int:
	parser = _build_parser()
	args = parser.parse_args(argv)
	if args.command == "rate":
		return _handle_rate(args)
	if args.command == "pwf":
		return _handle_pwf(args)
	if args.command == "table":
		return _handle_table(args)
	return 1


if __name__ == "__main__":
	sys.exit(main())