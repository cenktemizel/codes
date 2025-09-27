from __future__ import annotations

import math
from typing import Tuple


def _ensure_non_negative(name: str, value: float) -> None:
	if value < 0:
		raise ValueError(f"{name} must be non-negative; got {value}")


def _ensure_positive(name: str, value: float) -> None:
	if value <= 0:
		raise ValueError(f"{name} must be > 0; got {value}")


# -------------------------
# Linear Productivity Index
# -------------------------

def pi_rate(productivity_index: float, reservoir_pressure: float, flowing_pressure: float) -> float:
	"""Compute rate from linear PI model: q = J * (pr - pwf).

	Args:
		productivity_index: Productivity index J (rate per pressure).
		reservoir_pressure: Average reservoir pressure pr.
		flowing_pressure: Bottomhole flowing pressure pwf.

	Returns:
		Rate q.
	"""
	_ensure_positive("productivity_index", productivity_index)
	_ensure_non_negative("reservoir_pressure", reservoir_pressure)
	_ensure_non_negative("flowing_pressure", flowing_pressure)
	if flowing_pressure > reservoir_pressure:
		raise ValueError("flowing_pressure cannot exceed reservoir_pressure in PI model")
	return productivity_index * (reservoir_pressure - flowing_pressure)


def pi_pwf(productivity_index: float, reservoir_pressure: float, rate: float) -> float:
	"""Compute pwf from linear PI model: pwf = pr - q/J.

	Args:
		productivity_index: Productivity index J.
		reservoir_pressure: Average reservoir pressure pr.
		rate: Target rate q.

	Returns:
		Flowing pressure pwf.
	"""
	_ensure_positive("productivity_index", productivity_index)
	_ensure_non_negative("reservoir_pressure", reservoir_pressure)
	_ensure_non_negative("rate", rate)
	pwf = reservoir_pressure - rate / productivity_index
	if pwf < 0:
		pwf = 0.0
	return pwf


# ------
# Vogel
# ------

def vogel_qmax_from_point(reservoir_pressure: float, rate: float, flowing_pressure: float) -> float:
	"""Calibrate Vogel q_max from a single test point and pr.

	Vogel: q/q_max = 1 - 0.2 x - 0.8 x^2, where x = pwf/pr.
	Given (q, pwf, pr), q_max = q / f(x).
	"""
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_non_negative("rate", rate)
	_ensure_non_negative("flowing_pressure", flowing_pressure)
	if flowing_pressure > reservoir_pressure:
		raise ValueError("flowing_pressure cannot exceed reservoir_pressure in Vogel model")
	x = flowing_pressure / reservoir_pressure
	f = 1.0 - 0.2 * x - 0.8 * x * x
	if f <= 0:
		raise ValueError("Unphysical test point for Vogel (f <= 0). Check inputs.")
	return rate / f


def vogel_rate_from_qmax(reservoir_pressure: float, q_max: float, flowing_pressure: float) -> float:
	"""Compute rate from Vogel given q_max and pr."""
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_positive("q_max", q_max)
	_ensure_non_negative("flowing_pressure", flowing_pressure)
	if flowing_pressure > reservoir_pressure:
		raise ValueError("flowing_pressure cannot exceed reservoir_pressure in Vogel model")
	x = flowing_pressure / reservoir_pressure
	f = 1.0 - 0.2 * x - 0.8 * x * x
	return max(0.0, q_max * f)


def vogel_pwf_from_qmax(reservoir_pressure: float, q_max: float, rate: float) -> float:
	"""Compute pwf from Vogel given q_max and pr.

	Solve normalized quadratic: 0.8 y^2 + 0.2 y + (q/q_max - 1) = 0, y in [0, 1].
	"""
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_positive("q_max", q_max)
	_ensure_non_negative("rate", rate)
	alpha = 0.8
	beta = 0.2
	gamma = (rate / q_max) - 1.0
	discriminant = beta * beta - 4.0 * alpha * gamma
	if discriminant < 0:
		raise ValueError("No real solution for pwf; rate exceeds q_max or inputs invalid")
	sqrt_d = math.sqrt(discriminant)
	# We need y in [0, 1]; choose the physically meaningful root
	y1 = (-beta + sqrt_d) / (2.0 * alpha)
	y2 = (-beta - sqrt_d) / (2.0 * alpha)
	candidates = [y for y in (y1, y2) if 0.0 <= y <= 1.0]
	if not candidates:
		raise ValueError("No physical pwf in [0, pr] for given rate")
	y = max(candidates)
	return reservoir_pressure * y


def vogel_two_point_calibrate(q1: float, pwf1: float, q2: float, pwf2: float) -> Tuple[float, float]:
	"""Calibrate (pr, q_max) from two test points using Vogel.

	Solves for pr > max(pwf_i) such that q1/f(x1) == q2/f(x2),
	where f(x) = 1 - 0.2x - 0.8x^2 and x_i = pwf_i/pr.
	Returns (pr, q_max).
	"""
	_ensure_non_negative("q1", q1)
	_ensure_non_negative("q2", q2)
	_ensure_non_negative("pwf1", pwf1)
	_ensure_non_negative("pwf2", pwf2)
	pwf_max = max(pwf1, pwf2)
	if q1 == q2 and pwf1 == pwf2:
		raise ValueError("Two-point calibration requires distinct points")

	def f(x: float) -> float:
		return 1.0 - 0.2 * x - 0.8 * x * x

	def g(pr: float) -> float:
		x1 = pwf1 / pr
		x2 = pwf2 / pr
		if not (0.0 <= x1 <= 1.0 and 0.0 <= x2 <= 1.0):
			# Force sign to push search higher
			return -1.0
		return (q1 / f(x1)) - (q2 / f(x2))

	# Bracket a root for g(pr) = 0 with pr > pwf_max
	lower = max(pwf_max * 1.000001, pwf_max + 1e-6)
	upper = max(10.0 * lower, lower + 100.0)
	g_low = g(lower)
	g_up = g(upper)
	attempts = 0
	while g_low * g_up > 0 and attempts < 30:
		upper *= 2.0
		g_up = g(upper)
		attempts += 1
	if g_low * g_up > 0:
		raise ValueError("Failed to bracket reservoir pressure for two-point Vogel calibration")

	# Bisection search
	for _ in range(80):
		mid = 0.5 * (lower + upper)
		g_mid = g(mid)
		if abs(g_mid) < 1e-9:
			lower = upper = mid
			break
		if g_low * g_mid <= 0:
			upper = mid
			g_up = g_mid
		else:
			lower = mid
			g_low = g_mid

	pr = 0.5 * (lower + upper)
	x1 = pwf1 / pr
	x2 = pwf2 / pr
	f1 = f(x1)
	f2 = f(x2)
	qmax1 = q1 / f1
	qmax2 = q2 / f2
	q_max = 0.5 * (qmax1 + qmax2)
	return pr, q_max


# ----------
# Fetkovich
# ----------

def fetkovich_rate(productivity_index: float, reservoir_pressure: float, flowing_pressure: float, exponent_n: float) -> float:
	"""Power-law PI (Fetkovich) for undersaturated reservoirs: q = J * (pr - pwf)^n."""
	_ensure_positive("productivity_index", productivity_index)
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_non_negative("flowing_pressure", flowing_pressure)
	_ensure_positive("exponent_n", exponent_n)
	if flowing_pressure > reservoir_pressure:
		raise ValueError("flowing_pressure cannot exceed reservoir_pressure in Fetkovich power-law model")
	return productivity_index * math.pow(reservoir_pressure - flowing_pressure, exponent_n)


def fetkovich_pwf(productivity_index: float, reservoir_pressure: float, rate: float, exponent_n: float) -> float:
	"""Invert power-law PI to get pwf for a target rate."""
	_ensure_positive("productivity_index", productivity_index)
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_non_negative("rate", rate)
	_ensure_positive("exponent_n", exponent_n)
	if rate == 0:
		return reservoir_pressure
	delta_p = math.pow(rate / productivity_index, 1.0 / exponent_n)
	pwf = reservoir_pressure - delta_p
	return max(0.0, pwf)


def fetkovich_norm_rate(reservoir_pressure: float, q_max: float, flowing_pressure: float, exponent_n: float) -> float:
	"""Normalized Fetkovich: q = q_max * (1 - (pwf/pr)^n)."""
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_positive("q_max", q_max)
	_ensure_non_negative("flowing_pressure", flowing_pressure)
	_ensure_positive("exponent_n", exponent_n)
	if flowing_pressure > reservoir_pressure:
		raise ValueError("flowing_pressure cannot exceed reservoir_pressure in normalized Fetkovich")
	x = flowing_pressure / reservoir_pressure
	return max(0.0, q_max * (1.0 - math.pow(x, exponent_n)))


def fetkovich_norm_pwf(reservoir_pressure: float, q_max: float, rate: float, exponent_n: float) -> float:
	"""Invert normalized Fetkovich to get pwf for a target rate."""
	_ensure_positive("reservoir_pressure", reservoir_pressure)
	_ensure_positive("q_max", q_max)
	_ensure_non_negative("rate", rate)
	_ensure_positive("exponent_n", exponent_n)
	if rate > q_max:
		raise ValueError("rate cannot exceed q_max in normalized Fetkovich")
	if rate == 0:
		return reservoir_pressure
	y = 1.0 - (rate / q_max)
	pwf = reservoir_pressure * math.pow(y, 1.0 / exponent_n)
	return max(0.0, min(reservoir_pressure, pwf))