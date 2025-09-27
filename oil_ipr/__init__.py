"""oil_ipr: Oil well inflow performance relationships and CLI."""

__all__ = [
	"pi_rate",
	"pi_pwf",
	"vogel_qmax_from_point",
	"vogel_rate_from_qmax",
	"vogel_pwf_from_qmax",
	"vogel_two_point_calibrate",
	"fetkovich_rate",
	"fetkovich_pwf",
	"fetkovich_norm_rate",
	"fetkovich_norm_pwf",
]

__version__ = "0.1.0"

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