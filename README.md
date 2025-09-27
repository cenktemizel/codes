## Oil Well Inflow Performance Calculator (IPR)

Command-line tool to compute oil well inflow performance using:
- Darcy linear IPR (via productivity index `J`)
- Vogel IPR (solution-gas drive) with `q_max` from a test point or slope match
- Piecewise undersaturated method (linear above `P_b`, Vogel below)

File: `oil_ipr.py`

### Installation
- Requires Python 3.8+
- No external packages needed

### Usage
Run with `--help` to see options:

```bash
python3 oil_ipr.py --help
```

Key arguments:
- `--method`: `auto` | `darcy` | `vogel` | `piecewise` (default `auto`)
- `--pr`: average reservoir pressure, psi
- `--pb`: bubble point pressure, psi (needed for `piecewise` or `auto` with `Pr > Pb`)
- `--J`: productivity index, STB/day/psi
- or provide properties to compute `J`: `--k --h --mu --B --re --rw [--s]`
- Vogel calibration: `--q-test` and `--pwf-test` (or `J` for slope-matched `q_max`)
- Single-point: `--pwf`
- Curve generation: `--pmin --pmax --pstep`
- CSV export: `--csv path.csv`
- AOF printout: `--aof`

### Examples

1) Darcy single-point rate, given `J`:
```bash
python3 oil_ipr.py --method darcy --pr 3000 --J 1.2 --pwf 1500
```

2) Vogel with one test point, AOF at `Pwf=0`:
```bash
python3 oil_ipr.py --method vogel --pr 3200 --q-test 1200 --pwf-test 1200 --pwf 0 --aof
```

3) Piecewise curve with CSV export (Pr > Pb):
```bash
python3 oil_ipr.py --method piecewise --pr 4000 --pb 2500 --J 1.0 \
  --pmin 0 --pmax 4000 --pstep 1000 --csv ipr_piecewise.csv --aof
```

### Notes
- Units assumed (field): psi, md, ft, cP, rb/STB, STB/day
- `J` formula used when properties supplied:
  `J = 0.00708 * k * h / ( mu * B * ( ln(re/rw) - 0.75 + s ) )`
- Vogel fraction: `f = 1 - 0.2*y - 0.8*y^2`, where `y = Pwf/Pr`
- `q_max` from slope match: `q_max = J * Pr / 1.8`
- Piecewise: linear for `Pwf >= Pb`; below `Pb`, Vogel with `q_max` fitted at `Pwf=Pb`

### Disclaimer
Engineering correlations carry assumptions. Validate inputs and results with field data.
