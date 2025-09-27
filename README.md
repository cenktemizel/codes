# oil-ipr

A Python CLI to compute oil well inflow performance (IPR) using common models:

- Linear productivity index (PI)
- Vogel (single-point and two-point calibration)
- Fetkovich (power-law) and normalized Fetkovich

## Install (editable)

```bash
pip install -e .
```

## CLI usage

- Compute rate using linear PI:

```bash
oil-ipr rate pi --J 1.8 --pr 3000 --pwf 1500
```

- Compute rate using Vogel with a single test point (calibrates qmax):

```bash
oil-ipr rate vogel --pr 3000 --pwf 800 --q-test 1200 --pwf-test 1500
```

- Compute rate using Vogel two-point (calibrates pr and qmax from two points):

```bash
oil-ipr rate vogel2 --q1 1200 --pwf1 1500 --q2 800 --pwf2 2200 --pwf 1000
```

- Compute pwf for a target rate using normalized Fetkovich (needs qmax, n, pr):

```bash
oil-ipr pwf fetkovich-norm --pr 3000 --qmax 2500 --n 0.65 --q 900
```

- Export a pwf–q table to CSV for plotting:

```bash
oil-ipr table vogel --pr 3000 --q-test 1200 --pwf-test 1500 \
  --pwf-min 0 --pwf-max 3000 --steps 25 --csv ipr_vogel.csv
```

All pressures are psia (consistent units) and rates are field units (e.g., STB/D). Keep units consistent.