# Cosmic Strings

Small utilities for computing a cosmic-string angular power spectrum and loading bundled ACT DR6 and Planck 2018 data.

The repository is usable in two ways:

1. Directly from a cloned checkout.
2. As a lightweight installable package with bundled data files.

## What is included

- `strings.compute_cl(...)` for computing the model spectrum.
- `strings.ACTDR6` with attributes `.L`, `.CL`, `.ER`.
- `strings.PLANCK18` with attributes `.L`, `.CL`, `.ER`.
- Source notebooks under `notebooks/`.

The ACT and Planck spectra are shipped with the code, so imports do not depend on the repository layout once installed.

## Quick Start

Clone the repository and install it in editable mode:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

If you also want to run the notebooks:

```bash
pip install -r requirements.txt
```

## Basic Usage

```python
from strings import ACTDR6, PLANCK18, compute_cl

print(ACTDR6.L[:3])
print(ACTDR6.CL[:3])
print(ACTDR6.ER[:3])

cl_act = compute_cl(8.0e-4, 1.0, ell_arr=ACTDR6.L)
cl_planck = compute_cl(8.0e-4, 1.0, ell_arr=PLANCK18.L)
```

`ACTDR6.CL` and `ACTDR6.ER` are pre-scaled by

$$
\frac{\sqrt{2}\pi}{L^2(L+1)^2}
$$

and `PLANCK18.CL` and `PLANCK18.ER` are pre-scaled by

$$
\frac{2\pi}{L^2(L+1)^2}.
$$

## Project Layout

- `strings/`: importable code.
- `strings/_data/`: packaged data files used at runtime.
- `data/`: original data files kept in the repository.
- `notebooks/`: exploratory and analysis notebooks.

## Notes

- `compute_cl` is numerically expensive for large multipole ranges.
- The current implementation uses a thread pool internally and returns an array of `C_l` values for the requested `ell_arr` or `lmax` range.
- When `ell_arr` is passed, the returned array is aligned with that input ordering.

## Sharing This Repository

For another user, the simplest setup is:

```bash
git clone <repo-url>
cd cosmic_strings
python -m venv .venv
source .venv/bin/activate
pip install -e .
pip install -r requirements.txt
```

Then they can either import `strings` from Python or open the notebooks.