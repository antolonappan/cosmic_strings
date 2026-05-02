# Cosmic Strings

Small utilities for computing a cosmic-string angular power spectrum and loading bundled ACT DR6, Planck 2018, and Planck 2013 data.

The repository is usable in two ways:

1. Directly from a cloned checkout.
2. As a lightweight installable package with bundled data files.

## What is included

- `strings.compute_cl(...)` for computing the model spectrum.
- `strings.ACTDR6` with attributes `.L`, `.CL`, `.ER`.
- `strings.PLANCK18` with attributes `.L`, `.CL`, `.ER`.
- `strings.PLANCK13` with attributes `.L`, `.CL`, `.ER`.
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
from strings import ACTDR6, PLANCK13, PLANCK18, compute_cl

print(ACTDR6.L[:3])
print(ACTDR6.CL[:3])
print(ACTDR6.ER[:3])
print(PLANCK13.L[:3])

cl_act = compute_cl(8.0e-4, 1.0, ell_arr=ACTDR6.L)
cl_planck = compute_cl(8.0e-4, 1.0, ell_arr=PLANCK18.L)
cl_planck13 = compute_cl(8.0e-4, 1.0, ell_arr=PLANCK13.L)
```

## Fortran Backend

An experimental Fortran backend lives in `fortran_backend/`. It keeps plotting,
likelihoods, `emcee`, and GetDist in Python, but moves the expensive `C_L`
calculation into compiled Fortran.

From a notebook in `notebooks/`, use:

```python
import sys
from pathlib import Path

PROJECT_ROOT = Path.cwd()
if not (PROJECT_ROOT / "strings").exists():
    PROJECT_ROOT = PROJECT_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT / "fortran_backend"))
import cosmic_cl_backend as fcl

compute_cl_backend = fcl.compute_cl
cl = compute_cl_backend(
    1.0e-4,
    1.0,
    ell_arr=ACTDR6.L,
    n_k=48,
    n_chi=1024,
)
```

The first call builds `fortran_backend/libcosmic_strings_cl.dylib` with
`gfortran` if needed. You can also build it manually:

```bash
cd fortran_backend
/opt/anaconda3/bin/python3 build_backend.py
```

Use `n_k`, `n_chi`, `k_min`, and `k_max` to test convergence of the theory
amplitude before changing any data scaling.

`ACTDR6.CL` and `ACTDR6.ER` are pre-scaled by

$$
\frac{\sqrt{2}\pi}{L^2(L+1)^2}
$$

and `PLANCK18.CL` and `PLANCK18.ER` are pre-scaled by

$$
\frac{2\pi}{L^2(L+1)^2}.
$$

`PLANCK13.CL` and `PLANCK13.ER` are pre-scaled by

$$
10^{-8}\frac{2\pi}{L^2(L+1)^2}.
$$

For reproducing the Planck 2013 string limit, the notebook keeps the likelihood
in these loaded `CL` arrays, but bins the theory before comparison. The
Planck13 table is a binned bandpower table, so `notebooks/G_mu.ipynb` uses
center-inferred integer bins, computes the theory over the multipoles in each
bin, averages
`1e8 L^2(L+1)^2 C_L/(2*pi)` with `(2L+1)` weights, and converts the result back
to the same center-`CL` convention used by `PLANCK13.CL`.

For ACT and Planck18 the current data files only contain bin centers, not bin
edges. The notebook floors those centers to integer multipoles, matching the
old scipy notebook convention.

What is still missing for an exact reproduction is the full experiment
bandpower/window function. The center-inferred top-hat binning is kept because
it stays closest to Namikawa/Toshiya with the current data table and theory
curve.

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
