# Cosmic Strings

Utilities for computing the cosmic-string angular power spectrum.

## Quick Start

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

If you want to run the notebook as well:

```bash
pip install -r requirements.txt
```

## Basic Usage

```python
from strings import compute_cl

ell = range(2, 2001, 4)
cl = compute_cl(
    5.0e-5,
    1.0,
    ell_arr=ell,
    n_k=48,
    n_chi=1024,
)
```

`compute_cl` uses the compiled estimator. On the first call it builds
`strings/libcosmic_strings_cl.dylib` with `gfortran` if needed.

You can also build it manually:

```bash
python -m strings.build_backend
```

## Notebook

The repository keeps one notebook:

- `notebooks/power_spectrum.ipynb`

It computes theory spectra with the `strings` backend and plots them for a few
parameter choices.

## Numerical Controls

- `n_k`: Gauss-Legendre nodes in `log(k)`.
- `n_chi`: Gauss-Legendre nodes in comoving distance.
- `k_min`, `k_max`: integration bounds for `k`.

