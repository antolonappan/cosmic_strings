# Fortran CL Backend

This folder is a side-by-side migration path for the expensive cosmic-string
`C_L` calculation. The likelihoods, plotting, notebooks, `emcee`, and GetDist
can stay in Python.

## Notebook Usage

From a notebook in `notebooks/`:

```python
import sys
from pathlib import Path

sys.path.append(str(Path("..") / "fortran_backend"))
import cosmic_cl_backend as fcl

compute_cl_backend = fcl.compute_cl
```

Then replace calls like:

```python
cs.compute_cl(Gmu_ref, P_fixed, ell_arr=ell, N_k=30)
```

with:

```python
compute_cl_backend(
    Gmu_ref,
    P_fixed,
    ell_arr=ell,
    n_k=48,
    n_chi=1024,
)
```

The first call builds the Fortran shared library with `gfortran`. You can also
build it manually:

```bash
cd fortran_backend
python build_backend.py
```

## Numerical Knobs

- `n_k`: Gauss-Legendre nodes in `log(k)`.
- `n_chi`: Gauss-Legendre nodes in comoving distance.
- `k_min`, `k_max`: integration bounds for `k`.

The old Python code used `scipy.quad` for the `chi` integral and trapezoid
integration over log-spaced `k` samples. The Fortran backend instead performs
Gauss-Legendre integration in both `chi` and `log(k)`, which is meant to make
the theory amplitude easier to converge and diagnose. It is loaded from Python
with `ctypes`, so it does not need Meson or f2py.

If the curve is still offset from Namikawa/Toshiya, scan `n_k`, `n_chi`,
`k_min`, and `k_max` before changing the data scaling.
