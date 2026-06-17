from __future__ import annotations

import ctypes
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from numpy.ctypeslib import ndpointer


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "cosmic_strings_cl.f90"
LIB_NAME = "libcosmic_strings_cl.dylib" if sys.platform == "darwin" else "libcosmic_strings_cl.so"
LIB_PATH = ROOT / LIB_NAME


def build_extension(force: bool = False, verbose: bool = True) -> Path:
    """Build the Fortran shared library in this package."""

    if LIB_PATH.exists() and not force:
        return LIB_PATH

    compiler = os.environ.get("FC") or shutil.which("gfortran") or "/opt/homebrew/bin/gfortran"
    if not Path(compiler).exists() and shutil.which(compiler) is None:
        raise RuntimeError(
            "Could not find gfortran. Install gfortran or set FC=/path/to/gfortran."
        )

    command = [
        compiler,
        "-O3",
        "-shared",
        "-fPIC",
        "-ffree-line-length-none",
        "-o",
        str(LIB_PATH),
        str(SRC),
    ]

    stdout = None if verbose else subprocess.DEVNULL
    stderr = None if verbose else subprocess.DEVNULL
    subprocess.run(command, cwd=ROOT, check=True, stdout=stdout, stderr=stderr)

    if not LIB_PATH.exists():
        raise RuntimeError("gfortran completed, but no shared library was found")
    return LIB_PATH


def _load_library():
    if not LIB_PATH.exists():
        build_extension(force=False, verbose=True)

    lib = ctypes.CDLL(str(LIB_PATH))
    lib.compute_cl_array_c.argtypes = [
        ctypes.c_double,
        ctypes.c_double,
        ndpointer(dtype=np.int32, ndim=1, flags="C_CONTIGUOUS"),
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS"),
    ]
    lib.compute_cl_array_c.restype = None

    return lib


def compute_cl(
    G_mu,
    P,
    *,
    lmax=2000,
    ell_arr=None,
    n_k=48,
    n_chi=1024,
    k_min=1.0e-4,
    k_max=1.0e-1,
    progress=False,
    max_workers=None,
    N_k=None,
):
    """Compute C_L with the Fortran backend."""

    if N_k is not None:
        n_k = N_k

    if ell_arr is None:
        ell = np.arange(2, lmax + 1, dtype=np.int32)
    else:
        ell = np.asarray(np.floor(np.asarray(ell_arr, dtype=float) + 0.5), dtype=np.int32)

    lib = _load_library()
    cl = np.empty(len(ell), dtype=np.float64)
    lib.compute_cl_array_c(
        float(G_mu),
        float(P),
        ell,
        int(len(ell)),
        int(n_k),
        int(n_chi),
        float(k_min),
        float(k_max),
        cl,
    )
    return cl


def compute_cl_ell(
    ell,
    G_mu,
    P,
    *,
    n_k=48,
    n_chi=1024,
    k_min=1.0e-4,
    k_max=1.0e-1,
):
    """Compute one multipole with the Fortran backend."""

    return float(
        compute_cl(
            G_mu,
            P,
            ell_arr=np.array([ell]),
            n_k=n_k,
            n_chi=n_chi,
            k_min=k_min,
            k_max=k_max,
        )[0]
    )


__all__ = ["build_extension", "compute_cl", "compute_cl_ell"]