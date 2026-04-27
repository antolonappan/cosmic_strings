from importlib import resources
from pathlib import Path

import numpy as np


def _candidate_paths(filename):
    packaged = resources.files("strings") / "_data" / filename
    yield Path(packaged)

    repo_root = Path(__file__).resolve().parent.parent
    yield repo_root / "data" / filename


class _PowerSpectrum:
    def __init__(self, path, fac_func):
        data = np.loadtxt(path, comments="#")
        self.L  = data[:, 0]
        fac     = fac_func(self.L)
        self.CL = data[:, 1] * fac
        self.ER = data[:, 2] * fac


def _load_spectrum(filename, fac_func):
    for path in _candidate_paths(filename):
        if path.exists():
            return _PowerSpectrum(path, fac_func)

    raise FileNotFoundError(f"Could not find data file: {filename}")


ACTDR6  = _load_spectrum(
    "act_mv.txt",
    fac_func=lambda L: (np.sqrt(2) * np.pi) / (L**2 * (L + 1)**2),
)
PLANCK18 = _load_spectrum(
    "planck_TT.txt",
    fac_func=lambda L: (2 * np.pi) / (L**2 * (L + 1)**2),
)


__all__ = ["ACTDR6", "PLANCK18"]
