from .cl import compute_cl
from .data import ACTDR6, PLANCK13, PLANCK18
import numpy as np


def DR6_FACTOR(ell):
    return np.sqrt(ell) * (ell*(ell+1))**2 / (2 * np.pi)


__all__ = [
    "compute_cl",
    "ACTDR6",
    "PLANCK18",
    "PLANCK13",
    "DR6_FACTOR",
]
