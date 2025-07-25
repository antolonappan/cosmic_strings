import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, trapezoid
from scipy.special import spherical_jn,  gammaln
from tqdm import tqdm
from astropy.cosmology import Planck18
from numba import njit
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
import os


c_tilde = 0.23
z_rec = 1089
chi_star = Planck18.comoving_distance(z_rec).value


z_vals = np.linspace(0.01, z_rec, 1000)
chi_of_z_vals = Planck18.comoving_distance(z_vals).value
h_vals = (Planck18.H(z_vals).value / 299792.458)


@njit(cache=True)
def interp1d_jit(x_arr, y_arr, xi):
    idx = np.searchsorted(x_arr, xi)
    if idx == 0:
        return y_arr[0]
    if idx >= x_arr.shape[0]:
        return y_arr[-1]
    x0 = x_arr[idx-1]; x1 = x_arr[idx]
    y0 = y_arr[idx-1]; y1 = y_arr[idx]
    t = (xi - x0) / (x1 - x0)
    return y0 + t * (y1 - y0)


@njit(cache=True)
def z_of_chi_jit(chi):
    return interp1d_jit(chi_of_z_vals, z_vals, chi)

@njit(cache=True)
def H_of_chi_jit(chi):
    return interp1d_jit(chi_of_z_vals, h_vals, chi)


@njit(cache=True)
def delta1sq(k, chi, G_mu, P):
    gamma = np.sqrt(np.pi * np.sqrt(2.0) / (3.0 * c_tilde * P))
    v_sq = 0.5 * (1.0 - np.pi / (3.0 * gamma))
    v = np.sqrt(v_sq)


    z = z_of_chi_jit(chi)
    a = 1.0 / (1.0 + z)
    H = H_of_chi_jit(chi)

    xi = 1.0 / (H * gamma)
    prefactor = (16.0 * np.pi * G_mu)**2 * (np.sqrt(6.0 * np.pi) * v**2 /
                                           (12.0 * (1.0 - v**2)))
    scale = (a / (k * xi))**5
    er_term = math.erf(k * xi / (a * 2.0 * np.sqrt(6.0)))

    return prefactor * (4.0 * np.pi * k**3 * chi**2 * a**4 / H) * scale * er_term


def compute_cl_ell(l_val, G_mu, P, N_k=800):
    def inner_integrand(chi, k):
        return np.sqrt(delta1sq(k, chi, G_mu, P)) \
               * spherical_jn(l_val, k * chi) / chi

    def cl_integrand(k):
        inner, _ = quad(inner_integrand,
                        0, chi_star,
                        args=(k,),
                        limit=200,
                        epsabs=1e-8,
                        epsrel=1e-8)
        return (inner**2) / k

    k_vals = np.logspace(-4, -1, N_k)
    cl_values = [cl_integrand(k) for k in k_vals]
    Cl_result = trapezoid(cl_values, k_vals)

    log_prefactor = np.log(4.0 * np.pi) + gammaln(l_val) - gammaln(l_val + 2)
    prefactor = np.exp(log_prefactor)
    return prefactor * Cl_result / (2.0 * np.pi)

def compute_cls(G_mu, P, N_k=30, l_max=500, progress_bar=True, custom_l_arr=None):
    if custom_l_arr is not None:
        l_vals = custom_l_arr
    else:
        l_vals = np.arange(2, l_max + 1)
    cl_results = []
    for l in tqdm(l_vals, desc="Computing C_l", unit="ℓ", disable=not progress_bar):
        cl_results.append(compute_cl_ell(l, G_mu, P, N_k=N_k))
    return l_vals, np.array(cl_results)

def compute_cls_mt(G_mu, P, N_k=30, l_max=500, nworkers=None):
    if nworkers is None:
        nworkers = os.cpu_count() or 1

    l_vals = np.arange(2, l_max + 1)
    cl_results = np.empty_like(l_vals, dtype=float)


    def task(l):
        return l, compute_cl_ell(l, G_mu, P, N_k=N_k)

    with ThreadPoolExecutor(max_workers=nworkers) as executor:

        futures = {executor.submit(task, l): idx for idx, l in enumerate(l_vals)}
        for future in tqdm(as_completed(futures),
                           total=len(futures),
                           desc="Computing C_l"):
            l, cl_val = future.result()
            idx = futures[future]
            cl_results[idx] = cl_val

    return l_vals, cl_results


    


    

