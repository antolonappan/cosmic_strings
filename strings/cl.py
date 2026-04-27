import numpy as np
from concurrent.futures import ThreadPoolExecutor
from scipy.integrate import quad, cumulative_trapezoid, trapezoid
from scipy.interpolate import interp1d
from scipy.special import spherical_jn, erf, gammaln
from tqdm import tqdm
import os


c_tilde = 0.23
H0 = 70.0
h = H0 / 299792.458
Omega_m = 0.3
Omega_L = 0.7
chi_star = 14000.0


z_vals = np.linspace(0.01, 1100, 1000)
H_z = h * np.sqrt(Omega_m * (1 + z_vals)**3 + Omega_L)

chi_of_z = cumulative_trapezoid(1 / H_z, z_vals, initial=0)
z_of_chi_func = interp1d(chi_of_z, z_vals, bounds_error=False,
                         fill_value=(z_vals[0], z_vals[-1]))

def a_of_chi(chi):
    return 1 / (1 + z_of_chi_func(chi))

def H_of_chi(chi):
    z = z_of_chi_func(chi)
    return h * np.sqrt(Omega_m * (1 + z)**3 + Omega_L)

def delta1sq(k, chi, G_mu, P):
    gamma = np.sqrt(np.pi * np.sqrt(2) / (3 * c_tilde * P))
    v_sq = 0.5 * (1 - np.pi / (3 * gamma))
    v = np.sqrt(v_sq)

    xi = 1 / (H_of_chi(chi) * gamma)
    a = a_of_chi(chi)
    H = H_of_chi(chi)

    prefactor = (16 * np.pi * G_mu)**2 * (np.sqrt(6 * np.pi) * v**2 / (12 * (1 - v**2)))
    scale = (a / (k * xi))**5
    er_term = erf(k * xi / (a * 2 * np.sqrt(6)))

    return prefactor * (4 * np.pi * k**3 * chi**2 * a**4 / H) * scale * er_term


def compute_cl_ell(l_val, G_mu, P, N_k=30):

    def inner_integrand(chi, k):
        return np.sqrt(delta1sq(k, chi, G_mu, P)) * spherical_jn(l_val, k * chi) / chi

    def cl_integrand(k):
        inner, _ = quad(inner_integrand, 0, chi_star, args=(k,), limit=200)
        return (inner)**2 / k

    k_vals = np.logspace(-4, -1, N_k)
    cl_vals = [cl_integrand(k) for k in k_vals]

    Cl = trapezoid(cl_vals, k_vals)

    log_pref = np.log(4 * np.pi) + gammaln(l_val) - gammaln(l_val + 2)
    prefactor = np.exp(log_pref)
    return prefactor * Cl 


def compute_cl(G_mu, P, lmax=2000, ell_arr=None, progress=False, max_workers=None):
    if ell_arr is None:
        ell = np.arange(2, lmax + 1)
    else:
        ell = ell_arr
    ell = np.asarray(ell)

    def _compute_single(l_val):
        return compute_cl_ell(l_val, G_mu, P)
    
    if max_workers is None:
        max_workers = os.cpu_count()

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        cl_iter = executor.map(_compute_single, ell)
        if progress:
            Cl = list(tqdm(cl_iter, total=len(ell), desc="Computing Cl"))
        else:
            Cl = list(cl_iter)

    return np.array(Cl)