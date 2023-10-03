# Functions needed for the fit
# v0.1.1, 2021-04-07 - Copied from aemulus_eboss_v2

# Imports
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

import numpy as np


def invert_cov_matrix(cov_matrix, n):
    u, s, vh = np.linalg.svd(cov_matrix)
    inv_cov_matrix = np.asmatrix(vh).H * np.diag(1. / s) * np.asmatrix(u).H

    # Apply Hartlap factor (Missing before 2020-10-23)
    p = inv_cov_matrix.shape[0]
    inv_cov_matrix = inv_cov_matrix * float((n - p - 2.) / (n - 1.))

    return inv_cov_matrix


c = 299792.458
Omega_M_fid = 0.31
z_aemulus = 0.7


def E(z, Omega_M):
    return np.sqrt(Omega_M*np.power(1. + z, 3.) + (1. - Omega_M))


def D_M_integrand(z, Omega_M):
    return 1. / E(z, Omega_M)


def D_M(Omega_M, z):
    dH = c / 100.
    return dH * quad(D_M_integrand, 0, z, args=Omega_M)[0]


def D_H(Omega_M, z):
    return 1. / 100. / E(z, Omega_M)


def D_integrand(a, omegam, omegal):
    return ((omegam / (a**3.) / (omegam / (a**3.) + omegal))**0.55 - 1.) / a


def calc_D(z, omegam, omegal):
    a = 1. / (1. + z)
    integral = quad(D_integrand, 0., a, args=(omegam, omegal))
    return a*np.exp(integral[0])


def omegam(z, omegam0=0.315):
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
    l = 0.83 + 0.19 / (omegam0**0.8)
    beta = 0.98 + 0.01 / (omegam0**1.2)
    gamma = 0.64 + 0.09 / (omegam0**0.4)
    return l * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))


def ap_scaling(mono_pred, quad_pred, wp_pred, alpha_perp, alpha_para):
    alpha = alpha_para**(1./3.) * alpha_perp**(2./3.)
    epsilon = (alpha_para/alpha_perp)**(1./3.) - 1.

    # Cubic Spline in Display Space
    mono_pred_scaled = CubicSpline(seps, mono_pred*np.power((seps), 2.))
    quad_pred_scaled = CubicSpline(seps, quad_pred*np.power((seps), 2.))
    wp_pred_scaled = CubicSpline(np.log10(seps), np.log10(wp_pred))

    mono_pred_fid = mono_pred_scaled(alpha*seps) / np.power(alpha*seps, 2.)  + 0.4 * epsilon * (quad_pred_scaled(alpha*seps) / np.power(alpha*seps, 2.) + quad_pred_scaled(alpha*seps, 1) / (alpha*seps))
    quad_pred_fid = -4. * epsilon * mono_pred_scaled(alpha*seps) / np.power(alpha*seps, 2.) + (1. - 2. / 7. * epsilon) * quad_pred_scaled(alpha*seps) / np.power(alpha*seps, 2.) + 2. * epsilon / (alpha*seps) * mono_pred_scaled(alpha*seps, 1) + 4. / 7. * epsilon / (alpha*seps) * quad_pred_scaled(alpha*seps, 1)
    wp_pred_fid = np.power(10., wp_pred_scaled(np.log10(alpha_perp*seps)))

    return mono_pred_fid, quad_pred_fid, wp_pred_fid


def find_indices(input_dict):
    indices = []
    edges = [0, 9, 18, 27]
    measurements = ["mono", "quad", "wp"]
    scales = ["small", "intermediate", "large"]
    for i in range(len(measurements)):
        if measurements[i] in input_dict["measurements"]:
            indices += list(range(edges[i], edges[i+1]))

    remove_indices = []
    for i in range(len(scales)):
        if scales[i] not in input_dict["scales"]:
            for j in range(3):
                remove_indices += list(range(edges[j]+3*i, edges[j]+3*(i+1)))

    for remove_index in remove_indices:
        if remove_index in indices: indices.remove(remove_index)

    print("Indices:", indices)
    return indices

