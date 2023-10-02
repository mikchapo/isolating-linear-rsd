"""
Compare pairwise velocity distributions from particle catalogues.

v1.0.1, 2021-11-08 - Updated for PEP8 compliance
"""

# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn
import sys


def H_z(z, H0, omegam0):
    """
    Calculate the Hubble constant as a function of redshift.

    Input:
    ------
    z - Redshift at which to evaluate the Hubble constant
    H0 - Hubble constant at redshift 0
    omegam0 - Omegam matter at redshift 0

    Output:
    -------
    H - Hubble constant at redshift z
    """
    H = (H0 * np.sqrt(omegam0 * (1. + z) * (1. + z) * (1. + z) +
         (1. - omegam0)))
    return H


def interpolator_intergrator(k, z, d, interpolator, k_hunit, h):
    """
    Call a CAMB matter power interpolator for integration.

    Input:
    ------
    k - Wavenumber at which to evaluate the power spectrum
    z - Redshift at which to evaluate the power spectrum
    d - Separation at which to calculate the pairwise velocity
    interpolator - CAMB matter power interpolator instance that is being
                   integrated
    k_hunit - Boolean to determine whether k has units of [h/Mpc] (True) or
              [Mpc] (False)
    h - Dimensionless Hubble constant, needed for converting distance units

    Ouput:
    ------
    The integrand for the power specturm intrgral in the pairwise linear
    velocity calculation.
    """
    if k_hunit:
        return k * interpolator.P(z, k) * spherical_jn(1, d*k)
    else:
        # If using k_hunit then need to convert d from units of [Mpc/h]
        return k * interpolator.P(z, k) * spherical_jn(1, d*k/h)


def vp_lin_pred(ds, z, H0, omch2, ombh2, ns, sigma8, kmax=2.0, minkh=1e-4,
                maxkh=1, k_hunit=False, hubble_units=False):
    """
    Calculate the pairwise velocity linear theory prediction.

    See Eq.6 of Belloso et al. 2012 (https://arxiv.org/pdf/1204.5761.pdf) for
    the theoretical source. It should be noted that that equation gives a
    displacement, which needs to be scaled by aH(a) to give a velocity. The
    bias factor, b, is assumed to be 1.

    Input:
    ------
    ds - Separations at which to output velocity predictions as vp(s). Assumed
         to be a numpy.array
    z - Redshift at which to evaluate the pairwise velocity
    H0 - Hubble constant at z=0
    omch2 - Omega cold dark matter multiplied by h**2
    ombh2 - Omega baryon multiplied by h**2
    ns - The initial power spectrum slope
    sigma8 - RMS fluctuations in a sphere of 8 Mpc/h
    kmax - Maximum k scale for initializing the CAMB power spectrum.
           Insignificant effect on results
    minkh - Lower limit of the power spectrum integral. Theoretically should
            be 0. Minimally affects the results
    maxkh - Upper limit of the power spectrum integral. Theoretically should
            be infinity. Doesn't alter large separation slope, but increasing
            it moves the up turn to lower separations in a significant way
    k_hunit - Boolean to determine whether to use [h/Mpc] (True) or [Mpc]
              units for k
    hubble_units - Boolean to determine whether to use [(Mpc/h)^3] (True) or
                   [Mpc^3] units for P(k)

    Output:
    -------
    vp_lin - Numpy array of linear pairwise velocity predictions

    """
    # Default CAMB power spectrum normalization
    As_default = 2e-09

    # List of redshifts required for CAMB calulcations without duplicates, to
    # avoid CAMB error
    redshifts = list(dict.fromkeys([0., z]))

    # Initialize CAMB parameters for power spectrum calculation
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(ns=ns)
    pars.set_matter_power(redshifts=redshifts, kmax=kmax)
    pars.NonLinear = camb.model.NonLinear_none

    # Get the value of sigma8(z=0) with the default normalization in order to
    # calculate a new normalization
    results = camb.get_results(pars)
    sigma8_0_default = np.array(results.get_sigma8_0())
    As_rescale = (sigma8 / sigma8_0_default) ** 2.

    # Re-initalize parameters and get results with new normalization
    pars.InitPower.set_params(ns=ns, As=As_default*As_rescale)
    pars.set_matter_power(redshifts=redshifts, kmax=kmax)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)

    # Generate a spline-fit interpolator for the matter power spectrum to use
    # with a numerical integrator
    interp = camb.get_matter_power_interpolator(pars, zs=[z], nonlinear=False,
                                                hubble_units=hubble_units,
                                                k_hunit=k_hunit)
    vp_lin = np.zeros(ds.shape)
    # Loops through each separation, calculating v_pair(s)
    for i in range(ds.size):
        if k_hunit:
            vp_lin[i] = quad(interpolator_intergrator, minkh, maxkh, args=(z,
                             ds[i], interp, k_hunit, H0/100.))[0]
        else:
            # If not using k_hunit, rescales minkh and maxkh to give
            # consistent integration limits
            vp_lin[i] = quad(interpolator_intergrator, minkh*H0/100.,
                             maxkh*H0/100., args=(z, ds[i], interp, k_hunit,
                                                  H0/100.))[0]

    # Corrects any mismatched factors of h from the integral of the power
    # spectrum so that the integral gives units of Mpc. When multiplied by H
    # [km/s/Mpc] at the end this will give a velocity in km/s
    h_correction = (H0 / 100.)**(2 * k_hunit - 3 * hubble_units)

    fsigma8 = np.array(results.get_fsigma8())[0]
    sigma8_default = np.array(results.get_sigma8())[0]
    H = results.hubble_parameter(z)
    # Adds the missing constant factors
    vp_lin = (-H / (z + 1.) * vp_lin * fsigma8 / sigma8_default / np.pi**2. *
              h_correction)
    return vp_lin


def initialize_vp_plot():
    """Initialize a figure for plotting pairwise velocities."""
    plt.figure(figsize=(8., 6.), dpi=300)
    plt.axhline(y=0., linestyle="-", color="k")


def initialize_paper_plot(thesis=False):
    """Initialize a figure for plotting pairwise velocities."""

    if thesis:
        width = 6.375
        length = 7.

    else:
        width = 9.5
        length = 11.

    fig, axes = plt.subplots(3, 2, figsize=(width, length), dpi=300, sharex=True,
                             sharey=True)

    for a in axes:
        for ax in a:
            ax.axhline(y=0., linestyle="-", color="k")
            ax.set_xlim(0.01, 100.)
            ax.set_ylim(-1200., 400.)
            ax.set_xscale("log")

    axes[2, 0].set_xlabel(r"s [$h^{-1}$ Mpc]")
    axes[2, 1].set_xlabel(r"s [$h^{-1}$ Mpc]")
    axes[0, 0].set_ylabel(r"$v_p$ [km/s]")
    axes[1, 0].set_ylabel(r"$v_p$ [km/s]")
    axes[2, 0].set_ylabel(r"$v_p$ [km/s]")
    return fig, axes


def finish_vp_plot(output_path, legend=True):
    """
    Finish and save a pairwise velocity matplotlib.pyplot plot.

    Input:
    ------
    output_path - Path for saving the plot
    """
    plt.xlabel(r"s [$h^{-1}$ Mpc]")
    plt.ylabel(r"$v_p$ [km/s]")
    plt.xlim(0.01, 100.)
    plt.xscale("log")
    if legend:
        plt.legend()
    plt.savefig(output_path)


def finish_paper_plot(output_path, legend=True):
    """
    Finish and save a pairwise velocity matplotlib.pyplot plot.

    Input:
    ------
    output_path - Path for saving the plot
    """
    plt.tight_layout()
    plt.savefig(output_path)


def log_bin_centres(ll, ul, N):
    """
    Calculcate the edges and centres of log spaced separation bins.

    Input:
    ------
    ll - Lower limit of bins
    ul - Upper limit of bins
    N - Number of bins

    Output:
    -------
    bin_centres - Centres of the logarithmically spaced bins
    """
    bin_edges = np.logspace(np.log10(ll), np.log10(ul), N+1)
    log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
    bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
    return bin_centres


def plot_vp(vp_path, label, color="C0", N_bins=80, marker=".", linestyle="-"):
    """
    Plot pairwise velocities stored in a data file with errorbars.

    Input:
    ------
    vp_path - Path to the pairwise velocity data file. The columns of the file
              should be ordered:
              s  vp  var(vp) N_pairs
    label - Label for the plot legend
    color - Color for plot
    N_bins - Number of separation bins used, with limits assumed to be 0.01 -
             100 Mpc/h
    """
    bin_centres = log_bin_centres(0.01, 100., N_bins)

    vps = np.loadtxt(vp_path)
    vps[:, 1] = np.nan_to_num(vps[:, 1])
    # Errorbars equal to sqrt(var/N), may be underestimated
    plt.errorbar(bin_centres, vps[:, 1],
                 yerr=np.sqrt(vps[:, 2] / vps[:, 3]), color=color,
                 label=label, marker=marker, linestyle=linestyle)


def plot_vp_paper(vp_path, axes, m, n, color="C0", label=None, N_bins=80,
                  marker=".", linestyle="-"):
    """
    Plot pairwise velocities stored in a data file with errorbars.

    Input:
    ------
    vp_path - Path to the pairwise velocity data file. The columns of the file
              should be ordered:
              s  vp  var(vp) N_pairs
    label - Label for the plot legend
    color - Color for plot
    N_bins - Number of separation bins used, with limits assumed to be 0.01 -
             100 Mpc/h
    """
    bin_centres = log_bin_centres(0.01, 100., N_bins)

    vps = np.loadtxt(vp_path)
    vps[:, 1] = np.nan_to_num(vps[:, 1])
    # Errorbars equal to sqrt(var/N), may be underestimated
    axes[m, n].errorbar(bin_centres, vps[:, 1],
                        yerr=np.sqrt(vps[:, 2] / vps[:, 3]), color=color,
                        label=label, marker=marker, linestyle=linestyle)


def plot_static(z, cosmo_params):
    """
    Calculate and plot the static solution.

    The static solution is the pairwise velocity required to
    maintain a constant proper distance in an expanding background.

    Input:
    ------
    z - Redshift at which to evaluate the static solution
    cosmo_params - Parameters describing the desired cosmology, ordered as:
                   H0   omch2   ombh2   ns  sigma8
    """
    H_seps = np.logspace(np.log10(0.01), np.log10(10.), 200)
    H = H_z(z, cosmo_params[0], (cosmo_params[1] + cosmo_params[2]) /
            (cosmo_params[0] / 100.)**2.)
    plt.plot(H_seps, -H*H_seps / (cosmo_params[0] / 100.), color="k",
             label="Static solution")


def plot_static_paper(z, cosmo_params, axes):
    """
    Calculate and plot the static solution.

    The static solution is the pairwise velocity required to
    maintain a constant proper distance in an expanding background.

    Input:
    ------
    z - Redshift at which to evaluate the static solution
    cosmo_params - Parameters describing the desired cosmology, ordered as:
                   H0   omch2   ombh2   ns  sigma8
    """
    H_seps = np.logspace(np.log10(0.01), np.log10(10.), 200)
    H = H_z(z, cosmo_params[0], (cosmo_params[1] + cosmo_params[2]) /
            (cosmo_params[0] / 100.)**2.)
    include_label = True
    for a in axes:
        for ax in a:
            if include_label:
                ax.plot(H_seps, -H*H_seps / (cosmo_params[0] / 100.), color="k",
                        label="Static solution")
                include_label = False

            else:
                ax.plot(H_seps, -H*H_seps / (cosmo_params[0] / 100.), color="k")


def plot_lin_pred(z, cosmo_params):
    """
    Calculate and plot the linear theory prediction.

    Input:
    ------
    z - Redshift at which to evaluate the linear theory prediction
    cosmo_params - Parameters describing the desired cosmology, ordered as:
                   H0   omch2   ombh2   ns  sigma8
    """
    vp_lin_seps = np.logspace(np.log10(1.), np.log10(100.), 40)
    vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1],
                         cosmo_params[2], cosmo_params[3], cosmo_params[4])
    plt.plot(vp_lin_seps, vp_lin, color="k", linestyle="--",
             label="Linear theory prediction")


def plot_lin_pred_paper(z, cosmo_params, axes):
    """
    Calculate and plot the linear theory prediction.

    Input:
    ------
    z - Redshift at which to evaluate the linear theory prediction
    cosmo_params - Parameters describing the desired cosmology, ordered as:
                   H0   omch2   ombh2   ns  sigma8
    """
    vp_lin_seps = np.logspace(np.log10(1.), np.log10(100.), 40)
    vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1],
                         cosmo_params[2], cosmo_params[3], cosmo_params[4])
    include_label = True
    for a in axes:
        for ax in a:
            if include_label:
                ax.plot(vp_lin_seps, vp_lin, color="k", linestyle="--",
                        label="Linear theory prediction")
                include_label = False

            else:
                ax.plot(vp_lin_seps, vp_lin, color="k", linestyle="--")

def pairwise_vel_plot(vp_paths, z, cosmo_params, output_path, static=True,
                      lin_pred=True):
    """
    Call all the components of a pairwise velocity plot.

    Input:
    ------
    vp_paths - A list of paths from which to load pairwise velocity data files
    z - Redshift of evaluation
    cosmo_params - Parameters describing the desired cosmology, ordered as:
                   H0   omch2   ombh2   ns  sigma8
    output_path - Path for saving the plot
    static - Boolean to determine whether to include the static solution
    lin_pred - Boolean to determine whether to include the linear theory
               prediction
    """
    initialize_vp_plot()
    for i, vp_path in enumerate(vp_paths):
        # As a basic choice for the label, get the file name from the path and
        # strip the file type
        label = vp_path.split('/')[-1][:-4]
        plot_vp(vp_path, label, color="C{}".format(i))
    if static:
        plot_static(z, cosmo_params)
    if lin_pred:
        plot_lin_pred(z, cosmo_params)
    finish_vp_plot(output_path)


# Change Log
# v1.0.0, 2021-10-18 - After figuring out the proper form, archived the
#                      original code as _v1 and cleaned up this version for
#                      use
# v0.1.0 - v1.0.0 - Numerous rounds of debugging, adding new features, and
#                   correcting factors
# v0.1.0, 2021-09-15 - Code started from plot_vp_mass_comp_abacus.py
