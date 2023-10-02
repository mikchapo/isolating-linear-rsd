# Compare pairwise velocity distributions in different mass bins for several TestBoxes
# v0.1.2, 2021-08-04 - Correct static solution

# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn
import sys


def H_z(z, H0, omegam0):
    return H0 * np.sqrt(omegam0 * (1. + z) * (1. + z) * (1. + z) + (1. - omegam0))


def omegam_calc(z, omegam0=0.315):
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def D_integrand(a, omegam, omegal):
    return ((omegam / (a**3.) / (omegam / (a**3.) + omegal))**0.55 - 1.) / a


def calc_D(z, omegam0, omegal):
    a = 1. / (1. + z)
    integral = quad(D_integrand, 0., a, args=(omegam0, omegal))
    return a*np.exp(integral[0])


def calc_v12_lin(ds, H0, ombh2, omch2, ns, sigma8, z):
    # Assumes ds is a numpy.array
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(ns=ns)
    pars.set_matter_power(redshifts=[z], kmax=2.0)

    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    kh, zcamb, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    sigma8_default = np.array(results.get_sigma8())
    omegam0 = (ombh2 + omch2) / (H0 / 100.)**2.
    simga8 = sigma8 * calc_D(z, omegam0, 1. - omegam0) / calc_D(0., omegam0, 1. - omegam0)
    pk = pk * (sigma8 / sigma8_default[0])**2.

    v12s = np.zeros(ds.shape)
    for i in range(ds.size):
        for j in range(kh.size-1):
            dk = kh[j+1] - kh[j]
            # log_dk = np.log10(kh[j+1]) - np.log10(kh[j])
            # k = np.power(10., (np.log10(kh[j]) + log_dk / 2.))
            v12s[i] += dk * kh[j] * pk[0][j] * spherical_jn(1, ds[i]*kh[j])
    v12s = -v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
    H = H_z(z, H0, omegam0)
    return H * v12s


test_cosmos = np.loadtxt("../data/cosmology_camb_test_box_full.dat")

z = 0.7
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
mass_bin_lims = np.arange(11.5, 14., 0.5)

H_seps = np.logspace(np.log10(0.01), np.log10(100.), 200)
v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)

for i in range(mass_bin_lims.size-1):
    plt.figure(num=2*i, figsize=(8., 6.))
    plt.axvline(x=7., linestyle="--", color="k")
    plt.axhline(y=0., linestyle="-", color="k")
    ymax = 0.
    ymin = 0.
    for j in range(6):
        tot_mvps = np.zeros(N_bins)
        tot_pairs = np.zeros(N_bins)
        for k in range(5):
            mvps = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox00{}-00{}/m200b_{:.2f}-{:.2f}_vp.dat".format(j, k, mass_bin_lims[i], mass_bin_lims[i+1]))
            tot_mvps = tot_mvps + np.nan_to_num(mvps[:, 1]) * mvps[:, 2]
            tot_pairs = tot_pairs + mvps[:, 2]
        tot_mvps = np.divide(tot_mvps, tot_pairs, out=np.zeros_like(tot_mvps), where=tot_pairs!=0)
        ymax = np.max((ymax, np.max(tot_mvps)))
        ymin = np.min((ymin, np.min(tot_mvps)))
        plt.plot(bin_centres, np.nan_to_num(tot_mvps), color="C{}".format(j), label="TestBox00{}".format(j), marker=".", linestyle="")

        H = H_z(z, test_cosmos[j, 3]*100., test_cosmos[j, 0])
        plt.plot(H_seps, -H*H_seps, color="C{}".format(j))

        ombh2 = test_cosmos[j, 1] * test_cosmos[j, 3]**2.
        omch2 = (test_cosmos[j, 0] - test_cosmos[j, 1]) * test_cosmos[j, 3]**2.
        v12s = calc_v12_lin(v12_seps, test_cosmos[j, 3]*100., ombh2, omch2, test_cosmos[j, 4], test_cosmos[j, 2], z)
        plt.plot(v12_seps, v12s, color="C{}".format(j), linestyle="--")

    plt.xlabel(r"s [$h^{-1}$ Mpc]")
    plt.ylabel(r"$v_p$ [km/s]")
    plt.xlim(0.01, 100.)
    plt.ylim(ymin-50., ymax+50.)
    plt.xscale("log")
    plt.legend()
    plt.title("{:.2f}-{:.2f} Halo Mass".format(mass_bin_lims[i], mass_bin_lims[i+1]))
    plt.savefig("../output/plots/m200b_{:.2f}-{:.2f}_vp_mass_comp_plot.jpg".format(mass_bin_lims[i], mass_bin_lims[i+1]))


    plt.figure(num=2*i+1, figsize=(8., 6.))
    plt.axvline(x=7., linestyle="--", color="k")
    plt.axhline(y=0., linestyle="-", color="k")
    
    test2_tot_mvps = np.zeros(N_bins)
    test2_tot_pairs = np.zeros(N_bins)
    for k in range(5):
        test2_mvps = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox002-00{}/m200b_{:.2f}-{:.2f}_vp.dat".format(k, mass_bin_lims[i], mass_bin_lims[i+1]))
        test2_tot_mvps = test2_tot_mvps + np.nan_to_num(test2_mvps[:, 1]) * test2_mvps[:, 2]
        test2_tot_pairs = test2_tot_pairs + test2_mvps[:, 2]
    test2_tot_mvps = np.divide(test2_tot_mvps, test2_tot_pairs, out=np.zeros_like(test2_tot_mvps), where=test2_tot_pairs!=0)
    plt.plot(bin_centres, np.nan_to_num(test2_tot_mvps), color="C{}".format(2))

    for j in range(6):
        if j!=2:
            tot_mvps = np.zeros(N_bins)
            tot_pairs = np.zeros(N_bins)
            for k in range(5):
                mvps = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_TestBox00{}-00{}/m200b_{:.2f}-{:.2f}_vp.dat".format(j, k, mass_bin_lims[i], mass_bin_lims[i+1]))
                tot_mvps = tot_mvps + np.nan_to_num(mvps[:, 1]) * mvps[:, 2]
                tot_pairs = tot_pairs + mvps[:, 2]
            tot_mvps = np.divide(tot_mvps, tot_pairs, out=np.zeros_like(tot_mvps), where=tot_pairs!=0)
            norm = np.mean(test2_tot_mvps[-14:] / tot_mvps[-14:])
            plt.plot(bin_centres, np.nan_to_num(norm*tot_mvps), color="C{}".format(j), label="Normalized TestBox00{}".format(j))

    plt.xlabel(r"s [$h^{-1}$ Mpc]")
    plt.ylabel(r"$v_p$ [km/s]")
    plt.xlim(0.01, 100.)
    plt.xscale("log")
    plt.legend()
    plt.title("{:.2f}-{:.2f} Halo Mass, normalized to TestBox2 using 20-100 Mpc/h".format(mass_bin_lims[i], mass_bin_lims[i+1]))
    plt.savefig("../output/plots/m200b_{:.2f}-{:.2f}_vp_mass_comp_normed_plot.jpg".format(mass_bin_lims[i], mass_bin_lims[i+1]))


# Change Log
# v0.1.1, 2021-07-08 - Added static and linear solutions
# v0.1.0, 2021-07-02 - Code Started