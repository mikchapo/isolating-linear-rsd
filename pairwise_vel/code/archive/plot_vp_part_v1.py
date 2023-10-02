# Compare pairwise velocity distributions from particle catalogues
# v0.1.0, 2021-09-15 - Code started from plot_vp_mass_comp_abacus.py

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


def calc_v12_lin(ds, H0, ombh2, omch2, ns, sigma8, z, kmax=2.0, minkh=1e-4, maxkh=1, npoints=200):
    # Assumes ds is a numpy.array
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(ns=ns)
    pars.set_matter_power(redshifts=[z], kmax=kmax)

    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    kh, zcamb, pk = results.get_matter_power_spectrum(minkh=minkh, maxkh=maxkh, npoints=npoints)
    sigma8_default = np.array(results.get_sigma8())
    omegam0 = (ombh2 + omch2) / (H0 / 100.)**2.
    sigma8 = sigma8 * calc_D(z, omegam0, 1. - omegam0) / calc_D(0., omegam0, 1. - omegam0)
    pk = pk * (sigma8 / sigma8_default[0])**2.

    v12s = np.zeros(ds.shape)
    for i in range(ds.size):
        for j in range(kh.size):
            if j==kh.size-1:
                log_dk = np.log10(kh[j]) - np.log10(kh[j-1])
                k_next = np.power(10., (np.log10(kh[j]) + log_dk))
                dk = (k_next - kh[j-1]) / 2.
            else:
                dk = (kh[j+1] - kh[j-1]) / 2.
            v12s[i] += dk * kh[j] * pk[0][j] * spherical_jn(1, ds[i]*kh[j])
    # H = H_z(z, H0, omegam0)
    # v12s = -H * v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
    v12s = -H0 * v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
    return v12s


cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/cosmo_params.dat")

z = 0.7
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))

H_seps = np.logspace(np.log10(0.01), np.log10(100.), 200)
v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)

plt.figure(figsize=(8., 6.))
plt.axhline(y=0., linestyle="-", color="k")
ymax = 100.
ymin = -500.
H = H_z(z, cosmo_params[0], (cosmo_params[1] + cosmo_params[2])/(cosmo_params[0]/100.)**2.)
plt.plot(H_seps, -H*H_seps, color="k", label="Static solution")

ngc_mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_100000_no_growth_corr_vp_v5.dat")
ngc_mvps[:, 1] = np.nan_to_num(ngc_mvps[:, 1])
plt.errorbar(bin_centres, ngc_mvps[:, 1], yerr=np.sqrt(ngc_mvps[:, 2] / ngc_mvps[:, 3]), color="C1", label="Growth factor corrected z=49 IC particles /w corrected disp", marker=".", linestyle="-")

lr_mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_100000_low_redshift_vp_v4.dat")
lr_mvps[:, 1] = np.nan_to_num(lr_mvps[:, 1])
plt.errorbar(bin_centres, lr_mvps[:, 1], yerr=np.sqrt(lr_mvps[:, 2] / lr_mvps[:, 3]), color="C4", label="z=0.7 IC particles", marker=".", linestyle="-")

'''
mean_scaling = np.mean(lr_mvps[-20:] / ngc_mvps[-20:])
print("Mean amplitude different >10 Mpc/h is ", mean_scaling)
plt.errorbar(bin_centres, mean_scaling*ngc_mvps[:, 1], yerr=mean_scaling*np.sqrt(lr_mvps[:, 2] / lr_mvps[:, 3]), color="C2", label="Scaled z=49 IC particles", marker=".", linestyle="--")

z = 0.7
omch2_planck = 0.1199
ombh2_planck = 0.02222
omh2_planck = 0.14212
h_planck = 0.6726
ns_planck = 0.9652
sigma8_planck = 0.830
As_default = 2e-09
sigma8_default = 0.7917199
As_rescale = (sigma8_planck / sigma8_default) ** 2.

pars = camb.CAMBparams()
pars.set_cosmology(H0=h_planck*100., ombh2=ombh2_planck, omch2=omch2_planck)
pars.InitPower.set_params(ns=ns_planck, As=As_default*As_rescale)
pars.set_matter_power(redshifts=[0., 0.7, 49.], kmax=2.0)

pars.NonLinear = camb.model.NonLinear_none
results = camb.get_results(pars)
kh, zcamb, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints=200)
sigma8 = np.array(results.get_sigma8())
print("Camb sigma8(z=0)=", sigma8[2])

fsigma8_0p7 = np.array(results.get_fsigma8())[1]
fsigma8_49 = np.array(results.get_fsigma8())[0]
revised_scaling = fsigma8_49 / fsigma8_0p7 * sigma8[1] / sigma8[0]
print("Revised Scaling: ", revised_scaling)

plt.errorbar(bin_centres, revised_scaling*ngc_mvps[:, 1], yerr=mean_scaling*np.sqrt(lr_mvps[:, 2] / lr_mvps[:, 3]), color="C3", label="Growth Factor only z=49 IC particles", marker=".", linestyle=":")
'''

v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z)
plt.plot(v12_seps, v12s, color="k", linestyle="--", label="Linear theory prediction")


plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p$ [km/s]")
plt.xlim(0.01, 100.)
plt.ylim(ymin-50., ymax+50.)
plt.xscale("log")
plt.legend()
plt.savefig("../output/plots/abacus_part_vp_comp_plot_v5.jpg")


# Change Log

