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
    As_default = 2e-09

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(ns=ns)
    if z!=0:
        pars.set_matter_power(redshifts=[0, z], kmax=kmax)
    else:
        pars.set_matter_power(redshifts=[0], kmax=kmax)

    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    sigma8_0_default = np.array(results.get_sigma8_0())
    As_rescale = (sigma8 / sigma8_0_default) ** 2.

    pars.InitPower.set_params(ns=ns, As=As_default*As_rescale)
    if z!=0:
        pars.set_matter_power(redshifts=[0, z], kmax=kmax)
        pk_index = 1
    else:
        pars.set_matter_power(redshifts=[0], kmax=kmax)
        pk_index = 0
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)

    kh, zcamb, pk = results.get_matter_power_spectrum(minkh=minkh, maxkh=maxkh, npoints=npoints)

    v12s = np.zeros(ds.shape)
    for i in range(ds.size):
        for j in range(kh.size):
            if j==kh.size-1:
                log_dk = np.log10(kh[j]) - np.log10(kh[j-1])
                k_next = np.power(10., (np.log10(kh[j]) + log_dk))
                dk = (k_next - kh[j-1]) / 2.
            else:
                dk = (kh[j+1] - kh[j-1]) / 2.
            v12s[i] += dk * kh[j] * pk[pk_index][j] * spherical_jn(1, ds[i]*kh[j])
    # H = H_z(z, H0, omegam0)
    # v12s = -H * v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
    fsigma8 = np.array(results.get_fsigma8())[0]
    sigma8_default = np.array(results.get_sigma8())[0]
    print("fsigma8:", fsigma8)
    print("f camb:", fsigma8 / sigma8_default)
    print("fsigma8_approx:", omegam_calc(z, (omch2+ombh2)/(H0/100.)**2)**0.55)
    print("sigma8_default:", sigma8_default)
    v12s = -H0 * v12s * fsigma8 / sigma8_default / np.pi**2.
    return v12s

# def calc_v12_lin(ds, H0, ombh2, omch2, ns, sigma8, z, kmax=2.0, minkh=1e-4, maxkh=1, npoints=200):
#     # Assumes ds is a numpy.array
#     pars = camb.CAMBparams()
#     pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
#     pars.InitPower.set_params(ns=ns)
#     pars.set_matter_power(redshifts=[z], kmax=kmax)

#     pars.NonLinear = camb.model.NonLinear_none
#     results = camb.get_results(pars)
#     kh, zcamb, pk = results.get_matter_power_spectrum(minkh=minkh, maxkh=maxkh, npoints=npoints)
#     sigma8_0_default = np.array(results.get_sigma8_0())
#     omegam0 = (ombh2 + omch2) / (H0 / 100.)**2.
#     pk = pk * (sigma8 / sigma8_0_default[0])**2.

#     v12s = np.zeros(ds.shape)
#     for i in range(ds.size):
#         for j in range(kh.size):
#             if j==kh.size-1:
#                 log_dk = np.log10(kh[j]) - np.log10(kh[j-1])
#                 k_next = np.power(10., (np.log10(kh[j]) + log_dk))
#                 dk = (k_next - kh[j-1]) / 2.
#             else:
#                 dk = (kh[j+1] - kh[j-1]) / 2.
#             v12s[i] += dk * kh[j] * pk[0][j] * spherical_jn(1, ds[i]*kh[j])
#     # H = H_z(z, H0, omegam0)
#     # v12s = -H * v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
#     v12s = -H0 * v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
#     return v12s


cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/cosmo_params.dat")

z = float(sys.argv[1])
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))

v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)

plt.figure(figsize=(8., 6.))
plt.axhline(y=0., linestyle="-", color="k")

vel_scaling = 0.9957506478771345

# rescaling =  0.5784773960932065 / 0.5781453735795127

ngc_nplt_mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_ngc_nplt_vp_v6.dat".format(z))
ngc_nplt_mvps[:, 1] = vel_scaling * np.nan_to_num(ngc_nplt_mvps[:, 1])

ngc_mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_no_growth_corr_vp_v6.dat".format(z))
ngc_mvps[:, 1] = vel_scaling * np.nan_to_num(ngc_mvps[:, 1])
plt.errorbar(bin_centres, ngc_mvps[:, 1] - ngc_nplt_mvps[:, 1], yerr=np.sqrt(ngc_mvps[:, 2] / ngc_mvps[:, 3]), color="C0", label="ic_no_growth_corr", linestyle="-")
# plt.plot(bin_centres, ngc_mvps[:, 1] - ngc_nplt_mvps[:, 1], color="C0", label="ic_no_growth_corr", linestyle="-")

mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_32c_vp_v6.dat".format(z))
mvps[:, 1] = vel_scaling * np.nan_to_num(mvps[:, 1])
# plt.errorbar(bin_centres, mvps[:, 1] - ngc_nplt_mvps[:, 1], yerr=np.sqrt(mvps[:, 2] / mvps[:, 3]), color="C1", label="ic_32c", linestyle="--")
plt.plot(bin_centres, mvps[:, 1] - ngc_nplt_mvps[:, 1], color="C1", label="ic_32c", linestyle="--")

ngc_nplt_mvps_z0p7 = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_ngc_nplt_z0.7_vp_v6.dat".format(z))
ngc_nplt_mvps_z0p7[:, 1] = np.nan_to_num(ngc_nplt_mvps_z0p7[:, 1])
# plt.errorbar(bin_centres, ngc_nplt_mvps_z0p7[:, 1] - ngc_nplt_mvps[:, 1], yerr=np.sqrt(ngc_nplt_mvps_z0p7[:, 2] / ngc_nplt_mvps_z0p7[:, 3]), color="C2", label="ic_ngc_nplt_z0.7", linestyle=":")
plt.plot(bin_centres, ngc_nplt_mvps_z0p7[:, 1] - ngc_nplt_mvps[:, 1], color="C2", label="ic_ngc_nplt_z0.7", linestyle=":")

ngc_mvps_z0p7 = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_no_growth_corr_z0.7_vp_v6.dat".format(z))
ngc_mvps_z0p7[:, 1] = np.nan_to_num(ngc_mvps_z0p7[:, 1])
plt.errorbar(bin_centres, ngc_mvps_z0p7[:, 1] - ngc_nplt_mvps[:, 1], yerr=np.sqrt(ngc_mvps_z0p7[:, 2] / ngc_mvps_z0p7[:, 3]), color="C3", label="ic_no_growth_corr_z0.7", linestyle="-")
# plt.plot(bin_centres, ngc_mvps_z0p7[:, 1] - ngc_nplt_mvps[:, 1], color="C3", label="ic_no_growth_corr_z0.7", linestyle="-")

mvps_z0p7 = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_z0.7_vp_v6.dat".format(z))
mvps_z0p7[:, 1] = np.nan_to_num(mvps_z0p7[:, 1])
# plt.errorbar(bin_centres, mvps_z0p7[:, 1] - ngc_nplt_mvps[:, 1], yerr=np.sqrt(mvps_z0p7[:, 2] / mvps_z0p7[:, 3]), color="C4", label="ic_z0.7", linestyle="--")
plt.plot(bin_centres, mvps_z0p7[:, 1] - ngc_nplt_mvps[:, 1], color="C4", label="ic_z0.7", linestyle="--")



# v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z)
# plt.plot(v12_seps, v12s - , color="k", linestyle="--", label="Linear theory prediction")


plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p - v_p^R$ [km/s]")
# plt.xlim(0.01, 100.)
plt.xlim(1, 100.)
plt.ylim(-5, 5)
plt.xscale("log")
plt.legend()
plt.savefig("../output/plots/abacus_part_z{}_ic_settings_vp_comp_plot_v6_scaled.jpg".format(z))


# Change Log
