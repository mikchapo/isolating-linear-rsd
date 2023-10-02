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


def interpolator_intergrator(k, z, d, interpolator, k_hunit):
    if k_hunit:
        return k * interpolator.P(z, k) * spherical_jn(1, d*k)
    else:
        return k * interpolator.P(z, k) * spherical_jn(1, d*k/0.6726)

def calc_v12_lin(ds, H0, ombh2, omch2, ns, sigma8, z, kmax=2.0, minkh=1e-4, maxkh=1, npoints=200, use_trapz=False, use_interpolator=False, k_hunit=True, hubble_units=True):
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

    if use_interpolator:
        interpolator = camb.get_matter_power_interpolator(pars, zs=[z], nonlinear=False, hubble_units=hubble_units, k_hunit=k_hunit)
        v12s = np.zeros(ds.shape)
        for i in range(ds.size):
            if k_hunit:
                v12s[i] = quad(interpolator_intergrator, minkh, maxkh, args=(z, ds[i], interpolator, k_hunit))[0]
            else:
                v12s[i] = quad(interpolator_intergrator, minkh*0.6726, maxkh*0.6726, args=(z, ds[i], interpolator, k_hunit))[0]


    else:
        kh, zcamb, pk = results.get_matter_power_spectrum(minkh=minkh, maxkh=maxkh, npoints=npoints)

        if use_trapz:
            v12s = np.zeros(ds.shape)
            for i in range(ds.size):
                integrand = kh * pk[pk_index] * spherical_jn(1, ds[i]*kh)
                v12s[i] = np.trapz(integrand, x=kh)

        else:
            # kh = kh*0.6726
            kh = kh
            v12s = np.zeros(ds.shape)
            for i in range(ds.size):
                for j in range(kh.size):
                    if j==kh.size-1:
                        log_dk = np.log10(kh[j]) - np.log10(kh[j-1])
                        k_next = np.power(10., (np.log10(kh[j]) + log_dk))
                        dk = (k_next - kh[j-1]) / 2.
                    else:
                        dk = (kh[j+1] - kh[j-1]) / 2.
                    # v12s[i] += dk * kh[j] * pk[pk_index][j] * spherical_jn(1, ds[i]*kh[j]/0.6726)
                    v12s[i] += dk * kh[j] * pk[pk_index][j] * spherical_jn(1, ds[i]*kh[j])
    # H = H_z(z, H0, omegam0)
    # v12s = -H * v12s * omegam_calc(z, omegam0)**0.55 / np.pi**2.
    fsigma8 = np.array(results.get_fsigma8())[0]
    sigma8_default = np.array(results.get_sigma8())[0]
    print("fsigma8:", fsigma8)
    print("f camb:", fsigma8 / sigma8_default)
    print("fsigma8_approx:", omegam_calc(z, (omch2+ombh2)/(H0/100.)**2)**0.55)
    print("sigma8_default:", sigma8_default)
    # v12s = -H0 * v12s * fsigma8 / sigma8_default / np.pi**2.
    H = results.hubble_parameter(z)
    v12s = -H / (z + 1.) * v12s * fsigma8 / sigma8_default / np.pi**2.
    return v12s


def calc_v12_abacus(ds, H0, k_offset=0):
    camb_matterpower = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/camb_matterpower.dat")
    kh = camb_matterpower[:, 0]
    # pk = camb_matterpower[:, 1] * (0.0210839935761 / 0.83) ** 2.
    pk = camb_matterpower[:, 1]
    v12s = np.zeros(ds.shape)
    for i in range(ds.size):
        for j in range(kh.size - k_offset):
            if j==kh.size-1:
                log_dk = np.log10(kh[j]) - np.log10(kh[j-1])
                k_next = np.power(10., (np.log10(kh[j]) + log_dk))
                dk = (k_next - kh[j-1]) / 2.
            else:
                dk = (kh[j+1] - kh[j-1]) / 2.
            v12s[i] += dk * kh[j] * pk[j] * spherical_jn(1, ds[i]*kh[j])
    # v12s = -H0 * v12s / np.pi**2.
    v12s = -H0 * v12s / np.pi**2. / 100.
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
lr_ic_dir = sys.argv[2]
N = int(sys.argv[3])
plot_lr = True
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))

H_seps = np.logspace(np.log10(0.01), np.log10(100.), 200)
v12_seps = np.logspace(np.log10(1.), np.log10(100.), 40)

plt.figure(figsize=(8., 6.), dpi=300)
plt.axhline(y=0., linestyle="-", color="k")
ymax = 100.
ymin = -500.
H = H_z(z, cosmo_params[0], (cosmo_params[1] + cosmo_params[2])/(cosmo_params[0]/100.)**2.)
# plt.plot(H_seps, -H*H_seps, color="k", label="Static solution")

ngc_mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_{}_ic_no_growth_corr_vp_v6.dat".format(z, N))
ngc_mvps[:, 1] = np.nan_to_num(ngc_mvps[:, 1])
plt.errorbar(bin_centres, ngc_mvps[:, 1], yerr=np.sqrt(ngc_mvps[:, 2] / ngc_mvps[:, 3]), color="C1", label="Growth rate corrected z=49 IC particles /w corrected disp", marker=".", linestyle="-")

if plot_lr:
    lr_mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_{}_{}_vp_v6.dat".format(z, N, lr_ic_dir))
    lr_mvps[:, 1] = np.nan_to_num(lr_mvps[:, 1])
    plt.errorbar(bin_centres, lr_mvps[:, 1], yerr=np.sqrt(lr_mvps[:, 2] / lr_mvps[:, 3]), color="C4", label="z={} IC particles /w f scaling".format(z), marker=".", linestyle="-")


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

# v12s_trapz = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_trapz=True)
# plt.plot(v12_seps, v12s_trapz, color="C0", linestyle="-", label="Linear theory prediction, using trapz")


# if z == 0:
#     v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z)
#     plt.plot(v12_seps, v12s, color="k", linestyle="--", label="Linear theory prediction, maxkh=1")
#     v12s_maxkh = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, maxkh=0.2, use_interpolator=True)
#     plt.plot(v12_seps, v12s_maxkh, color="C2", linestyle=":", label="Linear theory prediction, maxkh=0.2")


v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=True, hubble_units=True, k_hunit=True)
plt.plot(v12_seps, v12s, color="C0", linestyle="-", label="Linear theory prediction, camb units k[h/Mpc] and P[(Mpc/h)^3]")

# v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=True, hubble_units=False, k_hunit=True)
# plt.plot(v12_seps, v12s, color="C2", linestyle="--", label="Linear theory prediction, k_hunit")

# v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=True, hubble_units=True, k_hunit=False)
# plt.plot(v12_seps, v12s, color="C3", linestyle="-.", label="Linear theory prediction, hubble_units")

v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=True, hubble_units=False, k_hunit=False)
plt.plot(v12_seps, v12s, color="C1", linestyle="--", label="Linear theory prediction, camb units k[1/Mpc] and P[Mpc^3]")

v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=True, hubble_units=True, k_hunit=True)
plt.plot(v12_seps, v12s / cosmo_params[0] * 100., color="C2", linestyle=":", label="Linear theory prediction divided by h, camb units k[h/Mpc] and P[(Mpc/h)^3]")

# v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=False, use_trapz=False)
# plt.plot(v12_seps, v12s, color="C5", linestyle=None, marker="o", label="Linear theory prediction, extra h")

# v12s_trapz = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, use_interpolator=True)
# plt.plot(v12_seps, v12s_trapz, color="C2", linestyle=":", label="Linear theory prediction, using interpolator")

# v12s_abacus_trunc = calc_v12_abacus(v12_seps, cosmo_params[0], k_offset=50)
# plt.plot(v12_seps, v12s_abacus_trunc, color="C2", linestyle="-.", label="Linear theory prediction, truncated Abacus power spectrum")

# v12s_abacus_full = calc_v12_abacus(v12_seps, cosmo_params[0])
# plt.plot(v12_seps, v12s_abacus_full, color="C3", linestyle=":", label="Linear theory prediction, full Abacus power spectrum")

# print("vp ratios from s=1-100 Mpc/h:")
# print(ngc_mvps[40:, 1] / v12s)
# print("Mean ratio from s=10-100 Mpc.h")
# print(np.mean(ngc_mvps[60:, 1] / v12s[20:]))

plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p$ [km/s]")
# plt.xlim(0.01, 100.)
if not plot_lr:
    plt.xlim(1., 100.)
    plt.ylim(-2., 2.)
else:
    plt.xlim(0.01, 100.)
plt.xscale("log")
plt.legend()
plt.savefig("../output/plots/abacus_part_z{}_{}_vp_comp_plot_v7_camb_units.jpg".format(z, N))


# Change Log

