# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from pdf2image import convert_from_path
from PIL import Image
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
    else:
        pars.set_matter_power(redshifts=[0], kmax=kmax)
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
            v12s[i] += dk * kh[j] * pk[0][j] * spherical_jn(1, ds[i]*kh[j])
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


v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)

# z = 0.7
# test_cosmos = np.loadtxt("../../../RSD/aemulus_fmax/cosmology_camb_test_box_full.dat")
# j = 0
# ombh2 = test_cosmos[j, 1] * test_cosmos[j, 3]**2.
# omch2 = (test_cosmos[j, 0] - test_cosmos[j, 1]) * test_cosmos[j, 3]**2.
# v12s = calc_v12_lin(v12_seps, test_cosmos[j, 3]*100., ombh2, omch2, test_cosmos[j, 4], test_cosmos[j, 2], z)

z = 0.
H0 = 73.
ommh2 = 0.25 * (H0 / 100.)**2.
ombh2 = 0.045 * (H0 / 100.)**2.
omch2 = ommh2 - ombh2
ns = 1.
sigma8 = 0.9
v12s = calc_v12_lin(v12_seps, H0, ombh2, omch2, ns, sigma8, z)

belloso_image = convert_from_path('../data/guo2010_vinfall_np0_z0_noiso_theory.pdf')[0]

im_width, im_height = belloso_image.size
print("Belloso size:", im_width, im_height)
print("My size:", 8*200, 6.5*200)

scaling = 1.

# hyperparams = [{'name': 'kmax', 'vals': [0.5, 1.0, 2.0, 4.0, 8.0]},
#                {'name': 'minkh', 'vals': [1e-5, 5e-5, 1e-4, 5e-4, 1e-3]},
#                {'name': 'maxkh', 'vals': [0.25, 0.5, 1, 2, 4]},
#                {'name': 'npoints', 'vals': [50, 100, 200, 400, 800]}]

# # hyperparams = [{'name': 'maxkh', 'vals': [0.25, 0.5, 1, 2, 4]}]

# linestyles = ["-", "--", "-.", (0, (3, 1, 1, 1)), ":"]

# for hyperparam in hyperparams:
#     plt.figure(figsize=(7.93, 6.38), dpi=200)
#     plt.axhline(y=0., linestyle="-", color="k")
#     for i in range(len(hyperparam["vals"])):
#         kwarg = {hyperparam['name']: hyperparam["vals"][i]}

#         if hyperparam["name"] == "maxkh":
#             print("Yes")
#             v12s = calc_v12_lin(v12_seps, H0, ombh2, omch2, ns, sigma8, z, kmax=2.*hyperparam["vals"][i], **kwarg)
#         else:
#             v12s = calc_v12_lin(v12_seps, H0, ombh2, omch2, ns, sigma8, z, **kwarg)

#         plt.plot(v12_seps, scaling*v12s, linestyle=linestyles[i], color="C{}".format(i), label="{}={}".format(hyperparam['name'], hyperparam["vals"][i]))
#     plt.xlabel(r"s [$h^{-1}$ Mpc]")
#     plt.ylabel(r"$v_p$ [km/s]")
#     plt.xlim(0.005, 100.)
#     plt.ylim(-350, 50.)
#     plt.xscale("log")
#     plt.tight_layout()
#     plt.legend()
#     plt.savefig("../output/plots/lin_theory_test_{}.png".format(hyperparam['name']))

# v12s_1000 = calc_v12_lin(v12_seps, H0, ombh2, omch2, ns, sigma8, z, npoints=1000)

plt.figure(figsize=(7.93, 6.38), dpi=200)
# plt.axhline(y=0., linestyle="-", color="k")
plt.plot(v12_seps, scaling*v12s, linestyle=":", color="k")
# plt.plot(v12_seps, scaling*v12s_1000, linestyle="--", color="b")
# plt.xlabel(r"s [$h^{-1}$ Mpc]")
# plt.ylabel(r"$v_p$ [km/s]")
plt.xlim(0.005, 100.)
plt.ylim(-350, 50.)
plt.xscale("log")
plt.tight_layout()
plt.savefig("../output/plots/lin_theory_test.png", transparent=True)

my_image = Image.open("../output/plots/lin_theory_test.png")
belloso_image.paste(my_image, (20, -12), my_image)
belloso_image.save("../output/plots/lin_theory_comp_v2.jpg")