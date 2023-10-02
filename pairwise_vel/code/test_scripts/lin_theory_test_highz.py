# Test the dependence of the high-z lin theory prediction on maxkh
# v0.1.0, 2021-09-15 - Code started from plot_vp_part.py

# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn
import sys


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
    print("sigma8_default:", sigma8_default)
    v12s = -H0 * v12s * fsigma8 / sigma8_default / np.pi**2.
    return v12s

cosmo_params = np.loadtxt("../../data/cosmo_params.dat")
z = 20.
v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)
maxkhs = [0.25, 0.5, 1, 2, 4]

plt.figure(figsize=(8., 6.))
plt.axhline(y=0., linestyle="-", color="k")

for maxkh in maxkhs:
    v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], z, maxkh=maxkh)
    plt.plot(v12_seps, v12s, linestyle="-", label="maxkh = {}".format(maxkh))

plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p$ [km/s]")
plt.xlim(0.01, 100.)
plt.xscale("log")
plt.legend()
plt.savefig("../../output/plots/test/lin_theory_test_z20.jpg")