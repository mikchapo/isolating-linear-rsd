import camb
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn


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
    print("Power Spectrum z:", zcamb)
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


cosmo_params = np.loadtxt("../../data/cosmo_params.dat")
v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)
zs = [0., 0.3, 0.7, 1., 1.5, 5., 10., 20., 49.]

plt.figure(figsize=(8., 6.))
plt.axhline(y=0., linestyle="-", color="k")

for i in range(len(zs)):
    v12s = calc_v12_lin(v12_seps, cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], zs[i])
    plt.plot(v12_seps, v12s, color="C{}".format(i), linestyle="-", label="z={}".format(zs[i]))

plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p$ [km/s]")
plt.xlim(0.01, 100.)
plt.xscale("log")
plt.legend()
plt.savefig("../../output/plots/test/vp_lin_z49.jpg")
