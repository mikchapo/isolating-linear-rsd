# Imports
import camb
import numpy as np
from scipy.special import spherical_jn


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

    v12s = -H0 * v12s * fsigma8 / sigma8_default / np.pi**2.
    return v12s