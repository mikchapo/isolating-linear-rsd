import camb
import numpy as np


cosmo_params = np.loadtxt("../data/cosmo_params.dat")

def calc_f(H0, ombh2, omch2, ns, sigma8, z, kmax=2.0, minkh=1e-4, maxkh=1, npoints=200):
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

    fsigma8 = np.array(results.get_fsigma8())[0]
    sigma8_default = np.array(results.get_sigma8())[0]
    print("fsigma8(z={}):".format(z), fsigma8)
    print("f(z={}) camb:".format(z), fsigma8 / sigma8_default)
    return fsigma8 / sigma8_default


v12s = calc_f(cosmo_params[0], cosmo_params[2], cosmo_params[1], cosmo_params[3], cosmo_params[4], 49.)