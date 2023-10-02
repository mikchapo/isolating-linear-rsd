# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from pdf2image import convert_from_path
from PIL import Image
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
    else:
        pars.set_matter_power(redshifts=[0], kmax=kmax)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)

    kh, zcamb, pk = results.get_matter_power_spectrum(minkh=minkh, maxkh=maxkh, npoints=npoints)
    print(pk.shape)
    print(pk)


v12_seps = np.logspace(np.log10(1.), np.log10(100.), 200)

z = 0.7
H0 = 73.
ommh2 = 0.25 * (H0 / 100.)**2.
ombh2 = 0.045 * (H0 / 100.)**2.
omch2 = ommh2 - ombh2
ns = 1.
sigma8 = 0.9
v12s = calc_v12_lin(v12_seps, H0, ombh2, omch2, ns, sigma8, z)