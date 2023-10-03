import camb
import numpy as np
import os
from struct import iter_unpack, unpack_from
import sys


def H_z(z, H0, omegam0):
    return H0 * np.sqrt(omegam0 * (1. + z) * (1. + z) * (1. + z) + (1. - omegam0))


def omegam(z, omegam0=0.315):
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
    l = 0.83 + 0.19 / (omegam0**0.8)
    beta = 0.98 + 0.01 / (omegam0**1.2)
    gamma = 0.64 + 0.09 / (omegam0**0.4)
    # print("lambda:", l, "beta:", beta, "gamma:", gamma)
    return l * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


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

print("Camb H(z=0.7)=", results.hubble_parameter(0.7))
print("Camb H(z=49)=", results.hubble_parameter(49.))
print("Camb fsig8(z=0.7)=", np.array(results.get_fsigma8())[1])
print("Camb fsig8(z=49)=", np.array(results.get_fsigma8())[0])

z = 0.7
omh2_planck = 0.14212
h_planck = 0.6726
omegam0_planck = omh2_planck / h_planck**2.
H_0p7 = H_z(0.7, h_planck*100., omegam0_planck)
H_49 = H_z(49., h_planck*100., omegam0_planck)

fsigma8_0p7 = fsigma8_approximate(0.7, sigma8=sigma8_planck, omegam0=omegam0_planck)
fsigma8_49 = fsigma8_approximate(49., sigma8=sigma8_planck, omegam0=omegam0_planck)

print("Approx. H(z=0.7)=", H_0p7)
print("Approx. H(z=49)=", H_49)
print("Approx. fsig8(z=0.7)=", fsigma8_0p7)
print("Approx. fsig8(z=49)=", fsigma8_49)

