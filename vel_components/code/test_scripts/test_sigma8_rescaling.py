import camb
import numpy as np
from scipy.integrate import quad


def D_integrand(a, omegam, omegal):
    return ((omegam / (a**3.) / (omegam / (a**3.) + omegal))**0.55 - 1.) / a


def calc_D(z, omegam0, omegal):
    a = 1. / (1. + z)
    integral = quad(D_integrand, 0., a, args=(omegam0, omegal))
    return a*np.exp(integral[0])


z1 = 0.
z2 = 0.7
ombh2 = 0.02222
omch2 = 0.1199
h = 0.6726
omegam0 = (ombh2 + omch2) / h**2
sigma8_1 = 0.83
sigma8_2 = sigma8_1 / calc_D(z1, omegam0, 1 - omegam0) * calc_D(z2, omegam0, 1
                                                                - omegam0)
print("My sigma8(z={}) =".format(z2), sigma8_2)

omh2 = 0.14212
ns = 0.9652
As_default = 2e-09
sigma8_default = 0.7917199
As_rescale = (sigma8_1 / sigma8_default) ** 2.

pars = camb.CAMBparams()
pars.set_cosmology(H0=h*100., ombh2=ombh2, omch2=omch2)
pars.InitPower.set_params(ns=ns, As=As_default*As_rescale)
pars.set_matter_power(redshifts=[z1, z2], kmax=2.0)

pars.NonLinear = camb.model.NonLinear_none
results = camb.get_results(pars)
sigma8 = np.array(results.get_sigma8())
print("Camb sigma8(z={}) =".format(z1), sigma8[1])
print("Camb sigma8(z={}) =".format(z2), sigma8[0])
