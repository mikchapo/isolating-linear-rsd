from scipy.integrate import odeint
from scipy.integrate import quad

import matplotlib.pyplot as plt
import numpy as np


def growth_factor(y, a, omegam0):
	delta, delta_prime = y
	H2_H02 = omegam0 / (a**3.) + (1. - omegam0)
	dyda = [delta_prime, ((omegam0 / 2 / (a**3.) / H2_H02 - 1.) * delta_prime + omegam0 * delta / 2. / (a**4.) / H2_H02) * 3 / a]
	return dyda


def fsigma8_numerical(zs, sigma8=0.811, omegam0=0.315):
	scale_factors = np.flip(1. / (1. + zs))
	y0 = [scale_factors[0], 1.]
	deltas = odeint(growth_factor, y0, scale_factors, args=(omegam0,))
	fsig8s = np.flip(sigma8 * scale_factors * deltas[:, 1] / deltas[-1, 0])
	return fsig8s


def omegam(z, omegam0=0.315):
	return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
	l = 0.83 + 0.19 / (omegam0**0.8)
	beta = 0.98 + 0.01 / (omegam0**1.2)
	gamma = 0.64 + 0.09 / (omegam0**0.4)
	print("lambda:", l, "beta:", beta, "gamma:", gamma)
	return l * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


fsig8_approx = fsigma8_approximate(0.7)
print("Approximate Result:", fsig8_approx)

zs = np.linspace(0., 100., 100001)
fisg8s_num_lin_full = 