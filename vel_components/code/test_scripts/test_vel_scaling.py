import numpy as np


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


omh2_planck = 0.14212
h_planck = 0.6726
omegam0_planck = omh2_planck * h_planck**2.
H_0p7 = H_z(0.7, h_planck*100., omegam0_planck)
H_49 = H_z(49., h_planck*100., omegam0_planck)
sigma8_planck = 0.830

fsigma8_0p7 = fsigma8_approximate(0.7, sigma8=sigma8_planck, omegam0=omegam0_planck)
fsigma8_49 = fsigma8_approximate(49., sigma8=sigma8_planck, omegam0=omegam0_planck)

vel_scaling = 50. / 1.7 * H_0p7 / H_49 * fsigma8_0p7 / fsigma8_49

print("H(0.7) =", H_0p7)
print("H(49) =", H_49)
print("fsig8(0.7) =", fsigma8_0p7)
print("fsig8(49) =", fsigma8_49)
print("Total Scaling:", vel_scaling)