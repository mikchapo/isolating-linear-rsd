# Compare pairwise velocity from two test sims
# v0.1.0, 2021-06-18 - Code Started

# Imports
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import sys


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


sim_box1 = sys.argv[1]
sim_box2 = sys.argv[2]

bin_centres = np.arange(0.5, 150.5, 1.)
mass_bin_lims = np.arange(10.75, 15.25, 0.25)

mass_bin = 10 # 13.00 < logM < 13.25

sigma8_1 = 6.944280e-01
h_1 = 6.32317e+01 / 100.
omegam0_1 = (2.32629e-02 + 1.07830e-01) / h_1**2.
sigma8_2 = 7.542319e-01
h_2 = 6.57317e+01 / 100.
omegam0_2 = (2.27629e-02 + 1.12830e-01) / h_1**2.
fsig8_1 = fsigma8_approximate(0.7, sigma8_1, omegam0_1)
fsig8_2 = fsigma8_approximate(0.7, sigma8_2, omegam0_2)
# gamma_f = fsig8_1 / fsig8_2

gamma_f = (omegam(0.7, omegam0_1) / omegam(0.7, omegam0_2))**0.55

print("Approx. gamma_f:", gamma_f)
# zs = np.linspace(0., 100., 100001)
# fsig8_1_num = fsigma8_numerical(zs, sigma8_1, omegam0_1)
# fsig8_2_num = fsigma8_numerical(zs, sigma8_2, omegam0_2)
# gamma_f_num = fsig8_1_num[700] / fsig8_2_num[700]
# print("Numerical gamma_f:", gamma_f_num)

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [2, 1]})
mvps1 = np.empty((5, 150))
mvps2 = np.empty((5, 150))
for i in range(5):
    mvps1[i, :] = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_{}-00{}/m200b_{:.2f}-{:.2f}_vp.dat".format(sim_box1, i, mass_bin_lims[mass_bin], mass_bin_lims[mass_bin+1]))[:, 1]
    mvps2[i, :] = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/aemulus_z0.70_{}-00{}/m200b_{:.2f}-{:.2f}_vp.dat".format(sim_box2, i, mass_bin_lims[mass_bin], mass_bin_lims[mass_bin+1]))[:, 1]
mean_mvps1 = np.mean(mvps1, axis=0)
std_mvps1 = np.std(mvps1, axis=0)
mean_mvps2 = np.mean(mvps2, axis=0)
std_mvps2 = np.std(mvps2, axis=0)
axes[0].errorbar(bin_centres, mean_mvps1, yerr=std_mvps1, label=sim_box1, color="C0")
axes[0].errorbar(bin_centres, mean_mvps2, yerr=std_mvps2, label=sim_box2, color="C1")
# axes[0].errorbar(bin_centres, gamma_f*mean_mvps2, yerr=gamma_f*std_mvps2, label=r"$\gamma_f$" + sim_box2, color="C2")
# axes[0].errorbar(bin_centres, gamma_f_num*mean_mvps2, yerr=gamma_f_num*std_mvps2, label=r"Num. $\gamma_f$" + sim_box2, color="C3")
# print(mean_mvps2/mean_mvps1)
ratio_err = mean_mvps2 / mean_mvps1 * np.sqrt(np.power(std_mvps1/mean_mvps1, 2.) + np.power(std_mvps2/mean_mvps2, 2.))
mean_ratio_lin = np.mean(mean_mvps2[7:]/mean_mvps1[7:])
mean_ratio_fit = np.mean(mean_mvps2[:60]/mean_mvps1[:60])
mean_ratio_quasi_lin = np.mean(mean_mvps2[7:60]/mean_mvps1[7:60])
print("Mean Ratio 7-150", mean_ratio_lin)
print("Mean Ratio 0-60", mean_ratio_fit)
print("Mean Ratio 7-60", mean_ratio_quasi_lin)
# axes[1].errorbar(bin_centres, mean_mvps2/mean_mvps1, yerr=ratio_err, color="C1")
axes[1].plot(bin_centres, mean_mvps2/mean_mvps1, color="C1")
axes[1].fill_between(bin_centres, mean_mvps2/mean_mvps1 - ratio_err, mean_mvps2/mean_mvps1 + ratio_err, color="C1", alpha=0.5)
axes[1].plot(bin_centres, gamma_f*mean_mvps2/mean_mvps1, color="C2")
axes[1].fill_between(bin_centres, gamma_f*mean_mvps2/mean_mvps1 - gamma_f*ratio_err, gamma_f*mean_mvps2/mean_mvps1 + gamma_f*ratio_err, color="C2", alpha=0.5)
# axes[1].errorbar(bin_centres, gamma_f*mean_mvps2 - mean_mvps1, yerr=gamma_f*std_mvps2, color="C2")
# axes[1].errorbar(bin_centres, gamma_f_num*mean_mvps2 - mean_mvps1, yerr=gamma_f_num*std_mvps2, color="C3")
axes[1].axhline(y=1, color="k")
axes[1].axhline(y=mean_ratio_lin, color="k", linestyle="--")
axes[0].axvline(x=7, linestyle="--", color="k")
axes[1].axvline(x=7, linestyle="--", color="k")
axes[1].set_xlabel(r"s [$h^{-1}$ Mpc]")
axes[0].set_ylabel(r"$v_p$ [km/s]")
axes[1].set_ylabel(r"$v_p^2~/~v_p^1$")
axes[1].set_xlim(0., 150.)
axes[0].legend()
plt.savefig("../output/plots/aemulus_z0.70_{}_{}_M_{}_comp.jpg".format(sim_box1, sim_box2, mass_bin))