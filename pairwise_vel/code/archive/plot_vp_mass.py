# Plot pairwise velocity
# v0.1.0, 2021-06-16 - Code Started

# Imports
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys


sim_box = sys.argv[1]

# bin_centres = np.arange(0.5, 150.5, 1.)
# mass_bin_lims = np.arange(10.75, 15.25, 0.25)
# colors = sns.color_palette("viridis", n_colors=17)

bin_edges = np.logspace(np.log10(0.01), np.log10(100.), 81)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
# mass_bin_lims = np.arange(10.5, 15.5, 0.5)
mass_bin_lims = np.arange(11.5, 14., 0.5)
colors = sns.color_palette("viridis", n_colors=mass_bin_lims.size-1)

plt.figure(figsize=(8., 6.))
plt.axvline(x=7., linestyle="--", color="k")
plt.axhline(y=0., linestyle="-", color="k")
for i in range(mass_bin_lims.size-1):
	mvps = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/pairwise_vel/{}/m200b_{:.2f}-{:.2f}_vp.dat".format(sim_box, mass_bin_lims[i], mass_bin_lims[i+1]))
	plt.plot(bin_centres, np.nan_to_num(mvps[:, 1]), color=colors[i])

# vps = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/{}_m200b_vp.dat".format(sim_box))
# plt.plot(bin_centres, vps[:, 1], color="k")
plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p$ [km/s]")
plt.xlim(0.01, 100.)
plt.xscale("log")
# plt.ylim(-50., 350.)
plt.savefig("../output/plots/{}_m200b_vp_mass_plot.jpg".format(sim_box))
