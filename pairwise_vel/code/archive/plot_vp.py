# Plot pairwise velocity
# v0.1.0, 2021-06-16 - Code Started

# Imports
import matplotlib.pyplot as plt
import numpy as np


vps = np.loadtxt("../../../projects/rrg-wperciva/mj3chapm/P2/aemulus_mocks/measurements/aemulus_z0.70_TestBox000-000_m200b_vp.dat")
bin_centres = np.arange(0.5, 150.5, 1.)

plt.figure(figsize=(8., 6.))
plt.plot(bin_centres, vps[:, 1])
plt.axvline(x=7.5, linestyle="--", color="k")
plt.xlabel(r"s [$h^{-1}$ Mpc]")
plt.ylabel(r"$v_p$ [km/s]")
plt.xlim(0., 150.)
plt.savefig("../output/plots/aemulus_z0.70_TestBox000-000_m200b_vp_plot.jpg")
