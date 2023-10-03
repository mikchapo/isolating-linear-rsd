# Check scaling parameters
# v0.1.0, 2021-10-05 - Code started

# Imports
import matplotlib.pyplot as plt
import numpy as np


N_rows = 50

orig_cat = np.loadtxt("../data/halo_lin_vel_ngc_nplt.dat")
zz_cat = np.loadtxt("../data/halo_new_velocity_test0", max_rows=N_rows)
lin_cat = np.loadtxt("../data/halo_new_velocity_test1", max_rows=N_rows)
nl_cat = np.loadtxt("../data/halo_new_velocity_test2", max_rows=N_rows)

for i in range(N_rows):
    no_scale_ratio = orig_cat[i, 7:10] / zz_cat[i, -3:]
    lin_scale_ratio = ((orig_cat[i, 7:10] - orig_cat[i, -3:]) +
                       orig_cat[i, -3:]*0.5) / lin_cat[i, -3:]
    nl_scale_ratio = ((orig_cat[i, 7:10] - orig_cat[i, -3:])*0.5 +
                      orig_cat[i, -3:]) / nl_cat[i, -3:]

    print()
    print("Galaxy {}".format(i))
    print(no_scale_ratio)
    print(lin_scale_ratio)
    print(nl_scale_ratio)

lin_mag_ratio = np.sqrt((orig_cat[:, 10]**2. + orig_cat[:, 11]**2. +
                        orig_cat[:, 12]**2.) /
                        (orig_cat[:, 7]**2. + orig_cat[:, 8]**2. +
                        orig_cat[:, 9]**2.))

mean_lin_mag_ratio = np.mean(lin_mag_ratio)

nl_mag_ratio = np.sqrt(((orig_cat[:, 7] - orig_cat[:, 10])**2. +
                        (orig_cat[:, 8] - orig_cat[:, 11])**2. +
                        (orig_cat[:, 9] - orig_cat[:, 12])**2.) /
                       (orig_cat[:, 7]**2. + orig_cat[:, 8]**2. +
                        orig_cat[:, 9]**2.))

mean_nl_mag_ratio = np.mean(nl_mag_ratio)

print("Mean linear component ratio to total:", mean_lin_mag_ratio)
print("Mean non-linear component ratio to total:", mean_nl_mag_ratio)

plt.figure(figsize=(8., 6.), dpi=300)
plt.hist(lin_mag_ratio, bins=np.logspace(np.log10(0.0001), np.log10(100.),
                                         100),
         color="C0")
plt.xlabel(r"$v_{lin} / v_{tot}$")
plt.ylabel("Count")
plt.xscale("log")
plt.yscale("log")
plt.savefig("../output/plots/lin_mag_ratio_hist.jpg")

plt.figure(figsize=(8., 6.), dpi=300)
plt.hist(nl_mag_ratio, bins=np.logspace(np.log10(0.0001), np.log10(100.), 100),
         color="C1")
plt.xlabel(r"$v_{nl} / v_{tot}$")
plt.ylabel("Count")
plt.xscale("log")
plt.yscale("log")
plt.savefig("../output/plots/nl_mag_ratio_hist.jpg")
