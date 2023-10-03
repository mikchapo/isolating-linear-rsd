# Compare smoothed and unsmoothed velocity distributions
# v0.2.1, 2023-06-08 - Updated for thesis

# Imports
import matplotlib.pyplot as plt
import numpy as np


mass_split = True
thesis = True
N_rows = "All"

full_unsmoothed_cat = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                 "AbacusCosmos_1100box_products/"
                                 "AbacusCosmos_1100box_00_products/"
                                 "AbacusCosmos_1100box_00_rockstar_halos/"
                                 "z0.700/halo_lin_vel.dat")

full_smoothed_cat = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                               "AbacusCosmos_1100box_products/"
                               "AbacusCosmos_1100box_00_products/"
                               "AbacusCosmos_1100box_00_rockstar_halos/z0.700/"
                               "halo_lin_vel_1100_tophat_5.0_smoothed_v2.dat")

N_unsm = full_unsmoothed_cat.shape[0]
N_sm = full_smoothed_cat.shape[0]
print("N unsmoothed: {}, N smoothed: {}".format(N_unsm, N_sm))
mass_bin_lims = np.arange(12., 16., 0.5)
mass_bin_centres = np.arange(12.25, 15.75, 0.5)
vel_stats = np.empty((5, 5, mass_bin_lims.size - 1))
labels = ["Total Vel.", "Unsmoothed Linear", "Unsmoothed Non-linear",
          "Smoothed Linear", "Smoothed Non-linear"]
colors = ["k", "C0", "C0", "C1", "C1"]
linestyles = ["-", "-", "--", "-", "--"]

for i in range(mass_bin_lims.size - 1):
    print("Starting mass bin {:.2f}-{:.2f}".format(mass_bin_lims[i],
                                                   mass_bin_lims[i+1]))

    unsmoothed_cat = full_unsmoothed_cat[(full_unsmoothed_cat[:, 1] >=
                                          10**mass_bin_lims[i]) * (
                                          full_unsmoothed_cat[:, 1] <
                                          10**mass_bin_lims[i+1])]
    smoothed_cat = full_smoothed_cat[(full_smoothed_cat[:, 1] >=
                                      10**mass_bin_lims[i]) * (
                                      full_smoothed_cat[:, 1] <
                                      10**mass_bin_lims[i+1])]

    N_bin_unsm = unsmoothed_cat.shape[0]
    N_bin_sm = smoothed_cat.shape[0]
    print("N unsmoothed: {}, N smoothed: {}".format(N_bin_unsm, N_bin_sm))

    tot_mag = np.sqrt(unsmoothed_cat[:, 7]**2. +
                      unsmoothed_cat[:, 8]**2. +
                      unsmoothed_cat[:, 9]**2.)
    unsmoothed_lin_mag = np.sqrt(unsmoothed_cat[:, 10]**2. +
                                 unsmoothed_cat[:, 11]**2. +
                                 unsmoothed_cat[:, 12]**2.)
    unsmoothed_nl_mag = np.sqrt((unsmoothed_cat[:, 7] - unsmoothed_cat[:, 10])**2. +
                                (unsmoothed_cat[:, 8] - unsmoothed_cat[:, 11])**2. +
                                (unsmoothed_cat[:, 9] - unsmoothed_cat[:, 12])**2.)

    smoothed_lin_mag = np.sqrt(smoothed_cat[:, 10]**2. + smoothed_cat[:, 11]**2. +
                               smoothed_cat[:, 12]**2.)
    smoothed_nl_mag = np.sqrt((smoothed_cat[:, 7] - smoothed_cat[:, 10])**2. +
                              (smoothed_cat[:, 8] - smoothed_cat[:, 11])**2. +
                              (smoothed_cat[:, 9] - smoothed_cat[:, 12])**2.)

    vel_lists = [tot_mag, unsmoothed_lin_mag, unsmoothed_nl_mag,
                 smoothed_lin_mag, smoothed_nl_mag]

    for j in range(5):
        vel_stats[j, 0, i] = np.quantile(vel_lists[j], 0.5)
        vel_stats[j, 1, i] = np.quantile(vel_lists[j], 0.25)
        vel_stats[j, 2, i] = np.quantile(vel_lists[j], 0.75)
        vel_stats[j, 3, i] = np.mean(vel_lists[j])
        vel_stats[j, 4, i] = np.std(vel_lists[j]) / np.sqrt(N_bin_sm)


if thesis:
    width = 6.375
    aspect_ratio = 3. / 4.

else:
   width = 4.75
   aspect_ratio = 4. / 4.75

plt.figure(figsize=(width, width * aspect_ratio), dpi=300)
for j in range(5):
    plt.errorbar(mass_bin_lims[:-1]+0.15+0.05*j, vel_stats[j, 0, :],
                 yerr=vel_stats[j, 1:3, :], label=labels[j], color=colors[j],
                 linestyle=linestyles[j], marker="o", capsize=5)
# plt.xlabel(r"Halo Mass $[M_\odot/h]$")
plt.xlabel(r"$\log(M_H)$")
plt.ylabel("Velocity [km/s]")
plt.legend()

if thesis:
    plt.savefig("../output/plots/smooth_vel_check/vel_stats_median_v2_thesis.png")

else:
    plt.savefig("../output/plots/smooth_vel_check/vel_stats_median_v2.jpg")

plt.figure(figsize=(width, width * aspect_ratio), dpi=300)
for j in range(5):
    plt.errorbar(mass_bin_lims[:-1]+0.15+0.05*j, vel_stats[j, 3, :],
                 yerr=vel_stats[j, 4, :], label=labels[j], color=colors[j],
                 linestyle=linestyles[j], marker="o", capsize=5)
    # plt.plot(mass_bin_lims[:-1], vel_stats[j, 3, :],
    #          color=colors[j], linestyle=linestyles[j], label=labels[j],
    #          marker="o")
# plt.xlabel(r"Halo Mass $[M_\odot/h]$")
plt.xlabel(r"$\log(M_H)$")
plt.ylabel("Velocity [km/s]")
plt.legend()

if thesis:
    plt.savefig("../output/plots/smooth_vel_check/vel_stats_mean_v2_thesis.png")

else:
    plt.savefig("../output/plots/smooth_vel_check/vel_stats_mean_v2.jpg")

# Change Log
# v0.2.0, 2022-05-26 - Code updated to do all mass bins at once
# v0.1.0, 2022-05-13 - Code started from check_halo_type_vel_dist.py
