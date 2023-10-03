# Compare smoothed and unsmoothed velocity distributions
# v0.1.0, 2022-05-13 - Code started from check_halo_type_vel_dist.py

# Imports
import matplotlib.pyplot as plt
import numpy as np


mass_split = True
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
                               "halo_lin_vel_smoothed.dat")

mass_bin_lims = np.arange(12., 14.5, 0.5)
if mass_split:
    for_length = mass_bin_lims.size - 1
else:
    for_length = 1

for i in range(for_length):
    if mass_split:
        print("Starting mass bin {:.2f}-{:.2f}".format(mass_bin_lims[i],
                                                       mass_bin_lims[i+1]))
        # unsmoothed_cat = full_unsmoothed_cat[(full_unsmoothed_cat[:, 1] / 0.6726 >=
        #                                   10**mass_bin_lims[i]) * (
        #                                   full_unsmoothed_cat[:, 1] / 0.6726 <
        #                                   10**mass_bin_lims[i+1])]
        unsmoothed_cat = full_unsmoothed_cat[(full_unsmoothed_cat[:, 1] >=
                                              10**mass_bin_lims[i]) * (
                                              full_unsmoothed_cat[:, 1] <
                                              10**mass_bin_lims[i+1])]
        smoothed_cat = full_smoothed_cat[(full_smoothed_cat[:, 1] >=
                                          10**mass_bin_lims[i]) * (
                                          full_smoothed_cat[:, 1] <
                                          10**mass_bin_lims[i+1])]
        print("Unsmoothed Size:", unsmoothed_cat.shape[0])
        print("Smoothed Size:", smoothed_cat.shape[0])
        print("First 10 smoothed:", full_smoothed_cat[:10, 1])

    else:
        unsmoothed_cat = full_unsmoothed_cat
        smoothed_cat = full_smoothed_cat

    unsmoothed_tot_mag = np.sqrt(unsmoothed_cat[:, 7]**2. +
                                 unsmoothed_cat[:, 8]**2. +
                                 unsmoothed_cat[:, 9]**2.)
    unsmoothed_lin_mag = np.sqrt(unsmoothed_cat[:, 10]**2. +
                                 unsmoothed_cat[:, 11]**2. +
                                 unsmoothed_cat[:, 12]**2.)
    unsmoothed_nl_mag = np.sqrt((unsmoothed_cat[:, 7] - unsmoothed_cat[:, 10])**2. +
                                (unsmoothed_cat[:, 8] - unsmoothed_cat[:, 11])**2. +
                                (unsmoothed_cat[:, 9] - unsmoothed_cat[:, 12])**2.)

    smoothed_tot_mag = np.sqrt(smoothed_cat[:, 7]**2. + smoothed_cat[:, 8]**2. +
                               smoothed_cat[:, 9]**2.)
    smoothed_lin_mag = np.sqrt(smoothed_cat[:, 10]**2. + smoothed_cat[:, 11]**2. +
                               smoothed_cat[:, 12]**2.)
    smoothed_nl_mag = np.sqrt((smoothed_cat[:, 7] - smoothed_cat[:, 10])**2. +
                              (smoothed_cat[:, 8] - smoothed_cat[:, 11])**2. +
                              (smoothed_cat[:, 9] - smoothed_cat[:, 12])**2.)

    mean_unsmoothed_tot_mag = np.mean(unsmoothed_tot_mag)
    mean_unsmoothed_lin_mag = np.mean(unsmoothed_lin_mag)
    mean_unsmoothed_nl_mag = np.mean(unsmoothed_nl_mag)
    mean_smoothed_tot_mag = np.mean(smoothed_tot_mag)
    mean_smoothed_lin_mag = np.mean(smoothed_lin_mag)
    mean_smoothed_nl_mag = np.mean(smoothed_nl_mag)

    # min_tot_mag = np.floor(np.log10(np.min((np.min(unsmoothed_tot_mag),
    #                        np.min(smoothed_tot_mag)))))
    # max_tot_mag = np.ceil(np.log10(np.max((np.max(unsmoothed_tot_mag),
    #                       np.max(smoothed_tot_mag)))))

    # min_lin_mag = np.floor(np.log10(np.min((np.min(unsmoothed_lin_mag),
    #                        np.min(smoothed_lin_mag)))))
    # max_lin_mag = np.ceil(np.log10(np.max((np.max(unsmoothed_lin_mag),
    #                       np.max(smoothed_lin_mag)))))
    min_tot_mag = np.min((np.min(unsmoothed_tot_mag),
                          np.min(smoothed_tot_mag)))
    max_tot_mag = np.max((np.max(unsmoothed_tot_mag),
                          np.max(smoothed_tot_mag)))

    min_lin_mag = np.min((np.min(unsmoothed_lin_mag),
                          np.min(smoothed_lin_mag)))
    max_lin_mag = np.max((np.max(unsmoothed_lin_mag),
                          np.max(smoothed_lin_mag)))

    min_nl_mag = np.min((np.min(unsmoothed_nl_mag),
                         np.min(smoothed_nl_mag)))
    max_nl_mag = np.max((np.max(unsmoothed_nl_mag),
                         np.max(smoothed_nl_mag)))

    plt.figure(figsize=(8., 6.), dpi=300)
    # plt.hist(unsmoothed_tot_mag, bins=np.logspace(min_tot_mag, max_tot_mag, 100),
    #          color="C0", label="Unsmoothed", alpha=0.7)
    # plt.hist(smoothed_tot_mag, bins=np.logspace(min_tot_mag, max_tot_mag, 100),
    #          color="C1", label="Smoothed", alpha=0.7)
    plt.hist(unsmoothed_tot_mag, bins=np.linspace(min_tot_mag, max_tot_mag, 100),
             color="C0", label="Unsmoothed", alpha=0.7)
    plt.hist(smoothed_tot_mag, bins=np.linspace(min_tot_mag, max_tot_mag, 100),
             color="C1", label="Smoothed", alpha=0.7)
    # plt.hist(unsmoothed_tot_mag, bins=100,
    #          color="C0", label="Unsmoothed", alpha=0.7)
    # plt.hist(smoothed_tot_mag, bins=100,
    #          color="C1", label="Smoothed", alpha=0.7)
    plt.axvline(mean_unsmoothed_tot_mag, color="C0", linestyle="-",
                label="Unsmoothed mean={:.2f}".format(mean_unsmoothed_tot_mag))
    plt.axvline(mean_smoothed_tot_mag, color="C1", linestyle="-",
                label="Smoothed mean={:.2f}".format(mean_smoothed_tot_mag))
    plt.title("Total velocity magnitude")
    plt.xlabel(r"$v_{tot}$")
    plt.ylabel("Count")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()
    if mass_split:
        plt.savefig("../output/plots/smooth_vel_check/tot_mag_dist_{}_"
                    "{:.2f}-{:.2f}_d3.jpg".format(N_rows, mass_bin_lims[i],
                                                  mass_bin_lims[i+1]))
    else:
        plt.savefig("../output/plots/smooth_vel_check/tot_mag_dist_"
                    "{}_d3.jpg".format(N_rows))

    plt.figure(figsize=(8., 6.), dpi=300)
    # plt.hist(unsmoothed_lin_mag, bins=np.logspace(min_lin_mag, max_lin_mag, 100),
    #          color="C0", label="Unsmoothed", alpha=0.7)
    # plt.hist(smoothed_lin_mag, bins=np.logspace(min_lin_mag, max_lin_mag, 100),
    #          color="C1", label="Smoothed", alpha=0.7)
    plt.hist(unsmoothed_lin_mag, bins=np.linspace(min_lin_mag, max_lin_mag, 100),
             color="C0", label="Unsmoothed", alpha=0.7)
    plt.hist(smoothed_lin_mag, bins=np.linspace(min_lin_mag, max_lin_mag, 100),
             color="C1", label="Smoothed", alpha=0.7)
    # plt.hist(unsmoothed_lin_mag, bins=100,
    #          color="C0", label="Unsmoothed", alpha=0.7)
    # plt.hist(smoothed_lin_mag, bins=100,
    #          color="C1", label="Smoothed", alpha=0.7)
    plt.axvline(mean_unsmoothed_lin_mag, color="C0", linestyle="-",
                label="Unsmoothed mean={:.2f}".format(mean_unsmoothed_lin_mag))
    plt.axvline(mean_smoothed_lin_mag, color="C1", linestyle="-",
                label="Smoothed mean={:.2f}".format(mean_smoothed_lin_mag))
    plt.title("Linear velocity magnitude")
    plt.xlabel(r"$v_{lin}$")
    plt.ylabel("Count")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()
    if mass_split:
        plt.savefig("../output/plots/smooth_vel_check/lin_mag_dist_{}_"
                    "{:.2f}-{:.2f}_d3.jpg".format(N_rows, mass_bin_lims[i],
                                                  mass_bin_lims[i+1]))
    else:
        plt.savefig("../output/plots/smooth_vel_check/lin_mag_dist_"
                    "{}_d3.jpg".format(N_rows))

    plt.figure(figsize=(8., 6.), dpi=300)
    # plt.hist(unsmoothed_nl_mag, bins=np.logspace(min_nl_mag, max_nl_mag, 100),
    #          color="C0", label="Unsmoothed", alpha=0.7)
    # plt.hist(smoothed_nl_mag, bins=np.logspace(min_nl_mag, max_nl_mag, 100),
    #          color="C1", label="Smoothed", alpha=0.7)
    plt.hist(unsmoothed_nl_mag, bins=np.linspace(min_nl_mag, max_nl_mag, 100),
             color="C0", label="Unsmoothed", alpha=0.7)
    plt.hist(smoothed_nl_mag, bins=np.linspace(min_nl_mag, max_nl_mag, 100),
             color="C1", label="Smoothed", alpha=0.7)
    # plt.hist(unsmoothed_nl_mag, bins=100,
    #          color="C0", label="Unsmoothed", alpha=0.7)
    # plt.hist(smoothed_nl_mag, bins=100,
    #          color="C1", label="Smoothed", alpha=0.7)
    plt.axvline(mean_unsmoothed_nl_mag, color="C0", linestyle="-",
                label="Unsmoothed mean={:.2f}".format(mean_unsmoothed_nl_mag))
    plt.axvline(mean_smoothed_nl_mag, color="C1", linestyle="-",
                label="Smoothed mean={:.2f}".format(mean_smoothed_nl_mag))
    plt.title("Non-linear velocity magnitude")
    plt.xlabel(r"$v_{nl}$")
    plt.ylabel("Count")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()
    if mass_split:
        plt.savefig("../output/plots/smooth_vel_check/nl_mag_dist_{}_"
                    "{:.2f}-{:.2f}_d3.jpg".format(N_rows, mass_bin_lims[i],
                                               mass_bin_lims[i+1]))
    else:
        plt.savefig("../output/plots/smooth_vel_check/nl_mag_dist_"
                    "{}_d3.jpg".format(N_rows))

# Change Log