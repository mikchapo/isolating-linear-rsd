# Compare FoF and Rockstar velocity distributions
# v0.1.0, 2021-12-03 - Code started from check_emulator_scaling.py

# Imports
import matplotlib.pyplot as plt
import numpy as np


mass_split = True
N_rows = "All"

full_rockstar_cat = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                               "AbacusCosmos_1100box_products/"
                               "AbacusCosmos_1100box_planck_products/"
                               "AbacusCosmos_1100box_planck_rockstar_halos/"
                               "z0.700/halo_lin_vel.dat")

full_FoF_cat = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_planck_products/"
                          "AbacusCosmos_1100box_planck_FoF_halos/z0.700/"
                          "halo_lin_vel.dat")

mass_bin_lims = np.arange(12., 14.5, 0.5)
if mass_split:
    for_length = mass_bin_lims.size - 1
else:
    for_length = 1

for i in range(for_length):
    if mass_split:
        print("Starting mass bin {:.2f}-{:.2f}".format(mass_bin_lims[i],
                                                       mass_bin_lims[i+1]))
        # rockstar_cat = full_rockstar_cat[(full_rockstar_cat[:, 1] / 0.6726 >=
        #                                   10**mass_bin_lims[i]) * (
        #                                   full_rockstar_cat[:, 1] / 0.6726 <
        #                                   10**mass_bin_lims[i+1])]
        rockstar_cat = full_rockstar_cat[(full_rockstar_cat[:, 1] >=
                                          10**mass_bin_lims[i]) * (
                                          full_rockstar_cat[:, 1] <
                                          10**mass_bin_lims[i+1])]
        FoF_cat = full_FoF_cat[(full_FoF_cat[:, 1] >= 10**mass_bin_lims[i]) *
                               (full_FoF_cat[:, 1] < 10**mass_bin_lims[i+1])]
        print("Rockstar Size:", rockstar_cat.shape[0])
        print("FoF Size:", FoF_cat.shape[0])
        print("First 10 FoF:", full_FoF_cat[:10, 1])

    else:
        rockstar_cat = full_rockstar_cat
        FoF_cat = full_FoF_cat

    rockstar_tot_mag = np.sqrt(rockstar_cat[:, 7]**2. +
                               rockstar_cat[:, 8]**2. +
                               rockstar_cat[:, 9]**2.)
    rockstar_lin_mag = np.sqrt(rockstar_cat[:, 10]**2. +
                               rockstar_cat[:, 11]**2. +
                               rockstar_cat[:, 12]**2.)
    rockstar_nl_mag = np.sqrt((rockstar_cat[:, 7] - rockstar_cat[:, 10])**2. +
                              (rockstar_cat[:, 8] - rockstar_cat[:, 11])**2. +
                              (rockstar_cat[:, 9] - rockstar_cat[:, 12])**2.)

    FoF_tot_mag = np.sqrt(FoF_cat[:, 7]**2. + FoF_cat[:, 8]**2. +
                          FoF_cat[:, 9]**2.)
    FoF_lin_mag = np.sqrt(FoF_cat[:, 10]**2. + FoF_cat[:, 11]**2. +
                          FoF_cat[:, 12]**2.)
    FoF_nl_mag = np.sqrt((FoF_cat[:, 7] - FoF_cat[:, 10])**2. +
                         (FoF_cat[:, 8] - FoF_cat[:, 11])**2. +
                         (FoF_cat[:, 9] - FoF_cat[:, 12])**2.)

    mean_rockstar_tot_mag = np.mean(rockstar_tot_mag)
    mean_rockstar_lin_mag = np.mean(rockstar_lin_mag)
    mean_rockstar_nl_mag = np.mean(rockstar_nl_mag)
    mean_FoF_tot_mag = np.mean(FoF_tot_mag)
    mean_FoF_lin_mag = np.mean(FoF_lin_mag)
    mean_FoF_nl_mag = np.mean(FoF_nl_mag)

    # min_tot_mag = np.floor(np.log10(np.min((np.min(rockstar_tot_mag),
    #                        np.min(FoF_tot_mag)))))
    # max_tot_mag = np.ceil(np.log10(np.max((np.max(rockstar_tot_mag),
    #                       np.max(FoF_tot_mag)))))

    # min_lin_mag = np.floor(np.log10(np.min((np.min(rockstar_lin_mag),
    #                        np.min(FoF_lin_mag)))))
    # max_lin_mag = np.ceil(np.log10(np.max((np.max(rockstar_lin_mag),
    #                       np.max(FoF_lin_mag)))))
    min_tot_mag = np.min((np.min(rockstar_tot_mag),
                          np.min(FoF_tot_mag)))
    max_tot_mag = np.max((np.max(rockstar_tot_mag),
                          np.max(FoF_tot_mag)))

    min_lin_mag = np.min((np.min(rockstar_lin_mag),
                          np.min(FoF_lin_mag)))
    max_lin_mag = np.max((np.max(rockstar_lin_mag),
                          np.max(FoF_lin_mag)))

    min_nl_mag = np.min((np.min(rockstar_nl_mag),
                         np.min(FoF_nl_mag)))
    max_nl_mag = np.max((np.max(rockstar_nl_mag),
                         np.max(FoF_nl_mag)))

    plt.figure(figsize=(8., 6.), dpi=300)
    # plt.hist(rockstar_tot_mag, bins=np.logspace(min_tot_mag, max_tot_mag, 100),
    #          color="C0", label="Rockstar", alpha=0.7)
    # plt.hist(FoF_tot_mag, bins=np.logspace(min_tot_mag, max_tot_mag, 100),
    #          color="C1", label="FoF", alpha=0.7)
    plt.hist(rockstar_tot_mag, bins=np.linspace(min_tot_mag, max_tot_mag, 100),
             color="C0", label="Rockstar", alpha=0.7)
    plt.hist(FoF_tot_mag, bins=np.linspace(min_tot_mag, max_tot_mag, 100),
             color="C1", label="FoF", alpha=0.7)
    # plt.hist(rockstar_tot_mag, bins=100,
    #          color="C0", label="Rockstar", alpha=0.7)
    # plt.hist(FoF_tot_mag, bins=100,
    #          color="C1", label="FoF", alpha=0.7)
    plt.axvline(mean_rockstar_tot_mag, color="C0", linestyle="-",
                label="Rockstar mean={:.2f}".format(mean_rockstar_tot_mag))
    plt.axvline(mean_FoF_tot_mag, color="C1", linestyle="-",
                label="FoF mean={:.2f}".format(mean_FoF_tot_mag))
    plt.title("Total velocity magnitude")
    plt.xlabel(r"$v_{tot}$")
    plt.ylabel("Count")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()
    if mass_split:
        plt.savefig("../output/plots/halo_type_tot_mag_dist_{}"
                    "{:.2f}-{:.2f}_d3.jpg".format(N_rows, mass_bin_lims[i],
                                               mass_bin_lims[i+1]))
    else:
        plt.savefig("../output/plots/halo_type_tot_mag_dist_"
                    "{}_d3.jpg".format(N_rows))

    plt.figure(figsize=(8., 6.), dpi=300)
    # plt.hist(rockstar_lin_mag, bins=np.logspace(min_lin_mag, max_lin_mag, 100),
    #          color="C0", label="Rockstar", alpha=0.7)
    # plt.hist(FoF_lin_mag, bins=np.logspace(min_lin_mag, max_lin_mag, 100),
    #          color="C1", label="FoF", alpha=0.7)
    plt.hist(rockstar_lin_mag, bins=np.linspace(min_lin_mag, max_lin_mag, 100),
             color="C0", label="Rockstar", alpha=0.7)
    plt.hist(FoF_lin_mag, bins=np.linspace(min_lin_mag, max_lin_mag, 100),
             color="C1", label="FoF", alpha=0.7)
    # plt.hist(rockstar_lin_mag, bins=100,
    #          color="C0", label="Rockstar", alpha=0.7)
    # plt.hist(FoF_lin_mag, bins=100,
    #          color="C1", label="FoF", alpha=0.7)
    plt.axvline(mean_rockstar_lin_mag, color="C0", linestyle="-",
                label="Rockstar mean={:.2f}".format(mean_rockstar_lin_mag))
    plt.axvline(mean_FoF_lin_mag, color="C1", linestyle="-",
                label="FoF mean={:.2f}".format(mean_FoF_lin_mag))
    plt.title("Linear velocity magnitude")
    plt.xlabel(r"$v_{lin}$")
    plt.ylabel("Count")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()
    if mass_split:
        plt.savefig("../output/plots/halo_type_lin_mag_dist_{}_"
                    "{:.2f}-{:.2f}_d3.jpg".format(N_rows, mass_bin_lims[i],
                                               mass_bin_lims[i+1]))
    else:
        plt.savefig("../output/plots/halo_type_lin_mag_dist_"
                    "{}_d3.jpg".format(N_rows))

    plt.figure(figsize=(8., 6.), dpi=300)
    # plt.hist(rockstar_nl_mag, bins=np.logspace(min_nl_mag, max_nl_mag, 100),
    #          color="C0", label="Rockstar", alpha=0.7)
    # plt.hist(FoF_nl_mag, bins=np.logspace(min_nl_mag, max_nl_mag, 100),
    #          color="C1", label="FoF", alpha=0.7)
    plt.hist(rockstar_nl_mag, bins=np.linspace(min_nl_mag, max_nl_mag, 100),
             color="C0", label="Rockstar", alpha=0.7)
    plt.hist(FoF_nl_mag, bins=np.linspace(min_nl_mag, max_nl_mag, 100),
             color="C1", label="FoF", alpha=0.7)
    # plt.hist(rockstar_nl_mag, bins=100,
    #          color="C0", label="Rockstar", alpha=0.7)
    # plt.hist(FoF_nl_mag, bins=100,
    #          color="C1", label="FoF", alpha=0.7)
    plt.axvline(mean_rockstar_nl_mag, color="C0", linestyle="-",
                label="Rockstar mean={:.2f}".format(mean_rockstar_nl_mag))
    plt.axvline(mean_FoF_nl_mag, color="C1", linestyle="-",
                label="FoF mean={:.2f}".format(mean_FoF_nl_mag))
    plt.title("Non-linear velocity magnitude")
    plt.xlabel(r"$v_{nl}$")
    plt.ylabel("Count")
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend()
    if mass_split:
        plt.savefig("../output/plots/halo_type_nl_mag_dist_{}_"
                    "{:.2f}-{:.2f}_d3.jpg".format(N_rows, mass_bin_lims[i],
                                               mass_bin_lims[i+1]))
    else:
        plt.savefig("../output/plots/halo_type_nl_mag_dist_"
                    "{}_d3.jpg".format(N_rows))

# Change Log
