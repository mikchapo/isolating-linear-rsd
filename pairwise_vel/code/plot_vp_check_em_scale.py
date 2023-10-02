# Compare total pairwise velocities for various component scalings
# v0.1.0, 2021-11-08 - Code started, copied from plot_vp_halos_ic_z_check.py

# Imports
import matplotlib.pyplot as plt
import numpy as np

from vp_plot_funcs import *


cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_planck_products/"
                          "AbacusCosmos_1100box_planck_FoF_halos/"
                          "info/cosmo_params.dat")

vel_components = [{"name": "v-lin", "label": "Linear Component", "ls": "--"},
                  {"name": "v-nl", "label": "Non-linear Component", "ls": ":"}]
labels = ["Total Velocity", "Tot. Vel., linear scaling = 0.5",
          "Tot. Vel., non-linear scaling = 0.5"]

z = 0.7
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
mass_bin_lims = np.arange(12., 14.5, 0.5)

linestyles = ["--", "-.", ":"]

for i in range(mass_bin_lims.size-1):
    initialize_vp_plot()
    ymax = 100.
    ymin = -250.

    for j, component in enumerate(vel_components):
        vp_path = ("/home/mj3chapm/scratch/abacus/"
                   "AbacusCosmos_1100box_products/"
                   "AbacusCosmos_1100box_planck_products/"
                   "AbacusCosmos_1100box_planck_FoF_halos/z0.700/"
                   "FoF_ngc_nplt_test_{}_{:.2f}-{:.2f}_vp"
                   ".dat".format(component["name"], mass_bin_lims[i],
                                 mass_bin_lims[i+1]))
        plot_vp(vp_path, component["label"], color="k", marker="",
                linestyle=component["ls"])

    for j in range(3):
        vp_path = ("../output/data_products/check_emulator_scaling_vp/"
                   "{}_v-tot_{:.2f}-{:.2f}_vp"
                   ".dat".format(j, mass_bin_lims[i],
                                 mass_bin_lims[i+1]))
        plot_vp(vp_path, labels[j], color="C{}".format(j), marker="")

    plt.ylim(ymin, ymax)
    plt.title("{:.2f}-{:.2f} Halo Mass".format(mass_bin_lims[i],
                                               mass_bin_lims[i+1]))
    output_path = ("../output/plots/"
                   "check_emulator_scaling_{:.2f}-{:.2f}"
                   ".jpg".format(mass_bin_lims[i], mass_bin_lims[i+1]))
    finish_vp_plot(output_path)


# Change Log
