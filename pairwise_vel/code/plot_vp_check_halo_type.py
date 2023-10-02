# Compare pairwise velocities for FoF and Rockstar halos
# v0.1.0, 2021-11-11 - Code started, copied from plot_vp_check_em_scale.py

# Imports
import matplotlib.pyplot as plt
import numpy as np

from vp_plot_funcs import *


cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_planck_products/"
                          "AbacusCosmos_1100box_planck_rockstar_halos/"
                          "info/cosmo_params.dat")

halo_types = ["FoF", "rockstar"]
vel_components = [{"name": "v-tot", "label": "Total Halo Velocity", "ls": "-"},
                  {"name": "v-lin", "label": "Linear Component", "ls": "--"},
                  {"name": "v-nl", "label": "Non-linear Component", "ls": ":"}]
labels = ["Total Velocity", "Tot. Vel., linear scaling = 0.5",
          "Tot. Vel., non-linear scaling = 0.5"]

z = 0.7
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
mass_bin_lims = np.arange(12., 14.5, 0.5)


for i in range(mass_bin_lims.size-1):
    initialize_vp_plot()
    ymax = 650.
    ymin = -850.

    plot_static(z, cosmo_params)

    for j, halo_type in enumerate(halo_types):
        for k, component in enumerate(vel_components):
            vp_path = ("/home/mj3chapm/scratch/abacus/"
                       "AbacusCosmos_1100box_products/"
                       "AbacusCosmos_1100box_planck_products/"
                       "AbacusCosmos_1100box_planck_{}_halos/z0.700/"
                       "halo_{}_{:.2f}-{:.2f}_vp_v2"
                       ".dat".format(halo_type, component["name"],
                                     mass_bin_lims[i],
                                     mass_bin_lims[i+1]))
            plot_vp(vp_path, "{} {}".format(halo_type, component["label"]),
                    color="C{}".format(j), marker="",
                    linestyle=component["ls"])

    plot_lin_pred(z, cosmo_params)

    plt.ylim(ymin, ymax)
    plt.title("{:.2f}-{:.2f} Halo Mass".format(mass_bin_lims[i],
                                               mass_bin_lims[i+1]))
    output_path = ("../output/plots/"
                   "check_halo_type_{:.2f}-{:.2f}_static_v2"
                   ".jpg".format(mass_bin_lims[i], mass_bin_lims[i+1]))
    finish_vp_plot(output_path)


# Change Log
