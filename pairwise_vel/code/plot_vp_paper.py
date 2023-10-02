# Plot halo pairwise velocity as a function of separation for paper
# v0.1.1, 2022-09-19 - Updated to be single figure, proportioned for paper

# Imports
import matplotlib.pyplot as plt
import numpy as np

from vp_plot_funcs import *


simname = "AbacusCosmos_1100box_planck"
boxid = "00-0"
single_path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                    "{}_{}_rockstar_halos".format(simname, simname, boxid,
                                                  simname, boxid))
cosmo_params = np.loadtxt("{}/info/cosmo_params.dat".format(single_path_root))
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "combined_pairwise_vel".format(simname))

vel_components = [{"name": "v-tot", "label": "Total Halo Velocity", "ls": "-"},
                  {"name": "v-lin", "label": "Linear Component", "ls": "--"},
                  {"name": "v-nl", "label": "Non-linear Component", "ls": ":"}]

z = 0.7
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
mass_bin_lims = np.arange(12., 15.5, 0.5)

fig, axes = initialize_paper_plot()
plot_static_paper(z, cosmo_params, axes)
for i in range(mass_bin_lims.size-1):
    m = i // 2
    n = i % 2
    for j, component in enumerate(vel_components):
        # vp_path = ("{}/z0.700/pairwise_vel/"
        #            "halo_{}_{:.2f}-{:.2f}_vp_v2_all-halos"
        #            ".dat".format(path_root, component["name"],
        #                          mass_bin_lims[i],
        #                          mass_bin_lims[i+1]))
        vp_path = ("{}/halo_{}_{:.2f}-{:.2f}_vp_v2_all-halos"
                   ".dat".format(path_root, component["name"],
                                 mass_bin_lims[i],
                                 mass_bin_lims[i+1]))
        plot_vp_paper(vp_path, component["label"], axes, m, n,
                      color="C0", marker="",
                      linestyle=component["ls"])

    for j, component in enumerate(vel_components[1:]):
        # vp_path = ("{}/z0.700/pairwise_vel/"
        #            "halo_{}_{:.2f}-{:.2f}_vp_v2_vel-sm"
        #            ".dat".format(path_root, component["name"],
        #                          mass_bin_lims[i],
        #                          mass_bin_lims[i+1]))
        vp_path = ("{}/halo_{}_{:.2f}-{:.2f}_vp_v2_vel-sm_all-halos"
                   ".dat".format(path_root, component["name"],
                                 mass_bin_lims[i],
                                 mass_bin_lims[i+1]))
        plot_vp_paper(vp_path, "Smoothed {}".format(component["label"]), axes,
                      m, n, color="C1", marker="",
                      linestyle=component["ls"])

    axes[m, n].set_title("{:.2f}-{:.2f} Halo Mass".format(mass_bin_lims[i],
                                                          mass_bin_lims[i+1]))

axes[0, 0].legend()
output_path = ("../output/plots/halo_vp_paper.jpg")

plot_lin_pred_paper(z, cosmo_params, axes)
finish_paper_plot(output_path)

# Change Log
# v0.1.0, 2022-08-09 - Code started, copied from plot_vp_check_halo_type.py
