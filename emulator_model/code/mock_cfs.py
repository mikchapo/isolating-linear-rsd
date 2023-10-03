# Redshift error correlation function check
# v0.2.0, 2022-02-04 - Updated for efficiency, removed transition method, added
#                      smoothed option

# Import
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns
import sys

from mock_meas import mock_meas


# User input parameters
# simname - [AbacusCosmos_1100box]
# boxid - [00, 01, 02, 03,...|00-0, 00-1, 00-2,...]
# halotype - [rockstar|FoF]
# scaling_method - [combined|split]
# scale_variable - [gamma_f|gamma_l|gamma_n]
# load_dvs - [True|False]
# smoothed_vel - [True|False]
simname = sys.argv[1]
boxid = sys.argv[2]
halotype = sys.argv[3]
scaling_method = sys.argv[4]
scale_variable = sys.argv[5]
load_dvs = sys.argv[6]
smoothed_vel = sys.argv[7]
N_grid = int(sys.argv[8])
smooth_type = sys.argv[9]
R_smooth = float(sys.argv[10])

load_dvs = True if load_dvs == "True" else False
smoothed_vel = True if smoothed_vel == "True" else False

redshift = 0.70
Lbox = 1100.
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
             "{}_{}_{}_halos/z0.700" .format(simname, simname, boxid, simname,
                                             boxid, halotype))
Path("{}/em_model_test".format(path_root)).mkdir(parents=False, exist_ok=True)
if smoothed_vel:
    real_cat_path = ("{}/halo_lin_vel_{}_{}_{}_smoothed_v2"
                     ".dat".format(path_root, N_grid, smooth_type, R_smooth))
    run_name = "{}_{}_{}_{}_smoothed_v2".format(scaling_method, N_grid,
                                                smooth_type, R_smooth)
else:
    real_cat_path = ("{}/halo_lin_vel.dat".format(path_root))
    run_name = scaling_method

if scaling_method == "split":
    run_name += "_{}".format(scale_variable)

bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))
colors = sns.color_palette("viridis")

if load_dvs:
    ref_dv = np.loadtxt("{}/em_model_test/combined_v2_1.0_red_dv"
                        ".dat".format(path_root))
else:
    output_root = "{}/em_model_test/combined_v2_1.0".format(path_root)
    ref_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift, Lbox,
                                scaling_method=scaling_method,
                                gamma_l=1.0)
ref_color = colors.pop(2)

fig_width = 10.
fig, axes = plt.subplots(2, 3,
                         figsize=(fig_width, fig_width / 9.5 * 4.),
                         sharex=True,
                         gridspec_kw={'hspace': 0,
                                      'height_ratios': [2, 1]},
                         dpi=100)

gamma_vals = [0.84, 0.92, 1.08, 1.16]
# gamma_vals = [0.84, 1.16]
for i, gamma_val in enumerate(gamma_vals):
    if load_dvs:
        red_dv = np.loadtxt("{}/em_model_test/{}_{}_red_dv"
                            ".dat".format(path_root, run_name, gamma_val))

    else:
        output_root = ("{}/em_model_test/{}_{}".format(path_root, run_name,
                                                       gamma_val))
        if scale_variable == "gamma_l":
            red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                        Lbox, scaling_method=scaling_method,
                                        gamma_l=gamma_val)

        elif scale_variable == "gamma_n":
            red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                        Lbox, scaling_method=scaling_method,
                                        gamma_n=gamma_val)

        else:
            red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                        Lbox, scaling_method=scaling_method,
                                        gamma_l=gamma_val)

    axes[0, 0].plot(seps, seps**2.*red_dv[:9, 1], color=colors[i],
                    marker=".")
    axes[1, 0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                    color=colors[i], marker=".")
    axes[0, 1].plot(seps, seps**2.*red_dv[9:18, 1],
                    label="{}={}".format(r"$\{}$".format(scale_variable),
                                         gamma_val),
                    color=colors[i], marker=".")
    # axes[1, 1].plot(seps,
    #                 (red_dv[9:18, 1] - ref_dv[9:18, 1])/ref_dv[9:18, 1],
    #                 color=colors[i], marker=".")
    axes[1, 1].plot(seps,
                    (red_dv[9:18, 1] - ref_dv[9:18, 1]),
                    color=colors[i], marker=".")
    axes[0, 2].plot(seps, red_dv[18:, 1], color=colors[i], marker=".")
    axes[1, 2].plot(seps,
                    (red_dv[18:, 1] - ref_dv[18:, 1])/ref_dv[18:, 1],
                    color=colors[i], marker=".")

    if len(gamma_vals) / (i+1) == 2:
        axes[0, 0].plot(seps, seps**2.*ref_dv[:9, 1], color=ref_color,
                        marker=".")
        axes[0, 1].plot(seps, seps**2.*ref_dv[9:18, 1],
                        label="{}={}".format(r"$\{}$".format(scale_variable),
                                             1.0) + " (Base)",
                        color=ref_color,
                        marker=".")
        axes[0, 2].plot(seps, ref_dv[18:, 1], color=ref_color, marker=".")

axes[0, 0].set_ylabel(r"$s^2\xi_0$")

axes[1, 0].axhline(y=0., color='k', linestyle='-')
axes[1, 0].set_xscale("log")
axes[1, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
axes[1, 0].set_xlabel(r"$s[h^{-1}Mpc]$")

axes[0, 1].set_ylabel(r"$s^2\xi_2$")

axes[1, 1].axhline(y=0., color='k', linestyle='-')
axes[1, 1].set_xscale("log")
# axes[1, 1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")
axes[1, 1].set_ylabel(r"$\xi_2 - \xi_2^b$")
axes[1, 1].set_xlabel(r"$s[h^{-1}Mpc]$")

axes[0, 2].set_yscale("log")
axes[0, 2].set_ylabel(r"$w_p$")
axes[0, 1].legend()

axes[1, 2].axhline(y=0., color='k', linestyle='-')
axes[1, 2].set_xscale("log")
axes[1, 2].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
axes[1, 2].set_xlabel(r"$r_p[h^{-1}Mpc]$")

if smoothed_vel:
    axes[0, 0].axvline(x=R_smooth, color="k", linestyle="--")
    axes[1, 0].axvline(x=R_smooth, color="k", linestyle="--")
    axes[0, 1].axvline(x=R_smooth, color="k", linestyle="--")
    axes[1, 1].axvline(x=R_smooth, color="k", linestyle="--")
    axes[0, 2].axvline(x=R_smooth, color="k", linestyle="--")
    axes[1, 2].axvline(x=R_smooth, color="k", linestyle="--")

if scale_variable == "gamma_l":
    if smoothed_vel:
        fig.suptitle("Change in correlation functions with linear velocity "
                     "scaling for {} smoothed, R={}, N={}, linear "
                     "velocities".format(smooth_type, R_smooth, N_grid))

    else:
        fig.suptitle("Change in correlation functions with linear velocity "
                     "scaling for unsmoothed linear velocities")

elif scale_variable == "gamma_n":
    if smoothed_vel:
        fig.suptitle("Change in correlation functions with non-linear "
                     "velocity scaling for {} smoothed, R={}, N={}, linear "
                     "velocities".format(smooth_type, R_smooth, N_grid))

    else:
        fig.suptitle("Change in correlation functions with non-linear "
                     "velocity scaling for unsmoothed linear velocities")

plt.tight_layout()
output_root = ("{}/"
               "em_model_test/{}".format(path_root, run_name))
plt.savefig("{}_cf_check.jpg".format(output_root))

# Change Log
# v0.1.1, 2022-02-04 - Copied from RSD to P2 for testing new emulator models
# v0.1.0, 2021-05-11 - Code started
