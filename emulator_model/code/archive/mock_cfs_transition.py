# Redshift error correlation function check
# v0.1.1, 2022-02-04 - Copied from RSD to P2 for testing new emulator models

# Import
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

from mock_meas import mock_meas


# User input parameters
# halotype - [rockstar|FoF]
# scaling_method - [combined|split|transition]
# scale_variable - [gamma_f|gamma_l|gamma_n]
# b - transition width
# c - transition location
# load_dvs - [True|False]
halotype = sys.argv[1]
scaling_method = sys.argv[2]
scale_variable = sys.argv[3]
tb = float(sys.argv[4])
tc = float(sys.argv[5])
load_dvs = sys.argv[6]

load_dvs = True if load_dvs == "True" else False
if scaling_method == "combined":
    run_name = scaling_method
elif scaling_method == "split":
    run_name = "{}_{}".format(scaling_method, scale_variable)
elif scaling_method == "transition":
    run_name = "{}_{}_b{}_c{}".format(scaling_method, scale_variable, tb, tc)

redshift = 0.70
Lbox = 1100.
real_cat_path = ("/home/mj3chapm/scratch/abacus/"
                 "AbacusCosmos_1100box_planck_products/"
                 "AbacusCosmos_1100box_planck_00-0_products/"
                 "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
                 "halo_lin_vel.dat".format(halotype))

bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))
colors = sns.color_palette("viridis")

# if load_dvs:
ref_dv = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                    "AbacusCosmos_1100box_planck_products/"
                    "AbacusCosmos_1100box_planck_00-0_products/"
                    "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
                    "em_model_test/planck_00-0_{}_halos"
                    "_combined_1.0_red_dv.dat".format(halotype, halotype))
# else:
#     output_root = ("/home/mj3chapm/scratch/abacus/"
#                    "AbacusCosmos_1100box_planck_products/"
#                    "AbacusCosmos_1100box_planck_00-0_products/"
#                    "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
#                    "em_model_test/planck_00-0_{}_halos"
#                    "_{}_1.0".format(halotype, halotype, run_name))
#     ref_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift, Lbox,
#                                 scaling_method=scaling_method,
#                                 gamma_l=1.0, b=b, c=c)
ref_color = colors.pop(2)

fig_width = 10.
fig, axes = plt.subplots(2, 3,
                         figsize=(fig_width, fig_width / 9.5 * 4.),
                         sharex=True,
                         gridspec_kw={'hspace': 0,
                                      'height_ratios': [2, 1]},
                         dpi=100)

if load_dvs:
    gamma_vals = [0.84, 0.92, 1.08, 1.16]
else:
    gamma_vals = [0.84, 1.16]
for i, gamma_val in enumerate(gamma_vals):
    output_root = ("/home/mj3chapm/scratch/abacus/"
                   "AbacusCosmos_1100box_planck_products/"
                   "AbacusCosmos_1100box_planck_00-0_products/"
                   "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
                   "em_model_test/planck_00-0_{}_halos"
                   "_{}_{}".format(halotype, halotype, run_name,
                                   gamma_val))
    if load_dvs:
        red_dv = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                            "AbacusCosmos_1100box_planck_products/"
                            "AbacusCosmos_1100box_planck_00-0_products/"
                            "AbacusCosmos_1100box_planck_00-0_{}_halos/"
                            "z0.700/em_model_test/planck_00-0_{}_halos"
                            "_{}_{}_red_dv.dat".format(halotype, halotype,
                                                       run_name,
                                                       gamma_val))

    else:
        if scale_variable == "gamma_l":
            red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                        Lbox, scaling_method=scaling_method,
                                        gamma_l=gamma_val, tb=tb, tc=tc)

        elif scale_variable == "gamma_n":
            red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                        Lbox, scaling_method=scaling_method,
                                        gamma_n=gamma_val, tb=tb, tc=tc)

        else:
            red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                        Lbox, scaling_method=scaling_method,
                                        gamma_l=gamma_val, tb=tb, tc=tc)

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

plt.tight_layout()
output_root = ("/home/mj3chapm/scratch/abacus/"
               "AbacusCosmos_1100box_planck_products/"
               "AbacusCosmos_1100box_planck_00-0_products/"
               "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
               "em_model_test/planck_00-0_{}_halos"
               "_{}".format(halotype, halotype, run_name))
plt.savefig("{}_cf_check.jpg".format(output_root))

# Change Log
# v0.1.0, 2021-05-11 - Code started
