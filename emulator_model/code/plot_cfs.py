"""
Plot gamma effects on halo correlation functions.

v0.2.2, 2023-06-08 - Updated paper plot to add thesis option
"""

# Import
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))


def plot_cfs(axes, red_dv, ref_dv, scale_variable, gamma_val, color,
             direct_comp=True, linestyle="-"):
    """Plot the correlation function."""
    if direct_comp:
        axes[0, 0].plot(seps, seps**2.*red_dv[:9, 1], color=color,
                        marker=".", linestyle=linestyle)
        axes[0, 1].plot(seps, seps**2.*red_dv[9:18, 1],
                        label="{}={}".format(r"$\{}$".format(scale_variable),
                                             gamma_val),
                        color=color, marker=".", linestyle=linestyle)
        axes[0, 2].plot(seps, red_dv[18:, 1], color=color, marker=".",
                        linestyle=linestyle)
    # axes[1, 1].plot(seps,
    #                 (red_dv[9:18, 1] - ref_dv[9:18, 1])/ref_dv[9:18, 1],
    #                 color=color, marker=".")
    axes[1, 0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                    color=color, marker=".", linestyle=linestyle)
    axes[1, 1].plot(seps,
                    (red_dv[9:18, 1] - ref_dv[9:18, 1]),
                    color=color, marker=".", linestyle=linestyle)
    axes[1, 2].plot(seps,
                    (red_dv[18:, 1] - ref_dv[18:, 1])/ref_dv[18:, 1],
                    color=color, marker=".", linestyle=linestyle)


def plot_cfs_relative(axes, red_dv, ref_dv, scale_variable, gamma_val, color,
                      linestyle="-", line_label=None):
    """Plot the realtive correlation function."""
    # axes[1, 1].plot(seps,
    #                 (red_dv[9:18, 1] - ref_dv[9:18, 1])/ref_dv[9:18, 1],
    #                 color=color, marker=".")
    if line_label is not None:
        axes[0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                     color=color, marker=".", linestyle=linestyle,
                     label=line_label)
    else:
        axes[0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                     color=color, marker=".", linestyle=linestyle)
    if linestyle == "-":
        axes[1].plot(seps, (red_dv[9:18, 1] - ref_dv[9:18, 1]),
                     color=color, marker=".", linestyle=linestyle,
                     label="{}={}".format(r"$\{}$".format(scale_variable),
                                          gamma_val))
    else:
        axes[1].plot(seps, (red_dv[9:18, 1] - ref_dv[9:18, 1]),
                     color=color, marker=".", linestyle=linestyle)
    axes[2].plot(seps, (red_dv[18:, 1] - ref_dv[18:, 1])/ref_dv[18:, 1],
                 color=color, marker=".", linestyle=linestyle)


def plot_cfs_relative_st(axes, red_dv, ref_dv, color, linestyle="-",
                         line_label=None):
    """Plot relative correlation functions."""
    # axes[1, 1].plot(seps,
    #                 (red_dv[9:18, 1] - ref_dv[9:18, 1])/ref_dv[9:18, 1],
    #                 color=color, marker=".")
    axes[0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                 color=color, marker=".", linestyle=linestyle)
    axes[1].plot(seps, (red_dv[9:18, 1] - ref_dv[9:18, 1]),
                 color=color, marker=".", linestyle=linestyle)
    if line_label is not None:
        axes[2].plot(seps, (red_dv[18:, 1] - ref_dv[18:, 1])/ref_dv[18:, 1],
                     color=color, marker=".", linestyle=linestyle,
                     label=line_label)
    else:
        axes[2].plot(seps, (red_dv[18:, 1] - ref_dv[18:, 1])/ref_dv[18:, 1],
                     color=color, marker=".", linestyle=linestyle)


def plot_cfs_relative_paper(axes, red_dv, ref_dv, color, linestyle="-",
                            line_label=None):
    """Plot relative correlation functions for the paper plot."""
    # axes[1, 1].plot(seps,
    #                 (red_dv[9:18, 1] - ref_dv[9:18, 1])/ref_dv[9:18, 1],
    #                 color=color, marker=".")
    if line_label is not None:
        axes[0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                     color=color, marker=".", linestyle=linestyle,
                     label=line_label)
    else:
        axes[0].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                     color=color, marker=".", linestyle=linestyle)
    axes[1].plot(seps, (red_dv[9:18, 1] - ref_dv[9:18, 1]),
                 color=color, marker=".", linestyle=linestyle)
    # axes[2].plot(seps, (red_dv[18:, 1] - ref_dv[18:, 1])/ref_dv[18:, 1],
    #              color=color, marker=".", linestyle=linestyle)


def plot_mono(axes, red_dv, ref_dv, scale_variable, gamma_val, color,
              direct_comp=True, linestyle="-"):
    """Plot the monopole."""
    if direct_comp:
        axes[0].plot(seps, seps**2.*red_dv[:9, 1], color=color,
                     marker=".", linestyle=linestyle,
                     label="{}={}".format(r"$\{}$".format(scale_variable),
                                          gamma_val))

    axes[1].plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                 color=color, marker=".", linestyle=linestyle)


def plot_mono_rel(red_dv, ref_dv, scale_variable, gamma_val, color,
                  linestyle="-", line_label=None):
    """Plot the relative difference in the monopole."""
    if line_label is not None:
        plt.plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                 color=color, marker=".", linestyle=linestyle,
                 label=line_label)
    else:
        plt.plot(seps, (red_dv[:9, 1] - ref_dv[:9, 1])/ref_dv[:9, 1],
                 color=color, marker=".", linestyle=linestyle)


def plot_cf_comp_halotype(scaling_method, scale_variable, tb=0.2, tc=0.5):
    """Compare the change in correlation functions for different halo type."""
    if scaling_method == "combined":
        run_name = scaling_method
    elif scaling_method == "split":
        run_name = "{}_{}".format(scaling_method, scale_variable)
    elif scaling_method == "transition":
        run_name = "{}_{}_b{}_c{}".format(scaling_method, scale_variable, tb,
                                          tc)

    fig_width = 10.
    fig, axes = plt.subplots(2, 3,
                             figsize=(fig_width, fig_width / 9.5 * 4.),
                             sharex=True,
                             gridspec_kw={'hspace': 0,
                                          'height_ratios': [2, 1]},
                             dpi=100)

    if scaling_method == "combined" or scaling_method == "split":
        gamma_vals = [0.84, 0.92, 1.08, 1.16]
        halotypes = ["rockstar", "FoF"]
        linestyles = ["-", "--"]
        include_direct = [True, False]
        ref_dvs = np.empty((2, 27, 2))
        for j, halotype in enumerate(halotypes):
            ref_dv = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                "AbacusCosmos_1100box_planck_products/"
                                "AbacusCosmos_1100box_planck_00-0_products/"
                                "AbacusCosmos_1100box_planck_00-0_{}_halos/"
                                "z0.700/em_model_test/planck_00-0_{}_halos"
                                "_combined_1.0_red_dv.dat".format(halotype,
                                                                  halotype))
            ref_dvs[j, :, :] = ref_dv
    else:
        gamma_vals = [0.84, 1.16]
        halotype = "rockstar"
        ref_dv = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                            "AbacusCosmos_1100box_planck_products/"
                            "AbacusCosmos_1100box_planck_00-0_products/"
                            "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
                            "em_model_test/planck_00-0_{}_halos"
                            "_combined_1.0_red_dv.dat".format(halotype,
                                                              halotype))

    colors = sns.color_palette("viridis")
    ref_color = colors.pop(int(len(gamma_vals) / 2))

    for i, gamma_val in enumerate(gamma_vals):
        if scaling_method == "combined" or scaling_method == "split":
            for j, halotype in enumerate(halotypes):
                red_dv = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                    "AbacusCosmos_1100box_planck_products/"
                                    "AbacusCosmos_1100box_planck"
                                    "_00-0_products/"
                                    "AbacusCosmos_1100box_planck"
                                    "_00-0_{}_halos/"
                                    "z0.700/em_model_test/planck_00-0_{}_halos"
                                    "_{}_{}_red_dv.dat".format(halotype,
                                                               halotype,
                                                               run_name,
                                                               gamma_val))

                plot_cfs(axes, red_dv, ref_dvs[j, :, :], scale_variable,
                         gamma_val, colors[i],
                         direct_comp=include_direct[j],
                         linestyle=linestyles[j])

        else:
            halotype = "rockstar"
            red_dv = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                "AbacusCosmos_1100box_planck_products/"
                                "AbacusCosmos_1100box_planck"
                                "_00-0_products/"
                                "AbacusCosmos_1100box_planck"
                                "_00-0_{}_halos/"
                                "z0.700/em_model_test/planck_00-0_{}_halos"
                                "_{}_{}_red_dv.dat".format(halotype,
                                                           halotype,
                                                           run_name,
                                                           gamma_val))

            plot_cfs(axes, red_dv, ref_dv, scale_variable, gamma_val,
                     colors[i], direct_comp=include_direct[j],
                     linestyle=linestyles[j])

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
    halotype = "rockstar"
    output_root = ("/home/mj3chapm/scratch/abacus/"
                   "AbacusCosmos_1100box_planck_products/"
                   "AbacusCosmos_1100box_planck_00-0_products/"
                   "AbacusCosmos_1100box_planck_00-0_{}_halos/z0.700/"
                   "em_model_test/planck_00-0_{}_halos"
                   "_{}".format(halotype, halotype, run_name))
    plt.savefig("{}_cf_check_v2.jpg".format(output_root))


def plot_cf_smooth_comp(scale_variable, scaling_method="split"):
    """Compare the change in correlation function from smoothing."""
    simname = "AbacusCosmos_1100box"
    boxid = "00"
    halotype = "rockstar"
    run_name = "{}_{}".format(scaling_method, scale_variable)

    fig_width = 10.
    fig, axes = plt.subplots(2, 3,
                             figsize=(fig_width, fig_width / 9.5 * 4.),
                             sharex=True,
                             gridspec_kw={'hspace': 0,
                                          'height_ratios': [2, 1]},
                             dpi=100)

    gamma_vals = [0.84, 0.92, 1.08, 1.16]
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_{}_halos/z0.700/em_model_test".format(simname, simname,
                                                              boxid, simname,
                                                              boxid, halotype))
    ref_dv = np.loadtxt("{}/combined_1.0_red_dv.dat".format(path_root))

    colors = sns.color_palette("viridis")
    ref_color = colors.pop(int(len(gamma_vals) / 2))

    for i, gamma_val in enumerate(gamma_vals):
        red_dv = np.loadtxt("{}/{}_{}_{}_red_dv"
                            ".dat".format(path_root, scaling_method,
                                          scale_variable, gamma_val))
        red_dv_smoothed = np.loadtxt("{}/{}_smoothed_{}_{}_red_dv"
                                     ".dat".format(path_root, scaling_method,
                                                   scale_variable, gamma_val))

        plot_cfs(axes, red_dv_smoothed, ref_dv, scale_variable,
                 gamma_val, colors[i],
                 direct_comp=True,
                 linestyle="-")
        plot_cfs(axes, red_dv, ref_dv, scale_variable,
                 gamma_val, colors[i],
                 direct_comp=False,
                 linestyle="--")

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
    halotype = "rockstar"
    output_root = ("{}/{}_{}".format(path_root, scaling_method,
                                     scale_variable))
    plt.savefig("{}_sm_usm_cf_check.jpg".format(output_root))


def plot_cf_smooth_comp_rel(scale_variable, scaling_method="split"):
    """Compare realtive change in correlation function for smoothing."""
    simname = "AbacusCosmos_1100box"
    boxid = "00"
    halotype = "rockstar"
    run_name = "{}_{}".format(scaling_method, scale_variable)

    fig_width = 15.
    fig, axes = plt.subplots(1, 3,
                             figsize=(fig_width, fig_width / 9.5 * 3.),
                             dpi=100)

    # gamma_vals = [0.84, 0.92, 1.08, 1.16]
    gamma_vals = [0.84, 1.16]
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_{}_halos/z0.700/em_model_test".format(simname, simname,
                                                              boxid, simname,
                                                              boxid, halotype))
    ref_dv = np.loadtxt("{}/combined_1.0_red_dv.dat".format(path_root))

    colors = sns.color_palette("viridis")
    ref_color = colors.pop(int(len(gamma_vals) / 2))
    colors.pop(1)
    colors.pop(1)

    for i, gamma_val in enumerate(gamma_vals):
        red_dv = np.loadtxt("{}/{}_{}_{}_red_dv"
                            ".dat".format(path_root, scaling_method,
                                          scale_variable, gamma_val))
        red_dv_smoothed = np.loadtxt("{}/{}_smoothed_{}_{}_red_dv"
                                     ".dat".format(path_root, scaling_method,
                                                   scale_variable, gamma_val))

        if i == 0:
            plot_cfs_relative(axes, red_dv_smoothed, ref_dv, scale_variable,
                              gamma_val, colors[i],
                              linestyle="--", line_label="Smoothed")
            plot_cfs_relative(axes, red_dv, ref_dv, scale_variable,
                              gamma_val, colors[i],
                              linestyle="-", line_label="Unsmoothed")
        else:
            plot_cfs_relative(axes, red_dv_smoothed, ref_dv, scale_variable,
                              gamma_val, colors[i],
                              linestyle="--")
            plot_cfs_relative(axes, red_dv, ref_dv, scale_variable,
                              gamma_val, colors[i],
                              linestyle="-")

    axes[0].axhline(y=0., color='k', linestyle='-')
    axes[0].axvline(x=5., color='k', linestyle='--')
    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes[0].set_xlabel(r"$s[h^{-1}Mpc]$")
    axes[0].legend()

    axes[1].axhline(y=0., color='k', linestyle='-')
    axes[1].axvline(x=5., color='k', linestyle='--')
    axes[1].set_xscale("log")
    # axes[1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")
    axes[1].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes[1].set_xlabel(r"$s[h^{-1}Mpc]$")
    axes[1].legend()

    axes[2].axhline(y=0., color='k', linestyle='-')
    axes[2].axvline(x=5., color='k', linestyle='--')
    axes[2].set_xscale("log")
    axes[2].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
    axes[2].set_xlabel(r"$r_p[h^{-1}Mpc]$")

    if scale_variable == "gamma_l":
        fig.suptitle("Relative change in correlation functions with linear "
                     "velocity scaling")
    elif scale_variable == "gamma_n":
        fig.suptitle("Relative change in correlation functions with "
                     "non-linear velocity scaling")

    plt.tight_layout()
    halotype = "rockstar"
    output_root = ("{}/{}_{}".format(path_root, scaling_method,
                                     scale_variable))
    plt.savefig("{}_sm_usm_cf_check_rel.jpg".format(output_root))


def plot_mono_smooth_comp(scale_variable, scaling_method="split"):
    """Compare change in correaltion function for monopole."""
    simname = "AbacusCosmos_1100box"
    boxid = "00"
    halotype = "rockstar"
    run_name = "{}_{}".format(scaling_method, scale_variable)

    fig, axes = plt.subplots(2, 1,
                             figsize=(4., 4.),
                             sharex=True,
                             gridspec_kw={'hspace': 0,
                                          'height_ratios': [2, 1]},
                             dpi=100)

    gamma_vals = [0.84, 0.92, 1.08, 1.16]
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_{}_halos/z0.700/em_model_test".format(simname, simname,
                                                              boxid, simname,
                                                              boxid, halotype))
    ref_dv = np.loadtxt("{}/combined_1.0_red_dv.dat".format(path_root))

    colors = sns.color_palette("viridis")
    ref_color = colors.pop(int(len(gamma_vals) / 2))

    for i, gamma_val in enumerate(gamma_vals):
        red_dv_smoothed = np.loadtxt("{}/{}_smoothed_{}_{}_red_dv"
                                     ".dat".format(path_root, scaling_method,
                                                   scale_variable, gamma_val))

        plot_mono(axes, red_dv_smoothed, ref_dv, scale_variable,
                  gamma_val, colors[i],
                  direct_comp=True,
                  linestyle="-")

        if len(gamma_vals) / (i+1) == 2:
            axes[0].plot(seps, seps**2.*ref_dv[:9, 1], color=ref_color,
                         label="{}={}".format(r"$\{}$".format(scale_variable),
                                              1.0) + " (Base)",
                         marker=".")

    axes[0].set_ylabel(r"$s^2\xi_0$")
    axes[0].legend()

    axes[1].axhline(y=0., color='k', linestyle='-')
    axes[1].set_xscale("log")
    axes[1].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes[1].set_xlabel(r"$s[h^{-1}Mpc]$")

    plt.tight_layout()
    halotype = "rockstar"
    output_root = ("{}/{}_{}".format(path_root, scaling_method, scale_variable))
    plt.savefig("{}_sm_usm_mono_check.jpg".format(output_root))

    fig = plt.figure(figsize=(4., 3.), dpi=100)

    gamma_vals = [0.84, 1.16]
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_{}_halos/z0.700/em_model_test".format(simname, simname,
                                                              boxid, simname,
                                                              boxid, halotype))
    ref_dv = np.loadtxt("{}/combined_1.0_red_dv.dat".format(path_root))

    colors = sns.color_palette("viridis")
    ref_color = colors.pop(int(len(gamma_vals) / 2))
    colors.pop(1)
    colors.pop(1)

    for i, gamma_val in enumerate(gamma_vals):
        red_dv = np.loadtxt("{}/{}_{}_{}_red_dv"
                            ".dat".format(path_root, scaling_method,
                                          scale_variable, gamma_val))
        red_dv_smoothed = np.loadtxt("{}/{}_smoothed_{}_{}_red_dv"
                                     ".dat".format(path_root, scaling_method,
                                                   scale_variable, gamma_val))

        if i == 0:
            plot_mono_rel(red_dv_smoothed, ref_dv, scale_variable,
                          gamma_val, colors[i],
                          linestyle="--", line_label="Smoothed")
            plot_mono_rel(red_dv, ref_dv, scale_variable,
                          gamma_val, colors[i],
                          linestyle="-", line_label="Unsmoothed")
        else:
            plot_mono_rel(red_dv_smoothed, ref_dv, scale_variable,
                          gamma_val, colors[i],
                          linestyle="--")
            plot_mono_rel(red_dv, ref_dv, scale_variable,
                          gamma_val, colors[i],
                          linestyle="-")

    plt.axhline(y=0., color='k', linestyle='-')
    plt.axvline(x=5., color='k', linestyle='--')
    plt.xscale("log")
    plt.ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    plt.xlabel(r"$s[h^{-1}Mpc]$")
    plt.legend()

    plt.tight_layout()
    halotype = "rockstar"
    output_root = ("{}/{}_{}".format(path_root, scaling_method,
                                     scale_variable))
    plt.savefig("{}_sm_usm_mono_check_rel.jpg".format(output_root))


def plot_cf_smooth_test_comp_rel(scale_variable, N_grid=1100,
                                 smooth_type="tophat", R_smooth=5.,
                                 scaling_method="split"):
    """Compare different smoothing methods using the correlation function."""
    simname = "AbacusCosmos_1100box"
    boxid = "00"
    halotype = "rockstar"

    fig_width = 10.
    fig, axes = plt.subplots(1, 3,
                             figsize=(fig_width, fig_width / 9.5 * 3.),
                             dpi=100)

    # gamma_vals = [0.84, 0.92, 1.08, 1.16]
    # gamma_vals = [0.84, 1.16]
    gamma_val = 1.16
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_{}_halos/z0.700/em_model_test".format(simname, simname,
                                                              boxid, simname,
                                                              boxid, halotype))
    ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv.dat".format(path_root))

    colors = sns.color_palette("viridis")

    if type(N_grid) is list:
        for i in range(len(N_grid)):
            run_name = "{}_{}_{}_{}_smoothed_v2_{}".format(scaling_method,
                                                           N_grid[i],
                                                           smooth_type,
                                                           R_smooth,
                                                           scale_variable)
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name, gamma_val))
            plot_cfs_relative_st(axes, red_dv, ref_dv, "C{}".format(i),
                                 linestyle="-",
                                 line_label=(r"$N_g$=" + str(N_grid[i])))
        axes[0].axvline(x=R_smooth, color='k', linestyle='--')
        axes[1].axvline(x=R_smooth, color='k', linestyle='--')
        axes[2].axvline(x=R_smooth, color='k', linestyle='--')

        fig.suptitle("Variable grid spacing, {} smoothed, R={}, {}=1.16 "
                     "rel. C.F. change".format(smooth_type, R_smooth,
                                               scale_variable))
        output_root = ("{}/variable_N_grid_{}_{}_{}_v2".format(path_root,
                                                               smooth_type,
                                                               R_smooth,
                                                               scale_variable))

    elif type(smooth_type) is list:
        for i in range(len(smooth_type)):
            run_name = ("{}_{}_{}_{}_smoothed_v2_"
                        "{}".format(scaling_method, N_grid, smooth_type[i],
                                    R_smooth[i], scale_variable))
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name, gamma_val))
            plot_cfs_relative_st(axes, red_dv, ref_dv, "C{}".format(i),
                                 linestyle="-",
                                 line_label="{}, R={}".format(smooth_type[i],
                                                              R_smooth[i]))
            axes[0].axvline(x=R_smooth[i], color='C{}'.format(i),
                            linestyle='--')
            axes[1].axvline(x=R_smooth[i], color='C{}'.format(i),
                            linestyle='--')
            axes[2].axvline(x=R_smooth[i], color='C{}'.format(i),
                            linestyle='--')
        fig.suptitle("N_g={}, kernel type comparrison, {}=1.16 "
                     "rel. C.F. change".format(N_grid, scale_variable))
        output_root = ("{}/{}_variable_type_smooth_{}"
                       "_v2".format(path_root, N_grid, scale_variable))

    elif type(R_smooth) is list:
        for i in range(len(R_smooth)):
            run_name = ("{}_{}_{}_{}_smoothed_v2_"
                        "{}".format(scaling_method, N_grid, smooth_type,
                                    R_smooth[i], scale_variable))
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name, gamma_val))
            plot_cfs_relative_st(axes, red_dv, ref_dv, "C{}".format(i),
                                 linestyle="-",
                                 line_label="R={}".format(R_smooth[i]))
            axes[0].axvline(x=R_smooth[i], color='C{}'.format(i),
                            linestyle='--')
            axes[1].axvline(x=R_smooth[i], color='C{}'.format(i),
                            linestyle='--')
            axes[2].axvline(x=R_smooth[i], color='C{}'.format(i),
                            linestyle='--')
        fig.suptitle("N_g={}, {} smoothed, variable smoothing scale, {}=0.84 "
                     "rel. C.F. change".format(N_grid, smooth_type,
                                               scale_variable))
        output_root = ("{}/{}_{}_variable_R_smooth_{}"
                       "_v2".format(path_root, N_grid, smooth_type,
                                    scale_variable))

    axes[0].axhline(y=0., color='k', linestyle='-')
    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes[0].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[1].axhline(y=0., color='k', linestyle='-')
    axes[1].set_xscale("log")
    # axes[1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")
    axes[1].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes[1].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[2].axhline(y=0., color='k', linestyle='-')
    axes[2].set_xscale("log")
    axes[2].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
    axes[2].set_xlabel(r"$r_p[h^{-1}Mpc]$")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig("{}_cf_check_rel.jpg".format(output_root))


def plot_cf_smooth_version_comp_rel(scale_variable, N_grid=1100,
                                    smooth_type="tophat", R_smooth=5.,
                                    scaling_method="split"):
    """Compare different smoothing methods."""
    simname = "AbacusCosmos_1100box"
    boxid = "00"
    halotype = "rockstar"

    fig_width = 10.
    fig, axes = plt.subplots(1, 3,
                             figsize=(fig_width, fig_width / 9.5 * 3.),
                             dpi=100)

    gamma_val = 0.84
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_{}_halos/z0.700/em_model_test".format(simname, simname,
                                                              boxid, simname,
                                                              boxid, halotype))

    ref_dv_v1 = np.loadtxt("{}/v1/combined_1.0_red_dv.dat".format(path_root))
    run_name_v1 = "{}_{}_{}_{}_smoothed_{}".format(scaling_method,
                                                   N_grid,
                                                   smooth_type, R_smooth,
                                                   scale_variable)
    red_dv_v1 = np.loadtxt("{}/v1/{}_{}_red_dv"
                           ".dat".format(path_root, run_name_v1, gamma_val))
    plot_cfs_relative_st(axes, red_dv_v1, ref_dv_v1, "C0",
                         linestyle="-",
                         line_label="v1")

    ref_dv_v2 = np.loadtxt("{}/combined_v2_1.0_red_dv.dat".format(path_root))
    run_name_v2 = "{}_{}_{}_{}_smoothed_v2_{}_v2".format(scaling_method,
                                                         N_grid,
                                                         smooth_type, R_smooth,
                                                         scale_variable)
    red_dv_v2 = np.loadtxt("{}/{}_{}_red_dv"
                           ".dat".format(path_root, run_name_v2, gamma_val))
    plot_cfs_relative_st(axes, red_dv_v2, ref_dv_v2, "C1",
                         linestyle="-",
                         line_label="v2")

    fig.suptitle("Version test, N_grid={}, {} smoothed, R={}, {}=0.84 "
                 "rel. C.F. change".format(N_grid, smooth_type, R_smooth,
                                           scale_variable))
    output_root = ("{}/{}_{}_{}_{}_version_test".format(path_root,
                                                        N_grid,
                                                        smooth_type,
                                                        R_smooth,
                                                        scale_variable))

    axes[0].axvline(x=R_smooth, color='k', linestyle='--')
    axes[1].axvline(x=R_smooth, color='k', linestyle='--')
    axes[2].axvline(x=R_smooth, color='k', linestyle='--')

    axes[0].axhline(y=0., color='k', linestyle='-')
    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes[0].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[1].axhline(y=0., color='k', linestyle='-')
    axes[1].set_xscale("log")
    # axes[1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")
    axes[1].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes[1].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[2].axhline(y=0., color='k', linestyle='-')
    axes[2].set_xscale("log")
    axes[2].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
    axes[2].set_xlabel(r"$r_p[h^{-1}Mpc]$")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig("{}_cf_check_rel.jpg".format(output_root))


def paper_comparison_single_plots():
    """Compare correaltion functions as separate plots for the paper."""
    simname = "AbacusCosmos_1100box_planck"
    boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
              "00-5", "00-6", "00-7", "00-8", "00-9",
              "00-10", "00-11", "00-12", "00-13", "00-14",
              "00-15", "00-16", "00-17", "00-18", "00-19"]
    halotype = "rockstar"
    scale_variables = ["gamma_l", "gamma_n"]
    N_grids = [550, 1375]
    N_grid_def = 1100
    smooth_type_def = "tophat"
    smooth_types = ["tophat", "gaussian"]
    R_smooths = [3., 7.]
    R_smooth_def = 5.
    scaling_method = "split"

    fig_width = 10.
    fig_list = []
    axes_list = []
    # for i in range(7):
    #     fig, axes = plt.subplots(2, 2,
    #                              figsize=(fig_width, fig_width * 0.85),
    #                              dpi=100, sharex=True)
    #     fig_list.append(fig)
    #     axes_list.append(axes)

    fig, axes1 = plt.subplots(2, 2,
                              figsize=(fig_width, fig_width * 0.85),
                              dpi=100, sharex=True)
    fig2, axes2 = plt.subplots(2, 2,
                               figsize=(fig_width, fig_width * 0.85),
                               dpi=100, sharex=True)
    fig3, axes3 = plt.subplots(2, 2,
                               figsize=(fig_width, fig_width * 0.85),
                               dpi=100, sharex=True)
    fig4, axes4 = plt.subplots(2, 2,
                               figsize=(fig_width, fig_width * 0.85),
                               dpi=100, sharex=True)
    fig5, axes5 = plt.subplots(2, 2,
                               figsize=(fig_width, fig_width * 0.85),
                               dpi=100, sharex=True)
    fig6, axes6 = plt.subplots(2, 2,
                               figsize=(fig_width, fig_width * 0.85),
                               dpi=100, sharex=True)
    fig7, axes7 = plt.subplots(2, 2,
                               figsize=(fig_width, fig_width * 0.85),
                               dpi=100, sharex=True)

    gamma_val = 1.2

    colors = sns.color_palette("viridis")
    print(len(colors))
    removed_color = colors.pop(1)
    ref_color = colors.pop(1)
    removed_color = colors.pop(1)
    linestyles = ["--", ":"]

    for i, scale_variable in enumerate(scale_variables):
        mono_comps = np.empty((20, 9))
        quad_comps = np.empty((20, 9))
        for j, boxid in enumerate(boxids):
            path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                         "{}_{}_products/{}_{}_{}_halos/z0.700/"
                         "em_model_test".format(simname, simname, boxid,
                                                simname, boxid, halotype))
            ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                ".dat".format(path_root))
            run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                        "{}".format(scaling_method, N_grid_def,
                                    smooth_type_def,
                                    R_smooth_def, scale_variable))
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name, gamma_val))

            mono_comps[j, :] = (red_dv[:9, 1] - ref_dv[:9, 1]) / ref_dv[:9, 1]
            quad_comps[j, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

        axes1[0, i].plot(seps, np.mean(mono_comps, axis=0),
                         color=ref_color, marker=".", linestyle="-",
                         label=(r"Tophat, $R_{TH}=5.0$, "
                                r"$L_g=1.0 h^{-1}{\rm Mpc}$ "
                                "(Fiducial)"))
        axes1[1, i].plot(seps, np.mean(quad_comps, axis=0),
                         color=ref_color, marker=".", linestyle="-")

        for j in range(20):
            axes2[0, i].plot(seps, mono_comps[j, :], alpha=0.3,
                             color=ref_color, linestyle="-")
            axes2[1, i].plot(seps, quad_comps[j, :], alpha=0.3,
                             color=ref_color, linestyle="-")

        axes2[0, i].errorbar(seps, np.mean(mono_comps, axis=0),
                             yerr=np.std(mono_comps, axis=0),
                             linewidth=3,
                             color=ref_color, linestyle="-",
                             label=(r"Tophat, $R_{TH}=5.0$, "
                                    r"$L_g=1.0 h^{-1}{\rm Mpc}$ "
                                    "(Fiducial)"))
        axes2[1, i].errorbar(seps, np.mean(quad_comps, axis=0),
                             yerr=np.std(quad_comps, axis=0),
                             linewidth=3,
                             color=ref_color, linestyle="-")

        for j in range(len(R_smooths)):
            mono_comps = np.empty((20, 9))
            quad_comps = np.empty((20, 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                            "{}".format(scaling_method, N_grid_def,
                                        smooth_type_def,
                                        R_smooths[j], scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = (red_dv[:9, 1] -
                                    ref_dv[:9, 1]) / ref_dv[:9, 1]
                quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            axes1[0, i].plot(seps, np.mean(mono_comps, axis=0),
                             color=colors[j], marker=".", linestyle="-",
                             label=(r"$R_{TH}$=" +
                                    str(R_smooths[j])))
            axes1[1, i].plot(seps, np.mean(quad_comps, axis=0),
                             color=colors[j], marker=".", linestyle="-")

            if j == 0:
                for k in range(20):
                    axes3[0, i].plot(seps, mono_comps[k, :], alpha=0.3,
                                     color=colors[j], linestyle="-")
                    axes3[1, i].plot(seps, quad_comps[k, :], alpha=0.3,
                                     color=colors[j], linestyle="-")

                axes3[0, i].errorbar(seps, np.mean(mono_comps, axis=0),
                                     yerr=np.std(mono_comps, axis=0),
                                     linewidth=3,
                                     color=colors[j], linestyle="-",
                                     label=(r"$R_{TH}$=" +
                                            str(R_smooths[j])))
                axes3[1, i].errorbar(seps, np.mean(quad_comps, axis=0),
                                     yerr=np.std(quad_comps, axis=0),
                                     linewidth=3,
                                     color=colors[j], linestyle="-")

            else:
                for k in range(20):
                    axes4[0, i].plot(seps, mono_comps[k, :], alpha=0.3,
                                     color=colors[j], linestyle="-")
                    axes4[1, i].plot(seps, quad_comps[k, :], alpha=0.3,
                                     color=colors[j], linestyle="-")

                axes4[0, i].errorbar(seps, np.mean(mono_comps, axis=0),
                                     yerr=np.std(mono_comps, axis=0),
                                     linewidth=3,
                                     color=colors[j], linestyle="-",
                                     label=(r"$R_{TH}$=" +
                                            str(R_smooths[j])))
                axes4[1, i].errorbar(seps, np.mean(quad_comps, axis=0),
                                     yerr=np.std(quad_comps, axis=0),
                                     linewidth=3,
                                     color=colors[j], linestyle="-")

        for j in range(len(N_grids)):
            grid_length = "{:.1f}".format(1100 / N_grids[j])
            mono_comps = np.empty((20, 9))
            quad_comps = np.empty((20, 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                            "{}".format(scaling_method, N_grids[j],
                                        smooth_type_def, R_smooth_def,
                                        scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = (red_dv[:9, 1] -
                                    ref_dv[:9, 1]) / ref_dv[:9, 1]
                quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            axes1[0, i].plot(seps, np.mean(mono_comps, axis=0),
                             color=ref_color, marker=".",
                             linestyle=linestyles[j],
                             label=(r"$L_g$=" + grid_length +
                                    r" $h^{-1}{\rm Mpc}$"))
            axes1[1, i].plot(seps, np.mean(quad_comps, axis=0),
                             color=ref_color, marker=".",
                             linestyle=linestyles[j])

            if j == 0:
                for k in range(20):
                    axes5[0, i].plot(seps, mono_comps[k, :], alpha=0.3,
                                     color=ref_color, linestyle="-")
                    axes5[1, i].plot(seps, quad_comps[k, :], alpha=0.3,
                                     color=ref_color, linestyle="-")

                axes5[0, i].errorbar(seps, np.mean(mono_comps, axis=0),
                                     yerr=np.std(mono_comps, axis=0),
                                     linewidth=3,
                                     color=ref_color, linestyle="-",
                                     label=(r"$L_g$=" + grid_length +
                                            r" $h^{-1}{\rm Mpc}$"))
                axes5[1, i].errorbar(seps, np.mean(quad_comps, axis=0),
                                     yerr=np.std(quad_comps, axis=0),
                                     linewidth=3,
                                     color=ref_color, linestyle="-")

            else:
                for k in range(20):
                    axes6[0, i].plot(seps, mono_comps[k, :], alpha=0.3,
                                     color=ref_color, linestyle="-")
                    axes6[1, i].plot(seps, quad_comps[k, :], alpha=0.3,
                                     color=ref_color, linestyle="-")

                axes6[0, i].errorbar(seps, np.mean(mono_comps, axis=0),
                                     yerr=np.std(mono_comps, axis=0),
                                     linewidth=3,
                                     color=ref_color, linestyle="-",
                                     label=(r"$L_g$=" + grid_length +
                                            r" $h^{-1}{\rm Mpc}$"))
                axes6[1, i].errorbar(seps, np.mean(quad_comps, axis=0),
                                     yerr=np.std(quad_comps, axis=0),
                                     linewidth=3,
                                     color=ref_color, linestyle="-")

        mono_comps = np.empty((20, 9))
        quad_comps = np.empty((20, 9))
        for j, boxid in enumerate(boxids):
            path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                         "{}_{}_products/{}_{}_{}_halos/z0.700/"
                         "em_model_test".format(simname, simname, boxid,
                                                simname, boxid, halotype))
            ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                ".dat".format(path_root))
            run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                        "{}".format(scaling_method, N_grid_def,
                                    "gaussian",
                                    2., scale_variable))
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name, gamma_val))

            mono_comps[j, :] = (red_dv[:9, 1] - ref_dv[:9, 1]) / ref_dv[:9, 1]
            quad_comps[j, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

        axes1[0, i].plot(seps, np.mean(mono_comps, axis=0),
                         color="r", marker=".", linestyle="-",
                         label=r"Gaussian, $\sigma$=2.0")
        axes1[1, i].plot(seps, np.mean(quad_comps, axis=0),
                         color="r", marker=".", linestyle="-")

        for j in range(20):
            axes7[0, i].plot(seps, mono_comps[j, :], alpha=0.3,
                             color="r", linestyle="-")
            axes7[1, i].plot(seps, quad_comps[j, :], alpha=0.3,
                             color="r", linestyle="-")

        axes7[0, i].errorbar(seps, np.mean(mono_comps, axis=0),
                             yerr=np.std(mono_comps, axis=0),
                             linewidth=3,
                             color="r", linestyle="-",
                             label=r"Gaussian, $\sigma$=2.0")
        axes7[1, i].errorbar(seps, np.mean(quad_comps, axis=0),
                             yerr=np.std(quad_comps, axis=0),
                             linewidth=3,
                             color="r", linestyle="-")

        axes1[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes1[0, i].axhline(y=0., color='k', linestyle='-')
        axes1[1, i].axhline(y=0., color='k', linestyle='-')
        # axes1[2, i].axhline(y=0., color='k', linestyle='-')
        axes1[1, i].set_xscale("log")
        axes1[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes1[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes1[0, i].set_ylim(-0.2, 0.2)
        axes1[1, i].set_ylim(-0.4, 0.4)
        # axes1[2, 0].set_ylim(-0.1, 0.1)

        axes2[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes2[0, i].axhline(y=0., color='k', linestyle='-')
        axes2[1, i].axhline(y=0., color='k', linestyle='-')
        # axes2[2, i].axhline(y=0., color='k', linestyle='-')
        axes2[1, i].set_xscale("log")
        axes2[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes2[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes2[0, i].set_ylim(-0.2, 0.2)
        axes2[1, i].set_ylim(-0.4, 0.4)
        # axes2[2, 0].set_ylim(-0.1, 0.1)

        axes3[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes3[0, i].axhline(y=0., color='k', linestyle='-')
        axes3[1, i].axhline(y=0., color='k', linestyle='-')
        # axes3[2, i].axhline(y=0., color='k', linestyle='-')
        axes3[1, i].set_xscale("log")
        axes3[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes3[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes3[0, i].set_ylim(-0.2, 0.2)
        axes3[1, i].set_ylim(-0.4, 0.4)
        # axes3[2, 0].set_ylim(-0.1, 0.1)

        axes4[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes4[0, i].axhline(y=0., color='k', linestyle='-')
        axes4[1, i].axhline(y=0., color='k', linestyle='-')
        # axes4[2, i].axhline(y=0., color='k', linestyle='-')
        axes4[1, i].set_xscale("log")
        axes4[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes4[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes4[0, i].set_ylim(-0.2, 0.2)
        axes4[1, i].set_ylim(-0.4, 0.4)
        # axes4[2, 0].set_ylim(-0.1, 0.1)

        axes5[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes5[0, i].axhline(y=0., color='k', linestyle='-')
        axes5[1, i].axhline(y=0., color='k', linestyle='-')
        # axes5[2, i].axhline(y=0., color='k', linestyle='-')
        axes5[1, i].set_xscale("log")
        axes5[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes5[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes5[0, i].set_ylim(-0.2, 0.2)
        axes5[1, i].set_ylim(-0.4, 0.4)
        # axes5[2, 0].set_ylim(-0.1, 0.1)

        axes6[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes6[0, i].axhline(y=0., color='k', linestyle='-')
        axes6[1, i].axhline(y=0., color='k', linestyle='-')
        # axes6[2, i].axhline(y=0., color='k', linestyle='-')
        axes6[1, i].set_xscale("log")
        axes6[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes6[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes6[0, i].set_ylim(-0.2, 0.2)
        axes6[1, i].set_ylim(-0.4, 0.4)
        # axes6[2, 0].set_ylim(-0.1, 0.1)

        axes6[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes6[0, i].axhline(y=0., color='k', linestyle='-')
        axes6[1, i].axhline(y=0., color='k', linestyle='-')
        # axes6[2, i].axhline(y=0., color='k', linestyle='-')
        axes6[1, i].set_xscale("log")
        axes6[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes6[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes6[0, i].set_ylim(-0.2, 0.2)
        axes6[1, i].set_ylim(-0.4, 0.4)
        # axes6[2, 0].set_ylim(-0.1, 0.1)

        axes7[0, i].set_title(r"$\{}={}$".format(scale_variable, gamma_val))
        axes7[0, i].axhline(y=0., color='k', linestyle='-')
        axes7[1, i].axhline(y=0., color='k', linestyle='-')
        # axes7[2, i].axhline(y=0., color='k', linestyle='-')
        axes7[1, i].set_xscale("log")
        axes7[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        # axes7[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes7[0, i].set_ylim(-0.2, 0.2)
        axes7[1, i].set_ylim(-0.4, 0.4)
        # axes7[2, 0].set_ylim(-0.1, 0.1)

    plot_names = ["all", "fiducial", "R-3", "R-7", "N-550", "N-1375",
                  "gaussian"]

    axes1[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes1[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes1[0, 0].legend()

    plt.figure(1)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[0]))

    axes2[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes2[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes2[0, 0].legend()

    plt.figure(2)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[1]))

    axes3[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes3[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes3[0, 0].legend()

    plt.figure(3)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[2]))

    axes4[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes4[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes4[0, 0].legend()

    plt.figure(4)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[3]))

    axes5[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes5[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes5[0, 0].legend()

    plt.figure(5)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[4]))

    axes6[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes6[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes6[0, 0].legend()

    plt.figure(6)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[5]))

    axes7[0, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes7[1, 0].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes7[0, 0].legend()

    plt.figure(7)
    plt.tight_layout()
    output_root = ("{}/paper_comp".format(path_root))
    plt.savefig("{}_cf_check_rel_{}.jpg".format(output_root,
                                                plot_names[6]))


def paper_comparison(thesis=False):
    """Test different smoothing methods using relative CF changes."""
    simname = "AbacusCosmos_1100box_planck"
    # boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
    #           "00-5", "00-6", "00-7", "00-8", "00-9",
    #           "00-10", "00-11", "00-12", "00-13", "00-14",
    #           "00-15", "00-16", "00-17", "00-18", "00-19"]
    boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
              "00-5", "00-6", "00-7", "00-8", "00-9",
              "00-10", "00-11", "00-12", "00-14",
              "00-15", "00-16", "00-17", "00-18", "00-19"]
    halotype = "rockstar"
    scale_variables = ["gamma_l", "gamma_n"]
    N_grids = [550, 1375]
    N_grid_def = 1100
    smooth_type_def = "tophat"
    smooth_types = ["tophat", "gaussian"]
    R_smooths = [3., 7.]
    R_smooth_def = 5.
    scaling_method = "split"

    if thesis:
        fig_width = 6.375
        fig, axes = plt.subplots(1, 2, figsize=(fig_width, fig_width * 0.9),
                                 dpi=300)

    else:
        fig_width = 10.
    # fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_width * 0.85),
    #                          dpi=300, sharex=True)
        fig, axes = plt.subplots(1, 2, figsize=(fig_width, fig_width * 0.65),
                                 dpi=300)

    gamma_vals = [1.2, 0.8]
    linestyles = ["-", "--"]

    for i, scale_variable in enumerate(scale_variables):
        for j, gamma_val in enumerate(gamma_vals):
            mono_comps = np.empty((len(boxids), 9))
            # quad_comps = np.empty((len(boxids), 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_all-halos_{}".format(scaling_method,
                                                     scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = ((red_dv[:9, 1] - ref_dv[:9, 1]) /
                                    ref_dv[:9, 1])
                # quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            if i == 0 and j == 0:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="k", marker=".",
                                linestyle=linestyles[j],
                                label="Unsmoothed")
            elif i == 1:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="k", marker=".",
                                linestyle=linestyles[j],
                                label=(r"$\gamma=$" + str(gamma_val)))
            else:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="k", marker=".",
                                linestyle=linestyles[j])
            # axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
            #                 color="k", marker=".",
            #                 linestyle=linestyles[j])

            mono_comps = np.empty((len(boxids), 9))
            # quad_comps = np.empty((len(boxids), 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                            "{}".format(scaling_method, N_grid_def,
                                        smooth_type_def,
                                        R_smooth_def, scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = ((red_dv[:9, 1] - ref_dv[:9, 1]) /
                                    ref_dv[:9, 1])
                # quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            if i == 0 and j == 0:
                # axes[i].plot(seps, np.mean(mono_comps, axis=0),
                #                 color="C0", marker=".",
                #                 linestyle=linestyles[j],
                #                 label=(r"Tophat, $R_{TH}=5.0$, "
                #                        r"$L_g=1.0 h^{-1}{\rm Mpc}$ "
                #                        "(Fiducial)"))
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C0", marker=".",
                                linestyle=linestyles[j], linewidth=4.,
                                label="Default Smoothed")

            else:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C0", marker=".",
                                linestyle=linestyles[j], linewidth=4.)
            # axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
            #                 color="C0", marker=".",
            #                 linestyle=linestyles[j], linewidth=4.)

            for m in range(len(R_smooths)):
                mono_comps = np.empty((len(boxids), 9))
                # quad_comps = np.empty((len(boxids), 9))
                for k, boxid in enumerate(boxids):
                    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                                 "{}_{}_products/{}_{}_{}_halos/z0.700/"
                                 "em_model_test".format(simname, simname,
                                                        boxid, simname, boxid,
                                                        halotype))
                    ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                        ".dat".format(path_root))
                    run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                                "{}".format(scaling_method, N_grid_def,
                                            smooth_type_def,
                                            R_smooths[m], scale_variable))
                    red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                        ".dat".format(path_root, run_name,
                                                      gamma_val))

                    mono_comps[k, :] = (red_dv[:9, 1] -
                                        ref_dv[:9, 1]) / ref_dv[:9, 1]
                    # quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

                if i == 0 and j == 0:
                    axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                    color="C{}".format(m+1), marker=".",
                                    linestyle=linestyles[j], alpha=0.6,
                                    label=(r"$R_{TH}$=" +
                                           str(R_smooths[m]) +
                                           r" $h^{-1}{\rm Mpc}$"))
                else:
                    axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                    color="C{}".format(m+1), marker=".",
                                    linestyle=linestyles[j], alpha=0.6,)

                # axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
                #                 color="C{}".format(m+1), marker=".",
                #                 linestyle=linestyles[j], alpha=0.6,)

            for m in range(len(N_grids)):
                grid_length = "{:.1f}".format(1100 / N_grids[m])
                mono_comps = np.empty((len(boxids), 9))
                # quad_comps = np.empty((len(boxids), 9))
                for k, boxid in enumerate(boxids):
                    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                                 "{}_{}_products/{}_{}_{}_halos/z0.700/"
                                 "em_model_test".format(simname, simname,
                                                        boxid, simname, boxid,
                                                        halotype))
                    ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                        ".dat".format(path_root))
                    run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                                "{}".format(scaling_method, N_grids[m],
                                            smooth_type_def, R_smooth_def,
                                            scale_variable))
                    red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                        ".dat".format(path_root, run_name,
                                                      gamma_val))

                    mono_comps[k, :] = (red_dv[:9, 1] -
                                        ref_dv[:9, 1]) / ref_dv[:9, 1]
                    # quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

                if i == 0 and j == 0:
                    axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                    color="C{}".format(m+3), marker=".",
                                    linestyle=linestyles[j], alpha=0.6,
                                    label=(r"$L_g$=" + grid_length +
                                           r" $h^{-1}{\rm Mpc}$"))
                else:
                    axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                    color="C{}".format(m+3), marker=".",
                                    linestyle=linestyles[j], alpha=0.6,)

                # axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
                #                 color="C{}".format(m+3), marker=".",
                #                 linestyle=linestyles[j], alpha=0.6,)

            mono_comps = np.empty((len(boxids), 9))
            # quad_comps = np.empty((len(boxids), 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                            "{}".format(scaling_method, N_grid_def,
                                        "gaussian",
                                        2., scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = ((red_dv[:9, 1] - ref_dv[:9, 1]) /
                                    ref_dv[:9, 1])
                # quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            if i == 0 and j == 0 and not thesis:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C5", marker=".",
                                linestyle=linestyles[j], alpha=0.6,
                                label=(r"Gaussian, $\sigma=2.0$ "
                                       r"$h^{-1}{\rm Mpc}$"))

            elif i == 0 and j == 0:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C5", marker=".",
                                linestyle=linestyles[j], alpha=0.6,
                                label=(r"Gauss. $\sigma=2.0$ "
                                       r"$h^{-1}{\rm Mpc}$"))

            else:
                axes[i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C5", marker=".",
                                linestyle=linestyles[j], alpha=0.6,)

            # axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
            #                 color="C5", marker=".",
            #                 linestyle=linestyles[j], alpha=0.6,)

        axes[i].set_title(r"Scaling $\{}$".format(scale_variable))
        axes[i].axhline(y=0., color='k', linestyle='-')
        axes[i].legend()
        # axes[1, i].axhline(y=0., color='k', linestyle='-')
        # axes[2, i].axhline(y=0., color='k', linestyle='-')
        axes[i].set_xscale("log")
        axes[i].set_xlabel(r"$s\,[h^{-1}Mpc]$")
        axes[i].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
        # axes[1, i].set_ylabel(r"$\xi_2 - \xi_2^b$")

        # axes[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        # axes[2, 0].set_ylim(-0.1, 0.1)

    # axes[0].set_ylim(-0.1, 0.1)
    axes[0].set_ylim(-0.2, 0.25)
    axes[1].set_ylim(-0.2, 0.25)
    # axes[1, 0].set_ylim(-0.08, 0.08)
    # axes[1, 1].set_ylim(-0.6, 0.6)

    plt.tight_layout()
    if thesis:
        output_root = ("{}/thesis_comp".format(path_root))
        plt.savefig("{}_cf_check_rel.png".format(output_root))
    else:
        output_root = ("{}/paper_comp".format(path_root))
        plt.savefig("{}_cf_check_rel.jpg".format(output_root))


def presentation_comparison():
    """Test different smoothing methods using relative CF changes."""
    simname = "AbacusCosmos_1100box_planck"
    # boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
    #           "00-5", "00-6", "00-7", "00-8", "00-9",
    #           "00-10", "00-11", "00-12", "00-13", "00-14",
    #           "00-15", "00-16", "00-17", "00-18", "00-19"]
    boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
              "00-5", "00-6", "00-7", "00-8", "00-9",
              "00-10", "00-11", "00-12", "00-14",
              "00-15", "00-16", "00-17", "00-18", "00-19"]
    halotype = "rockstar"
    scale_variables = ["gamma_l", "gamma_n"]
    N_grids = [550, 1375]
    N_grid_def = 1100
    smooth_type_def = "tophat"
    smooth_types = ["tophat", "gaussian"]
    R_smooths = [3., 7.]
    R_smooth_def = 5.
    scaling_method = "split"

    fig_width = 10.
    fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_width * 0.85),
                             dpi=300, sharex=True)

    gamma_vals = [1.2, 0.8]
    linestyles = ["-", "--"]

    for i, scale_variable in enumerate(scale_variables):
        for j, gamma_val in enumerate(gamma_vals):
            mono_comps = np.empty((len(boxids), 9))
            quad_comps = np.empty((len(boxids), 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_all-halos_{}".format(scaling_method,
                                                     scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = ((red_dv[:9, 1] - ref_dv[:9, 1]) /
                                    ref_dv[:9, 1])
                quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            if i == 1 and j == 0:
                axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C1", marker=".",
                                linestyle=linestyles[j],
                                label="Unsmoothed")
            elif i == 0:
                axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C1", marker=".",
                                linestyle=linestyles[j],
                                label=(r"$\gamma=$" + str(gamma_val)))
            else:
                axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C1", marker=".",
                                linestyle=linestyles[j])
            axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
                            color="C1", marker=".",
                            linestyle=linestyles[j])

            mono_comps = np.empty((len(boxids), 9))
            quad_comps = np.empty((len(boxids), 9))
            for k, boxid in enumerate(boxids):
                path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                             "{}_{}_products/{}_{}_{}_halos/z0.700/"
                             "em_model_test".format(simname, simname, boxid,
                                                    simname, boxid, halotype))
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
                run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                            "{}".format(scaling_method, N_grid_def,
                                        smooth_type_def,
                                        R_smooth_def, scale_variable))
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))

                mono_comps[k, :] = ((red_dv[:9, 1] - ref_dv[:9, 1]) /
                                    ref_dv[:9, 1])
                quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

            if i == 1 and j == 0:
                # axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
                #                 color="C0", marker=".",
                #                 linestyle=linestyles[j],
                #                 label=(r"Tophat, $R_{TH}=5.0$, "
                #                        r"$L_g=1.0 h^{-1}{\rm Mpc}$ "
                #                        "(Fiducial)"))
                axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C0", marker=".",
                                linestyle=linestyles[j],
                                label="Smoothed")

            else:
                axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
                                color="C0", marker=".",
                                linestyle=linestyles[j])
            axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
                            color="C0", marker=".",
                            linestyle=linestyles[j])

        axes[0, i].set_title(r"Scaling $\{}$".format(scale_variable))
        axes[0, i].axhline(y=0., color='k', linestyle='-')
        axes[0, i].legend()
        axes[1, i].axhline(y=0., color='k', linestyle='-')
        # axes[2, i].axhline(y=0., color='k', linestyle='-')
        axes[1, i].set_xscale("log")
        axes[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")
        axes[0, i].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
        axes[1, i].set_ylabel(r"$\xi_2 - \xi_2^b$")

        # axes[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        # axes[2, 0].set_ylim(-0.1, 0.1)

    axes[0, 0].set_ylim(-0.1, 0.1)
    axes[0, 1].set_ylim(-0.25, 0.25)
    axes[1, 0].set_ylim(-0.08, 0.08)
    axes[1, 1].set_ylim(-0.6, 0.6)

    plt.tight_layout()
    output_root = ("{}/pres_comp".format(path_root))
    plt.savefig("{}_cf_check_rel.jpg".format(output_root))


# def paper_smoothing_comp():
#     """
#     Compare effects of scaling smoothed and unsmoothed linear velocities.

#     Loads and plots the relative change in the monopole and quadrupole when the
#     linear velocity is scaled by a factor of 0.8 and 1.2. Compares the change
#     from scaling the unsmoothed linear velocity and the fiducial smoothed
#     linear velocity.
#     """
#     simname = "AbacusCosmos_1100box_planck"
#     boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
#               "00-5", "00-6", "00-7", "00-8", "00-9",
#               "00-10", "00-11", "00-12", "00-13", "00-14",
#               "00-15", "00-16", "00-17", "00-18", "00-19"]
#     halotype = "rockstar"
#     scale_variables = ["gamma_l", "gamma_n"]
#     N_grid_def = 1100
#     smooth_type_def = "tophat"
#     R_smooth_def = 5.
#     scaling_method = "split"
#     gamma_vals = [1.2, 0.8]

#     linestyles = ["-", "--"]
#     fig_width = 10.
#     fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_width * 0.85),
#                              dpi=300, sharex=True)

#     for i, scale_variable in enumerate(scale_variables):
#         for j, gamma_val in enumerate(gamma_vals):
#             mono_comps = np.empty((len(boxids), 9))
#             quad_comps = np.empty((len(boxids), 9))
#             for k, boxid in enumerate(boxids):
#                 path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
#                              "{}_{}_products/{}_{}_{}_halos/z0.700/"
#                              "em_model_test".format(simname, simname, boxid,
#                                                     simname, boxid, halotype))
#                 ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
#                                     ".dat".format(path_root))
#                 run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
#                             "{}".format(scaling_method, N_grid_def,
#                                         smooth_type_def,
#                                         R_smooth_def, scale_variable))
#                 red_dv = np.loadtxt("{}/{}_{}_red_dv"
#                                     ".dat".format(path_root, run_name,
#                                                   gamma_val))

#                 mono_comps[k, :] = ((red_dv[:9, 1] - ref_dv[:9, 1]) /
#                                     ref_dv[:9, 1])
#                 quad_comps[k, :] = (red_dv[9:18, 1] - ref_dv[9:18, 1])

#             if i == 1 and j == 0:
#                 axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
#                                 color=ref_color, marker=".",
#                                 linestyle=linestyles[j], label="Fiducial")
#             elif i == 0:
#                 axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
#                                 color=ref_color, marker=".",
#                                 linestyle=linestyles[j],
#                                 label=(r"$\gamma=$" + str(gamma_val)))
#             else:
#                 axes[0, i].plot(seps, np.mean(mono_comps, axis=0),
#                                 color=ref_color, marker=".",
#                                 linestyle=linestyles[j])
#             axes[1, i].plot(seps, np.mean(quad_comps, axis=0),
#                             color=ref_color, marker=".",
#                             linestyle=linestyles[j])

#         axes[0, i].set_title(r"Scaling $\{}$".format(scale_variable))
#         axes[0, i].axhline(y=0., color='k', linestyle='-')
#         axes[0, i].legend()
#         axes[1, i].axhline(y=0., color='k', linestyle='-')
#         # axes[2, i].axhline(y=0., color='k', linestyle='-')
#         axes[1, i].set_xscale("log")
#         axes[1, i].set_xlabel(r"$r_p[h^{-1}Mpc]$")
#         axes[0, i].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
#         axes[1, i].set_ylabel(r"$\xi_2 - \xi_2^b$")

#         # axes[2, 0].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
#         # axes[2, 0].set_ylim(-0.1, 0.1)

#     axes[0, 0].set_ylim(-0.1, 0.1)
#     axes[0, 1].set_ylim(-0.25, 0.25)
#     axes[1, 0].set_ylim(-0.08, 0.08)
#     axes[1, 1].set_ylim(-0.6, 0.6)

#     plt.tight_layout()
#     output_root = ("{}/paper_comp".format(path_root))
#     plt.savefig("{}_cf_check_rel.jpg".format(output_root))

#     ###################################

#     run_name = "{}_{}".format(scaling_method, scale_variable)

#     for i, gamma_val in enumerate(gamma_vals):
#         red_dv = np.loadtxt("{}/{}_{}_{}_red_dv"
#                             ".dat".format(path_root, scaling_method,
#                                           scale_variable, gamma_val))
#         red_dv_smoothed = np.loadtxt("{}/{}_smoothed_{}_{}_red_dv"
#                                      ".dat".format(path_root, scaling_method,
#                                                    scale_variable, gamma_val))

#         if i == 0:
#             plot_cfs_relative(axes, red_dv_smoothed, ref_dv, scale_variable,
#                               gamma_val, colors[i],
#                               linestyle="--", line_label="Smoothed")
#             plot_cfs_relative(axes, red_dv, ref_dv, scale_variable,
#                               gamma_val, colors[i],
#                               linestyle="-", line_label="Unsmoothed")
#         else:
#             plot_cfs_relative(axes, red_dv_smoothed, ref_dv, scale_variable,
#                               gamma_val, colors[i],
#                               linestyle="--")
#             plot_cfs_relative(axes, red_dv, ref_dv, scale_variable,
#                               gamma_val, colors[i],
#                               linestyle="-")


paper_comparison(thesis=True)

# Change Log
# v0.2.1, 2023-01-30 - Removed quadrupole from thesis plot
# Many versions of different plots
# v0.2.0, 2022-07-29 - Added paper plot (Approx. since deleted by accident)
# v0.1.1, 2022-02-04 - Copied from mock_cfs.py to make dedicated plotting code
# v0.1.0, 2021-05-11 - Code started
