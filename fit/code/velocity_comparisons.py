"""
Compare effects of non-linear velocity parameters.

v0.1.1, 2023-07-15 - Added option for thesis
"""

# Imports
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import numpy as np
import sys

import fit_funcs as ff

emulator_name = "smoothed_emulator"
sys.path.insert(0, "/home/mj3chapm/P2/fit/{}".format(emulator_name))
# Importing the emulator takes 30-120 seconds
import GP_prediction as emulator

zz_seps = np.loadtxt("/home/mj3chapm/P2/fit/{}/mean_mps_seps"
                     ".dat".format(emulator_name))


# [omegam, omegab, sigma8_emulator, h, ns, w, logM_sat, alpha, logM_cut,
#  sigma_logM, v_bc, v_bs, c_vir, gamma_n, fmax, gamma_l]
default_params = [0.3089, 0.0486, 0.8159, 0.6774, 0.9667, -1.,
                  14.460281, 1.5636618, 13.206889, 0.33605492,
                  0., 1., 1., 1., 0.25026736, 1.]

mono_pred = emulator.mono_pre(default_params)
quad_pred = emulator.quad_pre(default_params) / (zz_seps**2.)
wp_pred = emulator.wp_pre(default_params)
base_pred = np.concatenate((mono_pred, quad_pred, wp_pred))

linestyles = ["-", "--"]


def start_plot_rel(wp=False, thesis=False):
    """Start plot of relative only."""
    if wp:
        if thesis:
            fig_width = 6.375
            aspect_ratio = 16. / 27.
            fig, axes = plt.subplots(1, 3, dpi=300,
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

        else:
            fig_width = 10.
            aspect_ratio = 4.5 / 9.5
            fig, axes = plt.subplots(1, 3, dpi=300,
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

    else:
        if thesis:
            fig_width = 6.375
            aspect_ratio = 16. / 27.
            fig, axes = plt.subplots(1, 2, dpi=300,
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

        else:
            fig_width = 7.
            aspect_ratio = 4.5 * 10. / 9.5 / 7.
            fig, axes = plt.subplots(1, 2, dpi=300,
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))
    return fig, axes


def start_plot_comb(wp=False, thesis=False):
    """Start plot of direct and relative."""
    if wp:
        if thesis:
            fig_width = 6.375
            aspect_ratio = 16. * 6 / 27. / 4.5
            fig, axes = plt.subplots(2, 3, dpi=300, sharex=True,
                                     gridspec_kw={'hspace': 0},
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

        else:
            fig_width = 10.
            aspect_ratio = 6. / 9.5
            fig, axes = plt.subplots(2, 3, dpi=300, sharex=True,
                                     gridspec_kw={'hspace': 0},
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

    else:
        if thesis:
            fig_width = 6.375
            aspect_ratio = 16. * 6 / 27. / 4.5
            fig, axes = plt.subplots(2, 2, dpi=300, sharex=True,
                                     gridspec_kw={'hspace': 0},
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

        else:
            fig_width = 7.
            aspect_ratio = 6. * 10. / 9.5 / 7.
            fig, axes = plt.subplots(2, 2, dpi=300, sharex=True,
                                     gridspec_kw={'hspace': 0},
                                     figsize=(fig_width,
                                              fig_width * aspect_ratio))

    return fig, axes


def plot_dir(ax, data, color, linestyle, label=None, wp=False):
    """Plot direct."""
    if wp:
        ax.plot(ff.seps, data, color=color, linestyle=linestyle,
                label=label, marker=".")

    else:
        ax.plot(ff.seps, ff.seps**2.*data, color=color, linestyle=linestyle,
                label=label, marker=".")


def plot_rel(ax, data, base, color, linestyle, label=None):
    """Plot relative."""
    ax.plot(ff.seps, (data - base)/base, color=color, linestyle=linestyle,
            label=label, marker=".")


def finish_plot_rel(axes, wp=False):
    """Finish plot of relative only."""
    axes[0].axhline(y=0., color='k', linestyle='-')
    axes[0].set_xscale("log")
    # axes[0].set_ylim([-10, 10])
    axes[0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
    axes[0].set_xlabel(r"$s[h^{-1}Mpc]$")
    axes[0].legend()

    axes[1].axhline(y=0., color='k', linestyle='-')
    axes[1].set_xscale("log")
    # axes[1].set_ylim([-10, 10])
    # axes[1].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes[1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")
    axes[1].set_xlabel(r"$s[h^{-1}Mpc]$")

    if wp:
        axes[1].axhline(y=0., color='k', linestyle='-')
        axes[1].set_xscale("log")
        axes[1].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes[1].set_xlabel(r"$r_{\perp}[h^{-1}Mpc]$")

    plt.tight_layout()


def finish_plot_comb(axes, wp=False, thesis=False):
    """Start plot of direct and relative."""
    axes[1, 0].axhline(y=0., color='k', linestyle='-')
    axes[1, 0].set_xscale("log")
    # axes[0].set_ylim([-10, 10])
    axes[1, 0].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[1, 1].axhline(y=0., color='k', linestyle='-')
    axes[1, 1].set_xscale("log")
    # axes[1, 1].set_ylim([-10, 10])
    # axes[1].set_ylabel(r"$\xi_2 - \xi_2^b$")
    axes[1, 1].set_xlabel(r"$s[h^{-1}Mpc]$")

    if wp:
        axes[0, 2].set_yscale("log")

        axes[1, 2].axhline(y=0., color='k', linestyle='-')
        axes[1, 2].set_xscale("log")
        axes[1, 2].set_xlabel(r"$r_{\perp}[h^{-1}Mpc]$")

    if thesis:
        axes[0, 1].legend(handlelength=1)

        axes[0, 0].set_title(r"$s^2\xi_0$")
        axes[0, 1].set_title(r"$s^2\xi_2$")
        axes[1, 0].set_ylabel(r"$(\xi - \xi^b)/\xi^b$")

        if wp:
            axes[0, 2].set_title(r"$w_p$")

    else:
        axes[0, 0].legend()

        axes[0, 0].set_ylabel(r"$s^2\xi_0$")
        axes[0, 1].set_ylabel(r"$s^2\xi_2$")

        axes[1, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
        axes[1, 1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")

        if wp:
            axes[0, 2].set_ylabel(r"$w_p$")
            axes[1, 2].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")

    plt.tight_layout()


def plot_variation(variations, name, rel_only=True, wp=False, thesis=False):
    """Plot a set of variations."""
    if rel_only:
        fig, axes = start_plot_rel(wp=wp)

    else:
        fig, axes = start_plot_comb(wp=wp, thesis=thesis)
        plot_dir(axes[0, 0], base_pred[:9], "k", "-", label="Base")
        plot_dir(axes[0, 1], base_pred[9:18], "k", "-")
        if wp:
            plot_dir(axes[0, 2], base_pred[18:], "k", "-", wp=True)

    for i, variation in enumerate(variations):
        var_params = default_params.copy()
        for j, val in enumerate(variation["values"]):
            var_params[variation["index"]] = val
            mono_pred = emulator.mono_pre(var_params)
            quad_pred = emulator.quad_pre(var_params) / (zz_seps**2.)
            wp_pred = emulator.wp_pre(var_params)

            if variation["label"] == r"$\Omega_m$":
                base_cosmo = FlatLambdaCDM(default_params[3]*100.,
                                           default_params[0], Tcmb0=2.7255)
                var_cosmo = FlatLambdaCDM(default_params[3]*100., val,
                                          Tcmb0=2.7255)

                D_M_base = base_cosmo.angular_diameter_distance(0.7)
                D_M_var = var_cosmo.angular_diameter_distance(0.7)

                alpha_perp = D_M_var / D_M_base
                alpha_para = 1.

                (mono_pred, quad_pred,
                 wp_pred) = ff.ap_scaling(mono_pred, quad_pred, wp_pred,
                                          alpha_perp, alpha_para)

            var_pred = np.concatenate((mono_pred, quad_pred, wp_pred))
            label = "{}={}".format(variation["label"], variation["values"][j])

            if rel_only:
                plot_rel(axes[0], var_pred[:9], base_pred[:9], "C{}".format(i),
                         linestyles[j], label=label)
                plot_rel(axes[1], var_pred[9:18], base_pred[9:18],
                         "C{}".format(i), linestyles[j])
                if wp:
                    plot_rel(axes[2], var_pred[18:], base_pred[18:],
                             "C{}".format(i), linestyles[j])

            else:
                if thesis:
                    plot_dir(axes[0, 0], var_pred[:9], "C{}".format(i),
                             linestyles[j])
                    plot_dir(axes[0, 1], var_pred[9:18], "C{}".format(i),
                             linestyles[j], label=label)

                else:
                    plot_dir(axes[0, 0], var_pred[:9], "C{}".format(i),
                             linestyles[j], label=label)
                    plot_dir(axes[0, 1], var_pred[9:18], "C{}".format(i),
                             linestyles[j])

                plot_rel(axes[1, 0], var_pred[:9], base_pred[:9],
                         "C{}".format(i), linestyles[j])
                plot_rel(axes[1, 1], var_pred[9:18], base_pred[9:18],
                         "C{}".format(i), linestyles[j])
                if wp:
                    plot_dir(axes[0, 2], var_pred[18:], "C{}".format(i),
                             linestyles[j], wp=True)
                    plot_rel(axes[1, 2], var_pred[18:], base_pred[18:],
                             "C{}".format(i), linestyles[j])

    if rel_only:
        finish_plot_rel(axes, wp=wp)

    else:
        finish_plot_comb(axes, wp=wp, thesis=thesis)

    if rel_only:
        if wp:
            plt.savefig("../output/plots/clustering_effects/"
                        "{}_vel_comp_wp_rel.png".format(name))

        else:
            plt.savefig("../output/plots/clustering_effects/"
                        "{}_vel_comp_rel.png".format(name))

    else:
        if wp:
            if thesis:
                plt.savefig("../output/plots/clustering_effects/"
                            "{}_vel_comp_wp_comb_thesis.png".format(name))

            else:
                plt.savefig("../output/plots/clustering_effects/"
                            "{}_vel_comp_wp_comb.png".format(name))

        else:
            if thesis:
                plt.savefig("../output/plots/clustering_effects/"
                            "{}_vel_comp_comb_thesis.png".format(name))

            else:
                plt.savefig("../output/plots/clustering_effects/"
                            "{}_vel_comp_comb.png".format(name))


################# NON-LINEAR VELOCITY VARIATION ###################

nl_variations = [{"label": r"$\gamma_n$", "values": [1.2, 0.8], "index": 13},
                 {"label": r"$v_{\rm bc}$", "index": 10, "values": [0.4]},
                 {"label": r"$v_{\rm bs}$", "index": 11, "values": [1.4, 0.6]}]

plot_variation(nl_variations, "nl", rel_only=False, thesis=True)


################# LINEAR VELOCITY VARIATION ###################

lin_variations = [{"label": r"$\gamma_l$", "index": 15, "values": [1.1, 0.9]},
                  {"label": r"$\Omega_m$", "index": 0,
                   "values": [0.32, 0.30]},
                  {"label": r"$\sigma_8$", "index": 2,
                   "values": [0.84, 0.79]}]

plot_variation(lin_variations, "lin", rel_only=False, wp=True, thesis=True)

# Change Log
# v0.1.0, 2023-06-19 - Code started with snippets from
#                      gammas_clustering_effects.py
