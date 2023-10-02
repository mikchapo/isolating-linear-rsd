# Plot pairwise velocities for the Abacus halos in mass bins
# v0.1.0, 2021-10-29 - Code started, copied from plot_vp_mass_comp_abacus.py

# Imports
import matplotlib.pyplot as plt
import numpy as np

from vp_plot_funcs import H_z, vp_lin_pred, plot_static, plot_lin_pred


fill_types = ["v-tot", "v-lin", "v-nl"]
cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_planck_products/"
                          "AbacusCosmos_1100box_planck_rockstar_halos/"
                          "info/cosmo_params.dat")

ic_type = "ic_ngc_nplt_z49.0"
z = 0.7
N_bins = 80
bin_edges = np.logspace(np.log10(0.01), np.log10(100.), N_bins+1)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
bin_centres = np.power(10., (np.log10(bin_edges[:-1]) + log_dsep / 2.))
mass_bin_lims = np.arange(12., 14.5, 0.5)

linestyles = ["-", "--", ":"]


for i in range(mass_bin_lims.size-1):
    fig, axes = plt.subplots(2, 1, figsize=(8., 8.), sharex=True,
                             gridspec_kw={'hspace': 0,
                             'height_ratios': [3, 1]}, dpi=300, num=2*i)
    # plt.figure(num=2*i, figsize=(8., 6.))
    axes[0].axhline(y=0., linestyle="-", color="k")
    # ymax = 0.
    # ymin = 0.
    ymax = 100.
    ymin = -500.

    # H_seps = np.logspace(np.log10(0.01), np.log10(5.), 200)
    # H = H_z(z, cosmo_params[0], (cosmo_params[1] + cosmo_params[2]) /
    #         (cosmo_params[0] / 100.)**2.)
    # axes[0].plot(H_seps, -H*H_seps / (cosmo_params[0] / 100.), color="k",
    #              label="Static solution")

    for j in range(3):
        mvps = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_planck_products/"
                          "AbacusCosmos_1100box_planck_rockstar_halos/z0.700/"
                          "halo_{}_{:.2f}-{:.2f}_vp"
                          ".dat".format(fill_types[j],
                                        mass_bin_lims[i], mass_bin_lims[i+1]))
        mvps[:, 1] = np.nan_to_num(mvps[:, 1])
        if fill_types[j] == "v-tot":
            mvps_tot = mvps[:, 1]

        if fill_types[j] == "v-lin":
            mvps_lin = mvps[:, 1]

        axes[0].errorbar(bin_centres, mvps[:, 1],
                         yerr=np.sqrt(mvps[:, 2] / mvps[:, 3]),
                         color="C{}".format(j), label=fill_types[j],
                         marker=".", linestyle=linestyles[j])

    vp_lin_seps = np.logspace(np.log10(1.), np.log10(100.), 40)
    vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1],
                         cosmo_params[2], cosmo_params[3], cosmo_params[4])
    axes[0].plot(vp_lin_seps, vp_lin, color="k", linestyle="--",
                 label="Linear theory prediction")

    axes[1].axhline(y=0., linestyle="-", color="k")
    axes[1].plot(bin_centres, np.nan_to_num((mvps_lin - mvps_tot) / mvps_tot),
                 color="C1", marker=".", linestyle="-")
    axes[1].set_ylabel(r"$(v_p^{lin} - v_p^{tot}) / v_p^{tot}$")
    axes[1].set_ylim(-0.25, 0.25)

    axes[0].set_xlabel(r"s [$h^{-1}$ Mpc]")
    axes[0].set_ylabel(r"$v_p$ [km/s]")
    axes[0].set_xlim(0.1, 100.)
    # axes[0].set_ylim(ymin-50., ymax+50.)
    axes[0].set_xscale("log")
    axes[0].legend()
    axes[0].set_title("{:.2f}-{:.2f} Halo Mass".format(mass_bin_lims[i],
                                                       mass_bin_lims[i+1]))
    plt.savefig("../output/plots/"
                "abacus_rockstar_{:.2f}-{:.2f}"
                ".jpg".format(mass_bin_lims[i], mass_bin_lims[i+1]))


# Change Log
