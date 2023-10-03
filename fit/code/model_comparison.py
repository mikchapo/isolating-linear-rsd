# Comparison of several best fit models to the data
# v0.2.1, 2023-06-08 - Added option for thesis plot

# Imports
import matplotlib.pyplot as plt
import numpy as np
import sys
import yaml

# Importing the Aemulus emulator takes 30-120 seconds
# from chi_square_funcs import likelihood_chi_square_calc
import fit_funcs as ff
# import GP.GP_MCMC.GP_prediction as aemulus


# Load Params
input_file = sys.argv[1]
with open(input_file) as file:
    input_dict = yaml.load(file, Loader=yaml.FullLoader)

sys.path.insert(0, "/home/mj3chapm/P2/fit/{}".format(input_dict["emulator"]))
# Importing the emulator takes 30-120 seconds
import GP_prediction as emulator


zz_seps = np.loadtxt("/home/mj3chapm/P2/fit/{}/mean_mps_seps"
                     ".dat".format(input_dict["emulator"]))


def emulator_preds(params, z_data, include_ap=True):
    D_M_fid = ff.D_M(ff.Omega_M_fid, z_data)
    D_H_fid = ff.D_H(ff.Omega_M_fid, z_data)
    mono_pred = emulator.mono_pre(params)
    quad_pred = emulator.quad_pre(params) / (zz_seps**2.)
    wp_pred = emulator.wp_pre(params)

    if include_ap:
        alpha_perp = ff.D_M(params[0], z_data) / D_M_fid
        alpha_para = ff.D_H(params[0], z_data) / D_H_fid
        mono_pred, quad_pred, wp_pred = ff.ap_scaling(mono_pred, quad_pred,
                                                      wp_pred, alpha_perp,
                                                      alpha_para)

    return mono_pred, quad_pred, wp_pred


eboss_data = np.loadtxt("/home/mj3chapm/P2/fit/data/dv_eBOSS_LRG_comb_v7_2_pip"
                        "_ang_ab.dat")[:, 1]
eboss_cm = np.loadtxt("/home/mj3chapm/P2/fit/data/eboss_jk_200_pip_ang_diag_"
                      "abacus_smooth-vel_cov_mat.dat")
eboss_stds = np.sqrt(np.diag(eboss_cm))


def model_comparison(model_preds, model_labels, model_zorders, output_path,
                     thesis):
    """Compare emulator models to data."""
    if thesis:
        fig_width = 6.375
        aspect_ratio = 16. / 27.
        height_ratios = [3, 1]

    else:
        fig_width = 10.
        aspect_ratio = 4 / 9.5
        height_ratios = [2, 1]

    fig, axes = plt.subplots(2, 3, sharex=True, dpi=300,
                             figsize=(fig_width, fig_width * aspect_ratio),
                             gridspec_kw={'hspace': 0,
                                          'height_ratios': height_ratios})

    axes[1, 0].axhline(y=0., color='k', linestyle='-')
    axes[1, 1].axhline(y=0., color='k', linestyle='-')
    axes[1, 2].axhline(y=0., color='k', linestyle='-')

    for i in range(len(model_preds)):
        axes[0, 0].plot(ff.seps, ff.seps**2.*model_preds[i][:9],
                        zorder=model_zorders[i])
        axes[1, 0].plot(ff.seps,
                        (model_preds[i][:9] - eboss_data[:9])/eboss_stds[:9],
                        zorder=model_zorders[i])
        axes[0, 1].plot(ff.seps, ff.seps**2.*model_preds[i][9:18],
                        zorder=model_zorders[i], label=model_labels[i])
        axes[1, 1].plot(ff.seps, (model_preds[i][9:18] -
                                  eboss_data[9:18]) / eboss_stds[9:18],
                        zorder=model_zorders[i])
        axes[0, 2].plot(ff.seps, model_preds[i][18:], zorder=model_zorders[i])
        axes[1, 2].plot(ff.seps, (model_preds[i][18:] -
                                  eboss_data[18:]) / eboss_stds[18:],
                        zorder=model_zorders[i])

    axes[0, 0].errorbar(ff.seps, ff.seps**2.*eboss_data[:9],
                        yerr=ff.seps**2.*eboss_stds[:9], fmt="k.",
                        label="eBOSS", zorder=15)
    axes[0, 1].errorbar(ff.seps, ff.seps**2.*eboss_data[9:18],
                        yerr=ff.seps**2.*eboss_stds[9:18], fmt="k.", zorder=15)
    axes[0, 2].errorbar(ff.seps, eboss_data[18:], yerr=eboss_stds[18:],
                        fmt="k.", zorder=15)

    axes[0, 1].legend(handlelength=1)

    axes[0, 2].set_yscale("log")

    axes[1, 0].set_xscale("log")
    axes[1, 0].set_ylim([-3, 3])
    axes[1, 0].axhspan(-1, 1, alpha=0.5, color='grey')
#    axes[1, 0].set_ylabel(r"$\frac{\xi_0^M - \xi_0^D}{\sigma^D(\xi_0)}$")
    axes[1, 0].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[1, 1].set_xscale("log")
    axes[1, 1].set_ylim([-3, 3])
    axes[1, 1].axhspan(-1, 1, alpha=0.5, color='grey')
#    axes[1, 1].set_ylabel(r"$\frac{\xi_2^M - \xi_2^D}{\sigma^D(\xi_2)}$")
    axes[1, 1].set_xlabel(r"$s[h^{-1}Mpc]$")

    axes[1, 2].set_xscale("log")
    axes[1, 2].set_ylim([-3, 3])
    axes[1, 2].axhspan(-1, 1, alpha=0.5, color='grey')
#    axes[1, 2].set_ylabel(r"$\frac{w_p^M - w_p^D}{\sigma^D(w_p)}$")
    axes[1, 2].set_xlabel(r"$r_p[h^{-1}Mpc]$")

    if thesis:
        axes[0, 0].set_title(r"$s^2\xi_0$")
        axes[0, 1].set_title(r"$s^2\xi_2$")
        axes[0, 2].set_title(r"$w_p$")
        axes[1, 0].set_ylabel(r"$(\xi^M - \xi^D)/\sigma^D(\xi)$")

        axes[1, 0].set_yticks([-2, 0, 2])
        axes[1, 1].set_yticks([-2, 0, 2])
        axes[1, 2].set_yticks([-2, 0, 2])

    else:
        axes[0, 0].set_ylabel(r"$s^2\xi_0$")
        axes[0, 1].set_ylabel(r"$s^2\xi_2$")
        axes[0, 2].set_ylabel(r"$w_p$")

        axes[1, 0].set_ylabel(r"$(\xi_0^M - \xi_0^D)/\sigma^D(\xi_0)$")
        axes[1, 1].set_ylabel(r"$(\xi_2^M - \xi_2^D)/\sigma^D(\xi_2)$")
        axes[1, 2].set_ylabel(r"$(w_p^M - w_p^D)/\sigma^D(w_p)$")

    plt.tight_layout()
    plt.savefig(output_path)


model_preds = []
model_labels = []
model_zorders = []
for run in input_dict["runs"]:
    run["params"][3] = run["params"][3] / 100.
    run["params"][1] = run["params"][1] / run["params"][3]**2.
    model_pred = np.concatenate(emulator_preds(run["params"], run["z"],
                                               include_ap=run["include_ap"]))
    # chi_square = likelihood_chi_square_calc(eboss_data, model_pred, eboss_cm,
    #                                         200)
    # print("{} chi square: {}".format(run["run_name"], chi_square))
    model_preds.append(model_pred)
    model_labels.append(run["run_name"])
    model_zorders.append(run["zorder"])

model_comparison(model_preds, model_labels, model_zorders,
                 input_dict["output_path"], input_dict["thesis"])

# Change Log
# v0.2.0, 2022-12-01 - Copied from RSD/fit/aemulus_fmax/
# v0.1.3, 2022-03-02 - Changed y-label of lower panel to one line to make it
#                      more readable
# v0.1.2, 2021-05-31 - Added model zorder
# v0.1.1, 2021-05-03 - Added chi square calculation and print
# v0.1.0, 2021-04-15 - Code started with snippet from v_bc_check.ipynb
