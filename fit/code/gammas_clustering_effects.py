# Test effect of fsigma8 on small scales
# v0.3.0, 2023-06-16 - Updated to compare gamma_n, v_bc, v_bs


# Imports
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys

import fit_funcs as ff


emulator_name = "smoothed_emulator"
sys.path.insert(0, "/home/mj3chapm/P2/fit/{}".format(emulator_name))
# Importing the emulator takes 30-120 seconds
import GP_prediction as emulator

zz_seps = np.loadtxt("/home/mj3chapm/P2/fit/{}/mean_mps_seps"
                     ".dat".format(emulator_name))


colors_prime = sns.color_palette("viridis")

# params = [0.32359125, 0.022694162 / 0.67135187**2., 0.88166519, 0.67135187, 0.98832471, 3.046, -1., 14.413626, 0.39430394, 12.42274,
#           1.0892899, 0.040618409, 0.97836545, 0.69443739, 0.91985128, 0.98916452]
# variations = [{"label": "v_bc", "values": [0.2, 0.4, 0.6, 0.8], "index": 11, "base_label": "v_bc=0.04 (Base)"},
#               {"label": "sigma8", "values": [0.80, 0.84, 0.92, 0.96], "index": 2, "base_label": "sigma8=0.88 (Base)"},
#               {"label": "gamma_f", "values": [0.76, 0.84, 1., 1.08], "index": 14, "base_label": "gamma_f=0.0.92 (Base)"}]


sample_params = [0.3089, 0.0486, 0.8159, 0.6774, 0.9667, -1., 15.05624872366684, 0.9368163, 12.828254473234338, 1.353408, 0., 1., 1., 1., 1., 1.]
# hod_schemes = [{"label": "test389", "hod_params": [15.05624872366684, 0.9368163, 12.828254473234338, 1.353408]},
#                {"label": "test579", "hod_params": [14.017709727072987, 0.650746, 12.720192391248867, 1.09601]},
#                {"label": "test666", "hod_params": [14.1014528634932, 0.3249109, 12.084130674162253, 0.3509227]},
#                {"label": "test825", "hod_params": [15.732143888202852, 0.7981138, 10.427414369926954, 1.344095]},
#                {"label": "test986", "hod_params": [15.948381126130892, 0.3821675, 10.666587518263391, 0.160913]},
#                {"label": "eboss_base", "hod_params": [1.4499227E+01, 1.3747280, 1.3397143E+01, 4.0994775E-01]}]
hod_schemes = [{"label": "eboss_base", "hod_params": [1.4499227E+01, 1.3747280, 1.3397143E+01, 4.0994775E-01]}]

for hod_scheme in hod_schemes:
    params = sample_params.copy()
    params[6:10] = hod_scheme["hod_params"]
    mono_pred = emulator.mono_pre(params)
    quad_pred = emulator.quad_pre(params) / (zz_seps**2.)
    wp_pred = emulator.wp_pre(params)
    base_pred = np.concatenate((mono_pred, quad_pred, wp_pred))

#     variations = [{"label": "v_bc", "values": [0.2, 0.4, 0.6, 0.8], "index": 11, "base_label": "v_bc=0. (Base)"},
#                   {"label": "sigma8", "values": [0.7359, 0.7759, 0.8559, 0.8959], "index": 2, "base_label": "sigma8=0.8159 (Base)"},
#                   {"label": "gamma_f", "values": [0.84, 0.92, 1.08, 1.16], "index": 14, "base_label": "gamma_f=1 (Base)"}]
    variations = [{"label": r"$\gamma_n$", "plain_label": "gamma_n", "values": [0.8, 0.9, 1.1, 1.2], "index": 13, "base_label": r"$\gamma_n=1$ (Base)"},
                  {"label": r"$v_{\rm bc}$", "plain_label": "v_bc", "values": [0.2, 0.4, 0.6, 0.8], "index": 10, "base_label": r"$v_{\rm bc}=0$ (Base)"},
                  {"label": r"$v_{\rm bs}$", "plain_label": "v_bs", "values": [0.6, 0.8, 1.2, 1.4], "index": 11, "base_label": r"$v_{\rm bs}=0$ (Base)"},
                  {"label": r"$\gamma_l$", "plain_label": "gamma_l", "values": [0.8, 0.9, 1.1, 1.2], "index": 15, "base_label": r"$\gamma_l=1$ (Base)"}]

    for variation in variations:
        colors = colors_prime.copy()
        fig_width = 10.
        fig, axes = plt.subplots(2, 3, figsize=(fig_width, fig_width / 9.5 * 4.), sharex=True, gridspec_kw={'hspace': 0, 'height_ratios': [2, 1]}, dpi=300)
        var_params = params.copy()
        for i in range(len(variation["values"])):
            var_params[variation["index"]] = variation["values"][i]
            mono_pred = emulator.mono_pre(var_params)
            quad_pred = emulator.quad_pre(var_params) / (zz_seps**2.)
            wp_pred = emulator.wp_pre(var_params)
            var_pred = np.concatenate((mono_pred, quad_pred, wp_pred))

            print("axes size:", axes.size)
            print("var_pred size:", var_pred.size)
            print("colors len:", len(colors))
            axes[0, 0].plot(ff.seps, ff.seps**2.*var_pred[:9], color=colors[i], marker=".")
            # axes[1, 0].plot(ff.seps, (ff.seps*var_pred[:9] - ff.seps*base_pred[:9]))
            axes[1, 0].plot(ff.seps, (var_pred[:9] - base_pred[:9])/base_pred[:9], color=colors[i], marker=".")
            axes[0, 1].plot(ff.seps, ff.seps**2.*var_pred[9:18], label="{}={}".format(variation["label"], variation["values"][i]), color=colors[i], marker=".")
            # axes[1, 1].plot(ff.seps, (var_pred[9:18] - base_pred[9:18]))
            axes[1, 1].plot(ff.seps, (var_pred[9:18] - base_pred[9:18])/base_pred[9:18], color=colors[i], marker=".")
            axes[0, 2].plot(ff.seps, var_pred[18:], color=colors[i], marker=".")
            # axes[1, 2].plot(ff.seps, (var_pred[18:] - base_pred[18:]))
            axes[1, 2].plot(ff.seps, (var_pred[18:] - base_pred[18:])/base_pred[18:], color=colors[i], marker=".")

            if len(variation["values"]) / (i+1) == 2:
                color = colors.pop(i+1)
                axes[0, 0].plot(ff.seps, ff.seps**2.*base_pred[:9], color=color, marker=".")
                axes[0, 1].plot(ff.seps, ff.seps**2.*base_pred[9:18], label=variation["base_label"], color=color, marker=".")
                axes[0, 2].plot(ff.seps, base_pred[18:], color=color, marker=".")

        axes[0, 0].set_ylabel(r"$s^2\xi_0$")

        axes[1, 0].axhline(y=0., color='k', linestyle='-')
        axes[1, 0].set_xscale("log")
        # axes[1, 0].set_ylim([-10, 10])
        # axes[1, 0].set_ylabel(r"$s\xi_0 - s\xi_0^b$")
        axes[1, 0].set_ylabel(r"$(\xi_0 - \xi_0^b)/\xi_0^b$")
        axes[1, 0].set_xlabel(r"$s[h^{-1}Mpc]$")

        axes[0, 1].set_ylabel(r"$s^2\xi_2$")

        axes[1, 1].axhline(y=0., color='k', linestyle='-')
        axes[1, 1].set_xscale("log")
        # axes[1, 1].set_ylim([-10, 10])
        # axes[1, 1].set_ylabel(r"$\xi_2 - \xi_2^b$")
        axes[1, 1].set_ylabel(r"$(\xi_2 - \xi_2^b)/\xi_2^b$")
        axes[1, 1].set_xlabel(r"$s[h^{-1}Mpc]$")

        axes[0, 2].set_yscale("log")
        axes[0, 2].set_ylabel(r"$w_p$")
        axes[0, 1].legend()

        axes[1, 2].axhline(y=0., color='k', linestyle='-')
        axes[1, 2].set_xscale("log")
        # axes[1, 2].set_ylim([-10, 10])
        # axes[1, 2].set_ylabel(r"$w_p - w_p^b$")
        axes[1, 2].set_ylabel(r"$(w_p - w_p^b)/w_p^b$")
        axes[1, 2].set_xlabel(r"$r_p[h^{-1}Mpc]$")

        plt.tight_layout()
        plt.savefig("../output/plots/clustering_effects/cf_check_{}_{}.png".format(hod_scheme["label"], variation["plain_label"]))


# Change Log
# v0.2.0, 2022-01-21 - Copied from
#                      RSD/fit/aemulus_fmax/fsig8_clustering_effects.py
