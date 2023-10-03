# Plot the eBOSS LRG correlation function for showing in the defence
# v0.1.0, 2023-08-31 - Copied from model_comparison.py

# Imports
import matplotlib.pyplot as plt
import numpy as np

import fit_funcs as ff


eboss_data = np.loadtxt("/home/mj3chapm/P2/fit/data/dv_eBOSS_LRG_comb_v7_2_pip"
                        "_ang_ab.dat")[:, 1]
eboss_cm = np.loadtxt("/home/mj3chapm/P2/fit/data/eboss_jk_200_pip_ang_diag_"
                      "abacus_smooth-vel_cov_mat.dat")
eboss_stds = np.sqrt(np.diag(eboss_cm))

# fig_width = 3
# aspect_ratio = 4 / 3

fig_width = 6.
aspect_ratio = 2. / 6.

fig, axes = plt.subplots(1, 3, dpi=300, sharex=True,
                         figsize=(fig_width, fig_width * aspect_ratio))

axes[0].errorbar(ff.seps, ff.seps**2.*eboss_data[:9],
                 yerr=ff.seps**2.*eboss_stds[:9], marker=".",
                 label="eBOSS", linestyle='')
axes[1].errorbar(ff.seps, ff.seps**2.*eboss_data[9:18],
                 yerr=ff.seps**2.*eboss_stds[9:18], marker=".", linestyle='')
axes[2].errorbar(ff.seps, eboss_data[18:], yerr=eboss_stds[18:],
                 marker=".", linestyle='')

# axes[0].errorbar(ff.seps, eboss_data[:9],
#                  yerr=ff.seps**2.*eboss_stds[:9], fmt="b.",
#                  label="eBOSS")
# axes[1].errorbar(ff.seps, eboss_data[9:18],
#                  yerr=ff.seps**2.*eboss_stds[9:18], fmt="b.")
# axes[2].errorbar(ff.seps, eboss_data[18:], yerr=eboss_stds[18:],
#                  fmt="b.")

# axes[0].set_xlabel(r"$s\ [h^{-1}\,{\rm Mpc}]$")
# axes[1].set_xlabel(r"$s\ [h^{-1}\,{\rm Mpc}]$")
# axes[2].set_xlabel(r"$r_{\perp}\ [h^{-1}\,{\rm Mpc}]$")

axes[2].set_xlabel(r"$r\ [h^{-1}\,{\rm Mpc}]$")

# axes[0].set_xscale("log")
# axes[1].set_xscale("log")
axes[2].set_xscale("log")

axes[2].set_yscale("log")

# axes[0].set_title(r"$s^2\xi_0$")
# axes[1].set_title(r"$s^2\xi_2$")
# axes[2].set_title(r"$w_p$")

axes[0].set_ylabel(r"$s^2\xi_0$")
axes[1].set_ylabel(r"$s^2\xi_2$")
axes[2].set_ylabel(r"$w_p$")

# axes[0].set_title(r"$\xi_0$")
# axes[1].set_title(r"$\xi_2$")
# axes[2].set_title(r"$w_p$")

plt.tight_layout()
plt.savefig("../output/plots/defence_cf_horizontal.png")

# Change Log
