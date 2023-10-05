"""
Code for calculating the jackknife covariance matrix.

v0.1.0, 2022-12-22 - Code started from jackknife_cov_mat.ipynb
"""

# Imports
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


def calc_cov_mat(data_matrices, plot=True,
                 output="../output/plots/data_cov_mat.jpg",
                 title="eBOSS LRG v7 jackknife covariance matrix",
                 rescale=False, ratio=1, textsize=18):
    """Calculate the covariance matrix of a set of data."""
    M = data_matrices.shape[0]
    print(M)
    n = data_matrices.shape[1]

    data_means = np.empty(n)
    for i in range(n):
        data_means[i] = np.mean(data_matrices[:, i, 1])

    cov_matrix = np.zeros((n, n))
    variances = np.zeros(n)
    for i in range(n):
        for j in range(n):
            for k in range(M):
                cov_matrix[i, j] += ((data_matrices[k, i, 1] - data_means[i]) *
                                     (data_matrices[k, j, 1] - data_means[j]))

                if i == j:
                    variances[i] += ((data_matrices[k, i, 1] - data_means[i]) *
                                     (data_matrices[k, i, 1] - data_means[i]))

#             if i == 0:
#                 print("Before Normalization:", cov_matrix[i, j])
#                 cov_matrix[i, j] = cov_matrix[i, j] / (M - 1.)
#                 print(" After Normalization:", cov_matrix[i, j])

            else:
                cov_matrix[i, j] = cov_matrix[i, j] / (M - 1.)

    if rescale:
        cov_matrix = cov_matrix*ratio
    variances = variances / (M - 1.)

    if plot:
        # vmin = np.min(cov_matrix)
        # vmax = np.max(cov_matrix)
        # print(vmin)
        # print(vmax)
        plt.figure(figsize=(12, 12))
        plt.imshow(cov_matrix, origin="lower")
        plt.plot([8.5, 8.5], [-0.5, 26.5], 'k--')
        plt.plot([17.5, 17.5], [-0.5, 26.5], 'k--')
        plt.plot([-0.5, 26.5], [8.5, 8.5], 'k--')
        plt.plot([-0.5, 26.5], [17.5, 17.5], 'k--')
        plt.text(4, -3, r"$\xi_0$", fontsize=textsize)
        plt.text(13, -3, r"$\xi_2$", fontsize=textsize)
        plt.text(22, -3, r"$w_p$", fontsize=textsize)
        plt.text(-3, 4, r"$\xi_0$", fontsize=textsize)
        plt.text(-3, 13, r"$\xi_2$", fontsize=textsize)
        plt.text(-3, 22, r"$w_p$", fontsize=textsize)
        plt.colorbar()
        plt.title(title, fontsize=textsize)
        plt.savefig(output)

    return cov_matrix


def calc_cor_mat(cov_matrix, plot=True,
                 output="../output/plots/data_cor_mat.jpg",
                 title="eBOSS LRG v7 jackknife correlation matrix",
                 textsize=18):
    """Calculate correlation matrix from covariance matrix."""
    cor_matrix = np.copy(cov_matrix)

    for i in range(cor_matrix.shape[0]):
        for j in range(cor_matrix.shape[1]):
            cor_matrix[i, j] = cor_matrix[i, j] / np.sqrt(cov_matrix[i, i] *
                                                          cov_matrix[j, j])

    if plot:
        plt.figure(figsize=(12, 12))
        plt.imshow(cor_matrix, origin="lower", vmin=-1., vmax=1.)
        plt.plot([8.5, 8.5], [-0.5, 26.5], 'k--')
        plt.plot([17.5, 17.5], [-0.5, 26.5], 'k--')
        plt.plot([-0.5, 26.5], [8.5, 8.5], 'k--')
        plt.plot([-0.5, 26.5], [17.5, 17.5], 'k--')
        plt.text(4, -3, r"$\xi_0$", fontsize=textsize)
        plt.text(13, -3, r"$\xi_2$", fontsize=textsize)
        plt.text(22, -3, r"$w_p$", fontsize=textsize)
        plt.text(-3, 4, r"$\xi_0$", fontsize=textsize)
        plt.text(-3, 13, r"$\xi_2$", fontsize=textsize)
        plt.text(-3, 22, r"$w_p$", fontsize=textsize)
        plt.colorbar()
        plt.title(title, fontsize=textsize)
        plt.savefig(output)

#     print(cor_matrix)
    return cor_matrix


def diagonal_smoothing(cov_mat, n=2, lims=[9, 18], save=False,
                       output_root="../output/"):
    """Smooth covariance matrix along diagonal."""
    cor_mat = calc_cor_mat(cov_mat, plot=False)
    sm_cor_mat = np.copy(cor_mat)
    tot_mat = np.zeros(cor_mat.shape)
    num_mat = np.zeros(cor_mat.shape)

    index_lists = [range(0, lims[0]), range(lims[0], lims[1]),
                   range(lims[1], cov_mat.shape[0])]
    for list1 in index_lists:
        for list2 in index_lists:
            for i in list1:
                for j in list2:
                    if i != j:
                        tot = cor_mat[i, j]
                        num = 1
                        for p in range(1, n):
                            if (i - p) >= list1[0] and (j - p) >= list2[0]:
                                tot += cor_mat[i-p, j-p]
                                num += 1

                            if (i + p) <= list1[-1] and (j + p) <= list2[-1]:
                                tot += cor_mat[i+p, j+p]
                                num += 1

                        tot_mat[i, j] = tot
                        num_mat[i, j] = num
                        sm_cor_mat[i, j] = tot / float(num)

    sm_cov_mat = np.copy(sm_cor_mat)

    for i in range(sm_cov_mat.shape[0]):
        for j in range(sm_cov_mat.shape[1]):
            sm_cov_mat[i, j] = sm_cov_mat[i, j] * np.sqrt(cov_mat[i, i] *
                                                          cov_mat[j, j])

    if save:
        np.savetxt(output_root + "cor_mat.dat", cor_mat)
        np.savetxt(output_root + "sm_cor_mat.dat", sm_cor_mat)
        np.savetxt(output_root + "tot_mat.dat", tot_mat)
        np.savetxt(output_root + "num_mat.dat", num_mat)
        np.savetxt(output_root + "sm.dat", sm_cov_mat)

    return sm_cov_mat


def combine_mats(mat1, mat2):
    """Combine two matrices for plotting."""
    comb_mat = np.copy(mat1)
    for i in range(mat1.shape[0]):
        for j in range(i+1, mat1.shape[0]):
            comb_mat[i, j] = mat2[i, j]

    return comb_mat


def plot_mat(ax, mat, title="", textsize=12, vmin=None, vmax=None, xlabel="",
             ylabel=""):
    """Plot matrix."""
    if vmin is not None:
        ais = ax.imshow(mat, origin="lower", vmin=vmin, vmax=vmax)

    else:
        ais = ax.imshow(mat, origin="lower")

    ax.plot([8.5, 8.5], [-0.5, 26.5], 'k--', linewidth=1)
    ax.plot([17.5, 17.5], [-0.5, 26.5], 'k--', linewidth=1)
    ax.plot([-0.5, 26.5], [8.5, 8.5], 'k--', linewidth=1)
    ax.plot([-0.5, 26.5], [17.5, 17.5], 'k--', linewidth=1)
    ax.text(2, -4, r"$\xi_0$", fontsize=textsize)
    ax.text(12, -4, r"$\xi_2$", fontsize=textsize)
    ax.text(22, -4, r"$w_p$", fontsize=textsize)
    ax.text(-5, 2, r"$\xi_0$", fontsize=textsize)
    ax.text(-5, 12, r"$\xi_2$", fontsize=textsize)
    ax.text(-5, 22, r"$w_p$", fontsize=textsize)

    ax.set_title(title, fontsize=textsize)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=textsize)

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=textsize)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    plt.colorbar(ais, cax=cax)


# Data PIP+ANG JK 200
jk_corr_funcs = np.empty((200, 27, 2))
for i in range(200):
    jk_corr_funcs[i] = np.loadtxt("../output/data_products/"
                                  "eBOSS_pip+ang_jk-corr_200/removed_region_{}"
                                  ".dat".format(i))

R_eff_ngc_assigned = np.loadtxt("../output/data_products/jk_pc/200/"
                                "RR_norm_eBOSS_LRG_NGC_v7_2_pip.dat",
                                skiprows=1)[0]
R_eff_sgc_assigned = np.loadtxt("../output/data_products/jk_pc/200/"
                                "RR_norm_eBOSS_LRG_SGC_v7_2_pip.dat",
                                skiprows=1)[0]
R_eff_ngc_full = np.loadtxt("../output/data_products/jk_pc/201/"
                            "RR_norm_eBOSS_LRG_NGC_v7_2_pip.dat",
                            skiprows=1)[0]
R_eff_sgc_full = np.loadtxt("../output/data_products/jk_pc/201/"
                            "RR_norm_eBOSS_LRG_SGC_v7_2_pip.dat",
                            skiprows=1)[0]
R_assigned = R_eff_ngc_assigned + R_eff_sgc_assigned
R_full = R_eff_ngc_full + R_eff_sgc_full
rand_ratio = R_assigned / R_full * 199. * 199. / 200.
# rand_ratio = R_exclusive / R_full

jk_data_cov_mat = calc_cov_mat(jk_corr_funcs, rescale=True, ratio=rand_ratio,
                               title="eBOSS LRG v7_2 PIP+ANG JK 200 Covariance"
                                     " Matrix")
np.savetxt("../output/data_products/data_jk_200_cov_mat_pip_ang_jk-corr.dat",
           jk_data_cov_mat)

diag_cov_mat = diagonal_smoothing(jk_data_cov_mat, save=True,
                                  output_root="../output/data_products/"
                                              "data_jk_200_cov_mat_pip_ang_"
                                              "jk-corr_diag_")

data_cov_mat_uncorr = np.loadtxt("../data/data_jk_200_cov_mat_pip_ang"
                                 ".dat")
diag_cov_mat_uncorr = np.loadtxt("../data/data_jk_200_cov_mat_pip_ang_"
                                 "diag_sm.dat")

data_cor_mat = calc_cor_mat(jk_data_cov_mat, plot=False)
diag_cor_mat = calc_cor_mat(diag_cov_mat, plot=False)
data_cor_mat_uncorr = calc_cor_mat(data_cov_mat_uncorr, plot=False)
diag_cor_mat_uncorr = calc_cor_mat(diag_cov_mat_uncorr, plot=False)

ylabel = "Original"
xlabel = "Alpha Corrected"

data_cov_comb = combine_mats(jk_data_cov_mat, data_cov_mat_uncorr)
data_cor_comb = combine_mats(data_cor_mat, data_cor_mat_uncorr)
diag_cor_comb = combine_mats(diag_cor_mat, diag_cor_mat_uncorr)

plt.figure(figsize=(5, 4.7), dpi=300)
ax = plt.gca()
plot_mat(ax, data_cov_comb, textsize=12, vmin=None, vmax=None, xlabel=xlabel,
         ylabel=ylabel, title="Correlation matrix comparison")
plt.savefig("../output/plots/data_cov_correction_comp.jpg")

plt.figure(figsize=(5, 4.7), dpi=300)
ax = plt.gca()
plot_mat(ax, data_cor_comb, textsize=12, vmin=-1, vmax=1., xlabel=xlabel,
         ylabel=ylabel, title="Correlation matrix comparison")
plt.savefig("../output/plots/data_cor_correction_comp.jpg")

plt.figure(figsize=(5, 4.7), dpi=300)
ax = plt.gca()
plot_mat(ax, diag_cor_comb, textsize=12, vmin=-1, vmax=1., xlabel=xlabel,
         ylabel=ylabel,
         title="Diagonally smoothed correlation matrix comparison")
plt.savefig("../output/plots/diag_cor_correction_comp.jpg")

data_error_corr = np.sqrt(np.diagonal(jk_data_cov_mat))
data_error_uncorr = np.sqrt(np.diagonal(data_cov_mat_uncorr))

bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))

fig, axes = plt.subplots(1, 3, figsize=(9.5, 3))
axes[0].plot(seps, ((data_error_uncorr[:9] - data_error_corr[:9]) /
                    data_error_corr[:9]))
axes[1].plot(seps, ((data_error_uncorr[9:18] - data_error_corr[9:18]) /
                    data_error_corr[9:18]))
axes[2].plot(seps, ((data_error_uncorr[18:] - data_error_corr[18:]) /
                    data_error_corr[18:]))
axes[0].axhline(0)
axes[1].axhline(0)
axes[2].axhline(0)
axes[0].set_xscale("log")
axes[1].set_xscale("log")
axes[2].set_xscale("log")
plt.savefig("../output/plots/error_comparison.jpg")

# Change Log
