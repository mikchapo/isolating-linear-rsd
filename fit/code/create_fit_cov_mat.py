"""
Code to create covariance matrix for MCMC fits.

v0.2.0, 2022-11-17 - Updated with test clustering
"""

import matplotlib.pyplot as plt
import numpy as np


def cov_to_cor(cov_matrix):
    """Calculate correlation matrix from covariance matrix."""
    cor_matrix = np.copy(cov_matrix)

    for i in range(cor_matrix.shape[0]):
        for j in range(cor_matrix.shape[1]):
            cor_matrix[i, j] = (cor_matrix[i, j] /
                                np.sqrt(cov_matrix[i, i] * cov_matrix[j, j]))

    return cor_matrix


def cor_to_cov(cor_matrix, cov_mat_diags):
    """Calculate covariance matrix from correlation matrix."""
    cov_matrix = np.copy(cor_matrix)

    for i in range(cov_matrix.shape[0]):
        for j in range(cov_matrix.shape[1]):
            cov_matrix[i, j] = (cov_matrix[i, j] *
                                np.sqrt(cov_mat_diags[i] * cov_mat_diags[j]))

    return cov_matrix


def combine_data_emulator(data_dv, output_path, volume_scaling=1.):
    """Add the emulator error to the data covariance matrix."""
    emulator_error_frac_mono = np.loadtxt("/home/mj3chapm/P2/fit/"
                                          "smoothed_emulator/"
                                          "emulator_error_fractional_mono.dat")
    emulator_error_frac_quad = np.loadtxt("/home/mj3chapm/P2/fit/"
                                          "smoothed_emulator/"
                                          "emulator_error_fractional_quad.dat")
    emulator_error_frac_wp = np.loadtxt("/home/mj3chapm/P2/fit/"
                                        "smoothed_emulator/"
                                        "emulator_error_fractional_wp.dat")
    emulator_error_frac = np.concatenate((emulator_error_frac_mono,
                                          emulator_error_frac_quad,
                                          emulator_error_frac_wp))
    eboss_cm_diag = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                               "data_jk_200_cov_mat_pip_ang_jk-corr_diag_sm.dat")
    cov_mat_diag = eboss_cm_diag * volume_scaling
    data_error = np.sqrt(np.diagonal(cov_mat_diag))
    cor_mat_diag = cov_to_cor(cov_mat_diag)

    emulator_error = emulator_error_frac * np.abs(data_dv[:, 1])
    comb_error = np.sqrt(np.power(data_error, 2.) +
                         np.power(emulator_error, 2.))
    comb_cov_mat = cor_to_cov(cor_mat_diag, np.power(comb_error, 2.))
    np.savetxt(output_path, comb_cov_mat)


calc_eboss = True
calc_uchuu = True
calc_test_hod = True

if calc_eboss:
    print(f"eBOSS Cov Mat.")
    data_dv = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                         "dv_eBOSS_LRG_comb_v7_2_pip_ang_ab.dat")
    output_path = ("/home/mj3chapm/P2/fit/data/"
                   "eboss_jk_200_pip_ang_jk-corr_diag_abacus_smooth-vel_"
                   "cov_mat.dat")
    combine_data_emulator(data_dv, output_path)

if calc_uchuu:
    n = 1.
    V_eff = (n / (1. + n))**2. * 8.
    volume_scaling = 1.28 / V_eff
    print(f"Uchuu Cov Mat., V_eff={V_eff}, Scaling={volume_scaling}")
    data_dv = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                         "uchuu_z0.70_n1_red_dv.dat")
    output_path = ("/home/mj3chapm/P2/fit/data/"
                   "uchuu_z0.70_n1_red_vol-scale_jk-corr_abacus_smooth-vel_"
                   "cov_mat.dat")
    combine_data_emulator(data_dv, output_path, volume_scaling=volume_scaling)

if calc_test_hod:
    n = 1.
    volume_factors = [5, 20, 60]
    for i in range(10):
        data_zz = np.loadtxt("/home/mj3chapm/P2/fit/data/test_clustering/"
                             f"CF_cosmo_0_HOD_{i}_mean.dat")
        data_dv = np.empty((27, 2))
        data_dv[:9, :] = data_zz[:, 2:4]
        data_dv[9:18, 0] = data_zz[:, 2]
        data_dv[9:18, 1] = data_zz[:, 4]
        data_dv[18:, :] = data_zz[:, :2]

        for volume_factor in volume_factors:
            V_eff = (n / (1. + n))**2. * 1.1**3. * volume_factor
            volume_scaling = 1.28 / V_eff
            print(f"Test HOD {i} Cov Mat., V. Factor {volume_factor}, "
                  f"V_eff={V_eff}, V. Scaling={volume_scaling}")
            output_path = ("/home/mj3chapm/P2/fit/data/test_clustering/"
                           f"CF_cosmo_0_HOD_{i}_mean_vol-factor-"
                           f"{volume_factor}_jk-corr_abacus_smooth-vel_cov_mat"
                           ".dat")
            combine_data_emulator(data_dv, output_path,
                                  volume_scaling=volume_scaling)

        output_path = ("/home/mj3chapm/P2/fit/data/test_clustering/"
                       f"CF_cosmo_0_HOD_{i}_mean_vol-factor-eBOSS_jk-corr"
                       "_abacus_smooth-vel_cov_mat.dat")
        combine_data_emulator(data_dv, output_path,
                              volume_scaling=1.)

for i in range(10):
    data_zz = np.loadtxt("/home/mj3chapm/P2/fit/data/test_clustering/"
                         f"CF_cosmo_0_HOD_{i}_mean.dat")
    data_dv = np.empty((27, 2))
    data_dv[:9, :] = data_zz[:, 2:4]
    data_dv[9:18, 0] = data_zz[:, 2]
    data_dv[9:18, 1] = data_zz[:, 4]
    data_dv[18:, :] = data_zz[:, :2]
    np.savetxt("/home/mj3chapm/P2/fit/data/test_clustering/"
               f"CF_cosmo_0_HOD_{i}_mean_dv.dat", data_dv)
