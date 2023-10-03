# Measure the correlation functions of mocks
# v0.3.0, 2022-02-04 - Copied from RSD to P2 for testing new emulator models

from Corrfunc.theory import DDrppi, DDsmu, wp
from scipy.integrate import quad
import datetime as dt
import numpy as np
import sys


def comoving_dist_integral(z):
    return 1. / np.sqrt(Omega_M * (1. + z) * (1. + z) * (1. + z) + Omega_L)


def gen_red_cat(real_cat_path, output_root, zs_real, boxsize, save_output,
                Omega_M=0.3089, zerr=False, zerr_std=91.8,
                scaling_method="combined", gamma_l=1., gamma_n=1., tc=0.5,
                tb=0.2):
    start_time = dt.datetime.now()
    print("Generating Redshift Catalogue, Start Time:", start_time)

    print("Output Root:", output_root)

    # Hardcoded to flat, matter only
    Omega_L = 1. - Omega_M
    c = 299792.458

    if isinstance(real_cat_path, str):
        real_cat = np.loadtxt(real_cat_path)
    elif isinstance(real_cat_path, np.ndarray):
        real_cat = real_cat_path
    else:
        raise ValueError("Unrecognized type of real_cat_path."
                         "Only str or np.ndarray are allowed.")

    print("Catalogue Loaded:", dt.datetime.now() - start_time)

    # Fudge factor for the moment because this emulator set seems to have ended
    # up with range 550 to 1650
    xmin = np.floor(np.min(real_cat[:, 4]))
    real_cat[:, 4] -= xmin
    real_cat[:, 5] -= xmin
    real_cat[:, 6] -= xmin

    red_cat = np.copy(real_cat)

    for i in range(real_cat.shape[0]):
        if scaling_method == "combined":
            vz = real_cat[i, 9] * gamma_l
            # print("Using combined (gamma_f) scaling, vz={}".format(vz))

        elif scaling_method == "split":
            # print("Using split (gamma_l/gamma_n) scaling")
            vz = (real_cat[i, 12] * gamma_l +
                  (real_cat[i, 9] - real_cat[i, 12]) * gamma_n)

        elif scaling_method == "transition":
            A = np.sqrt(((real_cat[i, 7] - real_cat[i, 10])**2. +
                         (real_cat[i, 8] - real_cat[i, 11])**2. +
                         (real_cat[i, 9] - real_cat[i, 12])**2.) /
                        (real_cat[i, 7]**2. +
                         real_cat[i, 8]**2. +
                         real_cat[i, 9]**2.))
            fA = (np.tanh((A - tc) / tb) + 1.) / 2.
            # print("b={}, c={}, A={}, fA={}, "
            #       "tanh((A-c)/b)={}, tanh(A-c)={}, tanh(A/b)={}, "
            #       "tanh(A)={}".format(b, c, A, fA, np.tanh((A - c) / b),
            #                           np.tanh(A - c), np.tanh(A/b),
            #                           np.tanh(A)))
            vz = real_cat[i, 9] * ((1. - fA) * gamma_l + fA * gamma_n)
            # print("Using transition (from gamma_l to gamma_n) scaling,"
            #       " A={}, fA={}, vz={}".format(A, fA, vz))

        else:
            vz = real_cat[i, 9]

        # Faizan's no integral method
        if zerr:
            zerr = np.random.normal(1.3, zerr_std)
            z_red = (real_cat[i, 6] + (vz + zerr) /
                     (Omega_M * (1. + zs_real)**3. + Omega_L)**0.5 / 100. *
                     (1. + zs_real))
            red_cat[i, 9] += zerr

        else:
            z_red = (real_cat[i, 6] + vz /
                     (Omega_M * (1. + zs_real)**3. + Omega_L)**0.5 / 100. *
                     (1. + zs_real))

        if z_red < 0.:
            z_red += boxsize

        elif z_red > boxsize:
            z_red -= boxsize

        if (i + 1) % 100000 == 0:
            print("zs_real:", zs_real, "vz:", vz, "z_real:",
                  real_cat[i, 6], "z_red:", z_red)
            print("%d Lines Complete:" % (i+1), dt.datetime.now() - start_time)

        red_cat[i, 6] = z_red

    fmt_list = ['%i', '%1.6e', '%f', '%f', '%f', '%f', '%f', '%f', '%f', '%f',
                '%f', '%f', '%f']

    if save_output:
        np.savetxt(output_root + "_redshift.dat", red_cat, fmt=fmt_list)
        print("Save Catalogue:", dt.datetime.now() - start_time)

    return real_cat, red_cat


def L_0(x):
    return 1.


def L_2(x):
    return 0.5 * (3. * x * x - 1.)


def L_4(x):
    return (35. * x * x * x * x - 30. * x * x + 3.) / 8.


L_polynomials = [L_0, L_2, L_4]


def corr_func_multipole_calc(ell, corr_func, mu_bins=100):
    L = L_polynomials[np.floor(ell / 2).astype(int)]

    corr_func_multipole = np.zeros((np.unique(corr_func[:, 0]).size, 2))
    corr_func_multipole[:, 0] = np.unique(corr_func[:, 0])
    for multi in corr_func_multipole:
        for corr in corr_func:
            if corr[0] == multi[0]:

                multi[1] += (1. / float(mu_bins)) * corr[4] * L(corr[3] - 0.5 /
                                                                float(mu_bins))

    corr_func_multipole[:, 1] = corr_func_multipole[:, 1] * (2. * ell + 1.)

    return corr_func_multipole


def calc_corrfuncs(red_cat, output_root, boxsize, save_output, real_cat=None,
                   N_bin_s=9):
    start_time = dt.datetime.now()
    print("Calculating Correlation Functions, Start Time:", start_time)

    N = red_cat.shape[0]

    # pimax = 160.
    pimax = 80.
    nthreads = 1
    bins = np.logspace(np.log10(0.1), np.log10(60.0), N_bin_s+1)

    mu_max = 1.
    nmu_bins = 100
    dmu = float(mu_max / nmu_bins)

    print("About to start wp")

    wp_results_red = wp(boxsize, pimax, nthreads, bins, red_cat[:, 4],
                        red_cat[:, 5], red_cat[:, 6])
    print("Red wp Calculated:", dt.datetime.now() - start_time)

    DDsmu_list_red = DDsmu(True, nthreads, bins, mu_max, nmu_bins,
                           red_cat[:, 4], red_cat[:, 5], red_cat[:, 6],
                           periodic=True, boxsize=boxsize)
    print("Red DDsmu Calculated:", dt.datetime.now() - start_time)

    DDsmu_results_red = np.empty((len(DDsmu_list_red), 6))
    for i in range(len(DDsmu_list_red)):
        DDsmu_results_red[i, 0] = DDsmu_list_red[i][0]
        DDsmu_results_red[i, 1] = DDsmu_list_red[i][1]
        DDsmu_results_red[i, 2] = DDsmu_list_red[i][2]
        DDsmu_results_red[i, 3] = DDsmu_list_red[i][3]
        DDsmu_results_red[i, 4] = DDsmu_list_red[i][4]
        DDsmu_results_red[i, 5] = DDsmu_list_red[i][5]

    RRsmu = np.copy(DDsmu_results_red)
    for i in range(RRsmu.shape[0]):
        RRsmu[i, 4] = (N * (N - 1.) / boxsize**3. * 4. / 3. * np.pi *
                       (RRsmu[i, 1]**3. - RRsmu[i, 0]**3.) * dmu)

    xi_red = np.copy(DDsmu_results_red)
    xi_red[:, 4] = DDsmu_results_red[:, 4] / RRsmu[:, 4] - 1.
    mono_red = corr_func_multipole_calc(0, xi_red)
    quad_red = corr_func_multipole_calc(2, xi_red)

    wp_meas_red = np.empty(len(wp_results_red))
    for i in range(len(wp_results_red)):
        wp_meas_red[i] = wp_results_red[i][3]

    bin_mins = np.copy(mono_red[:, 0])
    dv_red = np.empty((3*N_bin_s, 2))
    dv_red[:, 0] = np.concatenate((bin_mins, bin_mins, bin_mins))
    dv_red[:, 1] = np.concatenate((mono_red[:, 1], quad_red[:, 1],
                                   wp_meas_red))

    if save_output:
        np.savetxt(output_root + "_red_wp_results.dat", wp_results_red)
        np.savetxt(output_root + "_red_DDsmu.dat", DDsmu_results_red)
        np.savetxt(output_root + "_red_mono.dat", mono_red)
        np.savetxt(output_root + "_red_quad.dat", quad_red)
        np.savetxt(output_root + "_red_dv.dat", dv_red)
        np.savetxt(output_root + "_RRsmu.dat", RRsmu)

    if real_cat is not None:
        wp_results_real = wp(boxsize, pimax, nthreads, bins, real_cat[:, 4],
                             real_cat[:, 5], real_cat[:, 6])
        DDsmu_list_real = DDsmu(True, nthreads, bins, mu_max, nmu_bins,
                                real_cat[:, 4], real_cat[:, 5], real_cat[:, 6],
                                periodic=True, boxsize=boxsize)
        DDsmu_results_real = np.empty((len(DDsmu_list_real), 6))
        for i in range(len(DDsmu_list_real)):
            DDsmu_results_real[i, 0] = DDsmu_list_real[i][0]
            DDsmu_results_real[i, 1] = DDsmu_list_real[i][1]
            DDsmu_results_real[i, 2] = DDsmu_list_real[i][2]
            DDsmu_results_real[i, 3] = DDsmu_list_real[i][3]
            DDsmu_results_real[i, 4] = DDsmu_list_real[i][4]
            DDsmu_results_real[i, 5] = DDsmu_list_real[i][5]

        xi_real = np.copy(DDsmu_results_real)
        xi_real[:, 4] = DDsmu_results_real[:, 4] / RRsmu[:, 4] - 1.
        mono_real = corr_func_multipole_calc(0, xi_real)
        quad_real = corr_func_multipole_calc(2, xi_real)
        wp_meas_real = np.empty(len(wp_results_real))
        for i in range(len(wp_results_real)):
            wp_meas_real[i] = wp_results_real[i][3]

        dv_real = np.empty((3*N_bin_s, 2))
        dv_real[:, 0] = np.concatenate((bin_mins, bin_mins, bin_mins))
        dv_real[:, 1] = np.concatenate((mono_real[:, 1], quad_real[:, 1],
                                        wp_meas_real))

        if save_output:
            np.savetxt(output_root + "_real_wp_results.dat", wp_results_real)
            np.savetxt(output_root + "_real_DDsmu.dat", DDsmu_results_real)
            np.savetxt(output_root + "_real_mono.dat", mono_real)
            np.savetxt(output_root + "_real_quad.dat", quad_real)
            np.savetxt(output_root + "_real_dv.dat", dv_real)

    return dv_red, N


def cov_to_cor(cov_matrix):
    cor_matrix = np.copy(cov_matrix)

    for i in range(cor_matrix.shape[0]):
        for j in range(cor_matrix.shape[1]):
            cor_matrix[i, j] = cor_matrix[i, j] / np.sqrt(cov_matrix[i, i] *
                                                          cov_matrix[j, j])

    return cor_matrix


def cor_to_cov(cor_matrix, cov_mat_diags):
    cov_matrix = np.copy(cor_matrix)

    for i in range(cov_matrix.shape[0]):
        for j in range(cov_matrix.shape[1]):
            cov_matrix[i, j] = cov_matrix[i, j] * np.sqrt(cov_mat_diags[i] *
                                                          cov_mat_diags[j])

    return cov_matrix


def calc_cov_mat(path, dv, N, boxsize, save_output):
    print("N:", N)

    n = N / boxsize**3.
    print("n:", n)
    n_scaled = n * 1.e4
    print("n_scaled:", n_scaled)
    V_eff = (n_scaled / (1. + n_scaled))**2. * (boxsize / 1000.)**3.

    emulator_error_frac_mono = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                                          "emulator_error_fractional_mono.dat")
    emulator_error_frac_quad = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                                          "emulator_error_fractional_quad.dat")
    emulator_error_frac_wp = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                                        "emulator_error_fractional_wp.dat")
    emulator_error_frac = np.concatenate((emulator_error_frac_mono,
                                          emulator_error_frac_quad,
                                          emulator_error_frac_wp))
    emulator_error_mock = emulator_error_frac * np.abs(dv[:, 1])

    eboss_cm_diag = np.loadtxt("/home/mj3chapm/P2/fit/data/"
                               "data_jk_200_cov_mat_pip_ang_diag_sm.dat")
    mock_cm = eboss_cm_diag * 1.28 / V_eff

    data_error_mock = np.sqrt(np.diagonal(mock_cm))
    comb_error_mock = (np.sqrt(np.power(data_error_mock, 2.) +
                               np.power(emulator_error_mock, 2.)))

    mock_cor_mat = cov_to_cor(mock_cm)
    mock_comb_cov_mat = cor_to_cov(mock_cor_mat, np.power(comb_error_mock, 2.))

    if save_output:
        np.savetxt(path + "_red_cov_mat.dat", mock_comb_cov_mat)

    return mock_comb_cov_mat


def mock_meas(real_cat_path, output_root, zs_real, boxsize, save_output=True,
              zerr=False, zerr_std=91.8, calc_real_cat=False,
              scaling_method="combined", gamma_l=1., gamma_n=1., tb=0.2,
              tc=0.5):
    real_cat, red_cat = gen_red_cat(real_cat_path, output_root, zs_real,
                                    boxsize, save_output, zerr=zerr,
                                    zerr_std=zerr_std,
                                    scaling_method=scaling_method,
                                    gamma_l=gamma_l, gamma_n=gamma_n, tb=tb,
                                    tc=tc)
    if calc_real_cat:
        dv_red, N = calc_corrfuncs(red_cat, output_root, boxsize, save_output,
                                   real_cat=real_cat)
    else:
        dv_red, N = calc_corrfuncs(red_cat, output_root, boxsize, save_output)

    cov_mat = calc_cov_mat(output_root, dv_red, N, boxsize, save_output)
    return dv_red, cov_mat


# Change Log
# v0.2.3, 2021-05-10 - Updated zerr
# v0.2.1, 2021-03-10 - Updated to include save_output and return red_dv
