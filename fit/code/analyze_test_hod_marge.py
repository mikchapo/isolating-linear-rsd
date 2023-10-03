"""
Use the marge stats data to calculate how parameters fit within limits.

v0.1.0, 2023-01-31 - Code started from margestats
"""

# Imports
import matplotlib.pyplot as plt
import numpy as np

param_labels = [r"H_0", r"\alpha", r"c_{vir}", r"f_{max}", r"\gamma_l",
                r"\gamma_n", r"logM_{cut}", r"logM_{sat}", r"n_\mathrm{s}",
                r"\Omega_\mathrm{b} h^2", r"\Omega_\mathrm{m}", r"\sigma_8",
                r"\sigma_{8,a}", r"\sigma_{logM}", r"v_{bc}", r"v_{bs}",
                r"f\sigma_8"]
param_to_plot_indices = [3, 6, 12, 9, 14, 13, 7, 5, 4, 1, 0, 2, 8, 10, 11, 15]
plot_to_param_indices = [10, 9, 11, 0, 8, 7, 1, 6, 13, 3, 14, 15, 2, 5, 4, 16]

cosmology = np.array([0.314153, 0.02222, 0.83, 67.26, 0.9652])
hod_models = np.loadtxt("/home/mj3chapm/P2/fit/data/test_clustering/"
                        "HOD_parameters_test.dat")
fsig8_ref = 0.4605572945563802
path_root = ("/home/mj3chapm/projects/rrg-wperciva/mj3chapm/P2/chains/"
             "smoothed_emulator/test_clustering")

N_params = len(param_labels)

N_rows = 4
N_cols = 4
fig, axes = plt.subplots(N_rows, N_cols, figsize=(10., 10.), dpi=300)

bins = np.zeros((N_rows, N_cols, 10))

for i in range(10):
# for i in range(1):
    run_name = ("{}/test_hod_{}_vol-factor-20_default/"
                "test_hod_{}_vol-factor-20_default".format(path_root, i, i))
    marge_params = np.genfromtxt("{}.margestats".format(run_name),
                                 skip_header=3,
                                 usecols=(1, 2, 3, 4, 6, 7, 9, 10))
    hod_params = np.empty(11)
    hod_params[:4] = hod_models[i, :4]
    hod_params[0] = np.log10(hod_params[0])
    hod_params[2] = np.log10(hod_params[2])
    hod_params[4] = hod_models[i, 8]
    hod_params[5:8] = hod_models[i, 4:7]
    hod_params[8] = hod_models[i, 7]
    hod_params[9] = hod_models[i, 9]
    hod_params[10] = hod_params[9] * fsig8_ref
    expected = np.concatenate((cosmology, hod_params))
    for j in range(N_rows):
        for k in range(N_cols):
            index = j * N_rows + k

            print("index:", index)
            print("param:", param_labels[plot_to_param_indices[index]])
            print("expected:", expected[index])
            print("mean:", marge_params[plot_to_param_indices[index], 0])
            print("limits:", marge_params[plot_to_param_indices[index], 2:])
            print()

            if expected[index] >= marge_params[plot_to_param_indices[index], 0]:
                if expected[index] >= marge_params[plot_to_param_indices[index], 5]:
                    if expected[index] >= marge_params[plot_to_param_indices[index], 7]:
                        bins[j, k, i] = 7
                    else:
                        bins[j, k, i] = 6

                else:
                    if expected[index] >= marge_params[plot_to_param_indices[index], 3]:
                        bins[j, k, i] = 5
                    else:
                        bins[j, k, i] = 4

            else:
                if expected[index] < marge_params[plot_to_param_indices[index], 4]:
                    if expected[index] < marge_params[plot_to_param_indices[index], 6]:
                        bins[j, k, i] = 0
                    else:
                        bins[j, k, i] = 1

                else:
                    if expected[index] < marge_params[plot_to_param_indices[index], 2]:
                        bins[j, k, i] = 2
                    else:
                        bins[j, k, i] += 3

            print(bins)

print(bins[0, 0, :])

for j in range(N_rows):
    for k in range(N_cols):
        axes[j, k].set_xlim(0, 7)
        axes[j, k].hist(bins[j, k, :])
        axes[j, k].axvline(3.5, 0, 1)
        index = j * N_rows + k
        axes[j, k].set_title(param_labels[plot_to_param_indices[index]])


plt.savefig("/home/mj3chapm/P2/fit/output/plots/analyze_chains/"
            "smoothed_emulator/test_hod/test_hod_marge.jpg")
