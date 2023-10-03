# Redshift error correlation function check
# v0.3.1, 2022-09-13 - Added gamma=0.8 calculation

# Import
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import seaborn as sns
import sys

from mock_meas import mock_meas


# User input parameters
# simname - [AbacusCosmos_1100box]
# boxid - [00, 01, 02, 03,...|00-0, 00-1, 00-2,...]
# halotype - [rockstar|FoF]
# scaling_method - [combined|split]
# scale_variable - [gamma_f|gamma_l|gamma_n]
# load_dvs - [True|False]
# smoothed_vel - [True|False]
simname = sys.argv[1]
boxid = sys.argv[2]
halotype = sys.argv[3]
scaling_method = sys.argv[4]
scale_variable = sys.argv[5]
smoothed_vel = sys.argv[6] == "True"
N_grid = int(sys.argv[7])
smooth_type = sys.argv[8]
R_smooth = float(sys.argv[9])
subhalos = sys.argv[10] == "True"

redshift = 0.70
Lbox = 1100.
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
             "{}_{}_{}_halos/z0.700" .format(simname, simname, boxid, simname,
                                             boxid, halotype))
Path("{}/em_model_test".format(path_root)).mkdir(parents=False, exist_ok=True)
if smoothed_vel:
    if subhalos:
        real_cat_path = ("{}/halo_lin_vel_{}_{}_{}_smoothed_v2_all-halos"
                         ".dat".format(path_root, N_grid, smooth_type,
                                       R_smooth))
        run_name = "{}_{}_{}_{}_smoothed_v2_all-halos".format(scaling_method,
                                                              N_grid,
                                                              smooth_type,
                                                              R_smooth)
    else:
        real_cat_path = ("{}/halo_lin_vel_{}_{}_{}_smoothed_v2"
                         ".dat".format(path_root, N_grid, smooth_type,
                                       R_smooth))
        run_name = "{}_{}_{}_{}_smoothed_v2".format(scaling_method, N_grid,
                                                    smooth_type, R_smooth)
else:
    if subhalos:
        real_cat_path = ("{}/halo_lin_vel_all-halos.dat".format(path_root))
        run_name = scaling_method + "_all-halos"
    else:
        real_cat_path = ("{}/halo_lin_vel.dat".format(path_root))
        run_name = scaling_method

if scaling_method == "split":
    run_name += "_{}".format(scale_variable)

bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))

output_root = "{}/em_model_test/combined_v2_1.0".format(path_root)
ref_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift, Lbox,
                            scaling_method=scaling_method,
                            gamma_l=1.0)

if subhalos:
    gamma_vals = [0.8, 1.2]
    # gamma_vals = [0.8]
else:
    gamma_vals = [0.84, 0.92, 1.08, 1.16]

for i, gamma_val in enumerate(gamma_vals):
    output_root = ("{}/em_model_test/{}_{}".format(path_root, run_name,
                                                   gamma_val))
    if scale_variable == "gamma_l":
        red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                    Lbox, scaling_method=scaling_method,
                                    gamma_l=gamma_val)

    elif scale_variable == "gamma_n":
        red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                    Lbox, scaling_method=scaling_method,
                                    gamma_n=gamma_val)

    else:
        red_dv, cov_mat = mock_meas(real_cat_path, output_root, redshift,
                                    Lbox, scaling_method=scaling_method,
                                    gamma_l=gamma_val)

# Change Log
# v0.3.0, 2022-08-02 - Changed for test boxes and paper plot
# v0.2.0, 2022-02-04 - Updated for efficiency, removed transition method, added
#                      smoothed option
# v0.1.1, 2022-02-04 - Copied from RSD to P2 for testing new emulator models
# v0.1.0, 2021-05-11 - Code started
