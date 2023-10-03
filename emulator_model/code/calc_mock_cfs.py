# Redshift error correlation function check
# v0.2.0, 2022-02-04 - Updated for efficiency, removed transition method, added
#                      smoothed option

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
load_dvs = sys.argv[6]
smoothed_vel = sys.argv[7]
N_grid = int(sys.argv[8])
smooth_type = sys.argv[9]
R_smooth = float(sys.argv[10])

load_dvs = True if load_dvs == "True" else False
smoothed_vel = True if smoothed_vel == "True" else False

redshift = 0.70
Lbox = 1100.
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
             "{}_{}_{}_halos/z0.700" .format(simname, simname, boxid, simname,
                                             boxid, halotype))
Path("{}/em_model_test".format(path_root)).mkdir(parents=False, exist_ok=True)
if smoothed_vel:
    real_cat_path = ("{}/halo_lin_vel_{}_{}_{}_smoothed_v2"
                     ".dat".format(path_root, N_grid, smooth_type, R_smooth))
    run_name = "{}_{}_{}_{}_smoothed_v2".format(scaling_method, N_grid,
                                                smooth_type, R_smooth)
else:
    real_cat_path = ("{}/halo_lin_vel.dat".format(path_root))
    run_name = scaling_method

if scaling_method == "split":
    run_name += "_{}".format(scale_variable)

gamma_val = 1.16

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
# v0.1.1, 2022-02-04 - Copied from RSD to P2 for testing new emulator models
# v0.1.0, 2021-05-11 - Code started
