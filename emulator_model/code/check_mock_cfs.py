import numpy as np

simname = "AbacusCosmos_1100box_planck"
boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
          "00-5", "00-6", "00-7", "00-8", "00-9",
          "00-10", "00-11", "00-12", "00-13", "00-14",
          "00-15", "00-16", "00-17", "00-18", "00-19"]
halotype = "rockstar"
scale_variables = ["gamma_l", "gamma_n"]
N_grids = [550, 1375]
N_grid_def = 1100
smooth_type_def = "tophat"
smooth_types = ["tophat", "gaussian"]
R_smooths = [3., 7.]
R_smooth_def = 5.
scaling_method = "split"
gamma_val = 1.2

for i, scale_variable in enumerate(scale_variables):
    for j, boxid in enumerate(boxids):
        path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                     "{}_{}_products/{}_{}_{}_halos/z0.700/"
                     "em_model_test".format(simname, simname, boxid,
                                            simname, boxid, halotype))
        try:
            ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                ".dat".format(path_root))
        except FileNotFoundError:
            print("Missing Box {}, reference".format(boxid))

        run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                    "{}".format(scaling_method, N_grid_def,
                                smooth_type_def,
                                R_smooth_def, scale_variable))
        try:
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name, gamma_val))
        except FileNotFoundError:
            print("Missing Box {}, Fiducial, {}".format(boxid, scale_variable))

    for j in range(len(R_smooths)):
        for k, boxid in enumerate(boxids):
            path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                         "{}_{}_products/{}_{}_{}_halos/z0.700/"
                         "em_model_test".format(simname, simname, boxid,
                                                simname, boxid, halotype))
            try:
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
            except FileNotFoundError:
                print("Missing Box {}, reference".format(boxid))

            run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                        "{}".format(scaling_method, N_grid_def,
                                    smooth_type_def,
                                    R_smooths[j], scale_variable))
            try:
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))
            except FileNotFoundError:
                print("Missing Box {}, R={}, {}".format(boxid, R_smooths[j],
                                                        scale_variable))

    for j in range(len(N_grids)):
        grid_length = "{:.1f}".format(1100 / N_grids[j])
        for k, boxid in enumerate(boxids):
            path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                         "{}_{}_products/{}_{}_{}_halos/z0.700/"
                         "em_model_test".format(simname, simname, boxid,
                                                simname, boxid, halotype))
            try:
                ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                    ".dat".format(path_root))
            except FileNotFoundError:
                print("Missing Box {}, reference".format(boxid))

            run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                        "{}".format(scaling_method, N_grids[j],
                                    smooth_type_def, R_smooth_def,
                                    scale_variable))
            try:
                red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                    ".dat".format(path_root, run_name,
                                                  gamma_val))
            except FileNotFoundError:
                print("Missing Box {}, N={}, {}".format(boxid, N_grids[j],
                                                        scale_variable))

    for j, boxid in enumerate(boxids):
        path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                     "{}_{}_products/{}_{}_{}_halos/z0.700/"
                     "em_model_test".format(simname, simname, boxid,
                                            simname, boxid, halotype))
        try:
            ref_dv = np.loadtxt("{}/combined_v2_1.0_red_dv"
                                ".dat".format(path_root))
        except FileNotFoundError:
            print("Missing Box {}, reference".format(boxid))

        run_name = ("{}_{}_{}_{}_smoothed_v2_all-halos_"
                    "{}".format(scaling_method, N_grid_def,
                                "gaussian",
                                2., scale_variable))
        try:
            red_dv = np.loadtxt("{}/{}_{}_red_dv"
                                ".dat".format(path_root, run_name,
                                              gamma_val))
        except FileNotFoundError:
            print("Missing Box {}, gaussian, {}".format(boxid,
                                                        scale_variable))
