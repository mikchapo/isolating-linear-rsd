# Assigned smoothed linear velocities to halos
# v0.2.0, 2022-08-31 - Added check to make sure halo zeroing was happening
#                      correctly

# Imports
import datetime as dt
import numpy as np
import sys

from AbacusCosmos import Halos

start_time = dt.datetime.now()
print("Starting at", start_time)

# Get user input parameters
# boxid - ID of the simulation box being used
# halo_type - rockstar or FoF
sim_name = sys.argv[1]
boxid = sys.argv[2]
halo_type = sys.argv[3]
N_grid = int(sys.argv[4])
smooth_type = sys.argv[5]
R_smooth = float(sys.argv[6])
subhalos = sys.argv[7] == "True"

boxsize = 1100.


path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

# Load the halo catalogue
cat = Halos.make_catalog_from_dir(dirname="{}/{}_{}_{}_halos/"
                                          "z0.700".format(path_root, sim_name,
                                                          boxid, halo_type),
                                  load_subsamples=False, load_pids=False)
halos = cat.halos
del cat

print("Halos loaded, elapsed time", dt.datetime.now() - start_time)

if halo_type == "FoF":
    halos = halos[halos["subsamp_len"] > 0]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["N"] * 4.e10
    halo_array[:, 2] = halos["r90"]
    halo_array[:, 3] = halos["vcirc_max"]

# elif (sim_name == "AbacusCosmos_1100box" and boxid == "planck" and
#       halo_type == "rockstar"):
#     halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
#     halo_array = np.empty((halos.shape[0], 13))
#     halo_array[:, 1] = halos["N"] * 4.e10
#     halo_array[:, 2] = halos["r"]
#     halo_array[:, 3] = halos["vmax"]

elif halo_type == "rockstar":
    if subhalos:
        halos = halos[halos["subsamp_len"] > 0]
    else:
        halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["m"]
    halo_array[:, 2] = halos["r"]
    halo_array[:, 3] = halos["vmax"]

xmin = np.floor(np.min(halos["pos"][:, 0]))
ymin = np.floor(np.min(halos["pos"][:, 1]))
zmin = np.floor(np.min(halos["pos"][:, 2]))

print("Minimum positions are:", xmin, ymin, zmin)

print("Halo array created, elapsed time", dt.datetime.now() - start_time)

vel_grid = np.load("{}/vel_grid_{}_{}_{}_smoothed_v3"
                   ".npy".format(path_root, N_grid, smooth_type, R_smooth))
vel_grid = np.nan_to_num(vel_grid, copy=False)

l_grid = boxsize / N_grid

print("Velocity grid loaded, elapsed time", dt.datetime.now() - start_time)

for i in range(halos.shape[0]):
    # Abacus halo IDs are too long for the emulator code
    halo_array[i, 0] = i
    # Default box limits are [-550, +550], so shift to always be positive
    halo_array[i, 4] = halos[i]["pos"][0] - xmin
    halo_array[i, 5] = halos[i]["pos"][1] - ymin
    halo_array[i, 6] = halos[i]["pos"][2] - zmin
    halo_array[i, 7] = halos[i]["vel"][0]
    halo_array[i, 8] = halos[i]["vel"][1]
    halo_array[i, 9] = halos[i]["vel"][2]
    try:
        halo_array[i, 10] = vel_grid[int(halo_array[i, 4] // l_grid),
                                     int(halo_array[i, 5] // l_grid),
                                     int(halo_array[i, 6] // l_grid), 0]
        halo_array[i, 11] = vel_grid[int(halo_array[i, 4] // l_grid),
                                     int(halo_array[i, 5] // l_grid),
                                     int(halo_array[i, 6] // l_grid), 1]
        halo_array[i, 12] = vel_grid[int(halo_array[i, 4] // l_grid),
                                     int(halo_array[i, 5] // l_grid),
                                     int(halo_array[i, 6] // l_grid), 2]
    except IndexError:
        print("IndexError for halo {} with position ({}, {}, {}"
              ")".format(i, halo_array[i, 4], halo_array[i, 5],
                         halo_array[i, 6]))
        x_index = int(halo_array[i, 4] // l_grid)
        y_index = int(halo_array[i, 5] // l_grid)
        z_index = int(halo_array[i, 6] // l_grid)
        if x_index >= N_grid:
            x_index = -1
        if y_index >= N_grid:
            y_index = -1
        if z_index >= N_grid:
            z_index = -1
        halo_array[i, 10] = vel_grid[x_index, y_index, z_index, 0]
        halo_array[i, 11] = vel_grid[x_index, y_index, z_index, 1]
        halo_array[i, 12] = vel_grid[x_index, y_index, z_index, 2]

    if ((i + 1) % 100000) == 0:
        print("Completed {} halos, elapsed time".format(i + 1),
              dt.datetime.now() - start_time)

# Columns of the resulting catalogue are shown in this header, which has been
# removed for consistency with the existing emulator code
# header = ("#ID M200b R200b Vmax X Y Z VX VY VZ VX_LIN VY_LIN VZ_LIN")
if subhalos:
    filename = ("{}/{}_{}_{}_halos/z0.700/"
                "halo_lin_vel_{}_{}_{}_smoothed_v3_all-halos"
                ".dat".format(path_root, sim_name, boxid, halo_type, N_grid,
                              smooth_type, R_smooth))

else:
    filename = ("{}/{}_{}_{}_halos/z0.700/"
                "halo_lin_vel_{}_{}_{}_smoothed_v3"
                ".dat".format(path_root, sim_name, boxid, halo_type, N_grid,
                              smooth_type, R_smooth))

np.savetxt(filename, halo_array,
           fmt=["%i", "%.4e", "%.4f", "%.4f", "%.3f", "%.3f", "%.3f", "%.6f",
                "%.6f", "%.6f", "%.6f", "%.6f", "%.6f"])

# Output to see if positions have been properly zeroed
if boxid == "00-0" and (sim_name == "AbacusCosmos_1100box_planck" and
                        N_grid == 1100 and smooth_type == "tophat" and
                        R_smooth == 5.):
    # sample = halo_array[halo_array[:, 4] >= 49.]
    # sample = sample[sample[:, 4] < 50.]
    # sample = sample[sample[:, 5] >= 1091.]
    # sample = sample[sample[:, 5] < 1092.]
    # sample = sample[sample[:, 6] >= 1020.]
    # sample = sample[sample[:, 6] < 1021.]
    # print("Test cell for proper position range (49, 1091, 1020):")
    # print("Test halos within cell:")
    # print(sample)
    # print("Expected smoothed linear velocities: "
    #       "[160.86853854, 205.63755677, 273.31266465]")
    # expected = np.array([160.86853854, 205.63755677, 273.31266465, 10.])
    # for i in range(sample.shape[0]):
    #     print("Differences:", sample[i, 10:13] - expected)

    print("Test cell for proper position range (811, 573, 807):")
    print("Halo that should fit:", halo_array[0, :])

# Change Log
# v0.1.2, 2022-08-02 - Added option to include subhalos
# v0.1.1, 2022-06-02 - Fixed bug where positions were between 550 and 1650
# v0.1.0, 2022-04-22 - Code started with snippets from calc_lin_vel_smoothed.py
