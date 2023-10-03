# Calculate the linear velocity component of halos
# v0.4.1, 2022-08-10 - Added option to include subhalos

# Imports
import numpy as np
import os
from struct import unpack_from
import sys

from AbacusCosmos import Halos
from lin_vel_funcs import calc_vel_scaling


# Get user input parameters
# boxid - ID of the simulation box being used
# halo_type - rockstar or FoF
boxid = sys.argv[1]
# Here for testing, but should be fixed for convenience in the future
halo_type = sys.argv[2]
sim_name = sys.argv[3]
subhalos = sys.argv[4] == "True"

# Important global constants, that could be changed to parameters in the future
z_cat = 0.7
z_ic = 49.0
ic_type = "ngc_nplt"

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

# Load the cosmological parameters of the simulation box
info_path = ("{}/{}_{}_rockstar_halos/"
             "info".format(path_root, sim_name, boxid))

# Load the cosmological parameters of the simulation box
cosmo_params = np.loadtxt("{}/{}_{}_rockstar_halos/"
                          "info/cosmo_params.dat".format(path_root, sim_name,
                                                         boxid))

# Calculate the scaling that needs to be applied to the initial condition
# particle velocities
vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)

# Load the halo catalogue
cat = Halos.make_catalog_from_dir(dirname="{}/{}_{}_{}_halos/"
                                          "z0.700".format(path_root, sim_name,
                                                          boxid, halo_type),
                                  load_subsamples=True, load_pids=True)
halos = cat.halos

# Prepare array to hold the particle velocities
subsamples = cat.subsamples

if halo_type == "FoF":
    pids = subsamples["pid"]

elif halo_type == "rockstar":
    pids = np.array([part[0] for part in subsamples])

N_pids = pids.size
subsample_array = np.empty((N_pids, 4))
subsample_array[:, 0] = pids

# Sort the particle IDs so only one loop of the particles is needed
sort_map = np.argsort(pids)
sorted_pids = pids[sort_map]
last_sort_pid_index = 0
tot_particles = 0

ic_dir = ("{}/ic_{}_z{}".format(path_root, ic_type, z_ic))

N_ic_files = len(os.listdir(ic_dir))
for i in range(N_ic_files):
    print("Starting IC file {}".format(i))
    with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
        bdata = file.read()

        particle_data_first = unpack_from("3h6f", bdata, offset=0)
        particle_data_last = unpack_from("3h6f", bdata, offset=-32)
        N_particles = ((particle_data_last[0] -
                        particle_data_first[0]) * 1440**2 +
                       (particle_data_last[1] -
                        particle_data_first[1]) * 1440 +
                       (particle_data_last[2] -
                        particle_data_first[2]) + 1)

        for j in range(last_sort_pid_index, N_pids):
            if int(sorted_pids[j] - tot_particles) < N_particles:
                particle_data = unpack_from("3h6f", bdata,
                                            offset=int((sorted_pids[j] -
                                                        tot_particles)*32))

                subsample_array[sort_map[j], 1:4] = particle_data[-3:]

            else:
                print("Switched IC at {}".format(j))
                last_sort_pid_index = j
                break

    tot_particles += N_particles

if halo_type == "FoF":
    halos = halos[halos["subsamp_len"] > 0]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["N"] * 4.e10
    halo_array[:, 2] = halos["r90"]
    halo_array[:, 3] = halos["vcirc_max"]

elif (sim_name == "AbacusCosmos_1100box" and boxid == "planck" and
      halo_type == "rockstar"):
    halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["N"] * 4.e10
    halo_array[:, 2] = halos["r"]
    halo_array[:, 3] = halos["vmax"]

elif halo_type == "rockstar":
    if subhalos:
        halos = halos[halos["subsamp_len"] > 0]
    else:
        halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["m"]
    halo_array[:, 2] = halos["r"]
    halo_array[:, 3] = halos["vmax"]

else:
    raise ValueError("Unrecognized halo type!")

for i in range(halos.shape[0]):
    # Abacus halo IDs are too long for the emulator code
    halo_array[i, 0] = i
    # Default box limits are [-550, +550], so shift to always be positive
    halo_array[i, 4] = halos[i]["pos"][0] + 550.
    halo_array[i, 5] = halos[i]["pos"][1] + 550.
    halo_array[i, 6] = halos[i]["pos"][2] + 550.
    halo_array[i, 7] = halos[i]["vel"][0]
    halo_array[i, 8] = halos[i]["vel"][1]
    halo_array[i, 9] = halos[i]["vel"][2]
    subsamp_start = halos[i]["subsamp_start"]
    subsamp_end = halos[i]["subsamp_start"]+halos[i]["subsamp_len"]
    halo_array[i, 10] = np.mean(subsample_array[subsamp_start:subsamp_end,
                                1]) * vel_scaling
    halo_array[i, 11] = np.mean(subsample_array[subsamp_start:subsamp_end,
                                2]) * vel_scaling
    halo_array[i, 12] = np.mean(subsample_array[subsamp_start:subsamp_end,
                                3]) * vel_scaling

# Columns of the resulting catalogue are shown in this header, which has been
# removed for consistency with the existing emulator code
# header = ("#ID M200b R200b Vmax X Y Z VX VY VZ VX_LIN VY_LIN VZ_LIN")
if subhalos:
    filename = ("{}/{}_{}_{}_halos/z0.700/halo_lin_vel_all-halos"
                ".dat".format(path_root, sim_name, boxid, halo_type))

else:
    filename = ("{}/{}_{}_{}_halos/z0.700/halo_lin_vel"
                ".dat".format(path_root, sim_name, boxid, halo_type))

np.savetxt(filename, halo_array,
           fmt=["%i", "%.4e", "%.4f", "%.4f", "%.3f", "%.3f", "%.3f", "%.6f",
                "%.6f", "%.6f", "%.6f", "%.6f", "%.6f"])


# Change Log
# v0.4.0, 2021-12-03 - Added other simulation types and changed FoF back to
#                      M_sun/h mass units
# v0.3.2, 2021-11-18 - Fixed Rockstar mask, bad ID assignment, and harded
#                      number of initial conditions files
# v0.3.1, 2021-11-11 - Updated for new calc_vel_scaling version and Rockstar
#                      halos
# v0.2.3, 2021-11-05 - Changed ID and position fields to match emulator code
# v0.2.2, 2021-11-05 - Changed save format to reduce file size and removed
#                      header
# v0.2.1, 2021-10-25 - Updated to be PEP8 compliant and more streamlined
# v0.1.1, 2021-08-12 - Updated the velocity scaling
# v0.1.0, 2021-07-20 - Code started
