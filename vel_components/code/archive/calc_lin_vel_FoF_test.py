# Calculate the linear velocity component of halos
# v0.3.2, 2021-11-17 - Testing FoF run problems

# Imports
import datetime as dt
import numpy as np
import os
from struct import unpack_from
import sys

from AbacusCosmos import Halos
from lin_vel_funcs import calc_vel_scaling


start_time = dt.datetime.now()

print("Starting preliminaries at", start_time)

# Get user input parameters
# boxid - ID of the simulation box being used
# halo_type - rockstar or FoF
boxid = "planck"
# Here for testing, but should be fixed for convenience in the future
halo_type = "FoF"

# Important global constants, that could be changed to parameters in the future
z_cat = 0.7
z_ic = 49.0
ic_type = "ngc_nplt"

# Load the cosmological parameters of the simulation box
info_path = ("/home/mj3chapm/scratch/abacus/"
             "AbacusCosmos_1100box_products/"
             "AbacusCosmos_1100box_{}_products/"
             "AbacusCosmos_1100box_{}_rockstar_halos/"
             "info".format(boxid, boxid))

# Load the cosmological parameters of the simulation box
cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_{}_products/"
                          "AbacusCosmos_1100box_{}_rockstar_halos/"
                          "info/cosmo_params.dat".format(boxid, boxid))

# Calculate the scaling that needs to be applied to the initial condition
# particle velocities
vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)

print("Loading halos, elapsed time", dt.datetime.now() - start_time)

# Load the halo catalogue
cat = Halos.make_catalog_from_dir(dirname="/home/mj3chapm/scratch/abacus/"
                                          "AbacusCosmos_1100box_products/"
                                          "AbacusCosmos_1100box_{}_products/"
                                          "FoF_test/".format(boxid),
                                  load_subsamples=True, load_pids=True)
halos = cat.halos

# Prepare array to hold the particle velocities
subsamples = cat.subsamples

if halo_type == "FoF":
    pids = subsamples["pid"]

elif halo_type == "rockstar":
    pids = np.array([part[0] for part in subsamples])

subsample_array = np.empty((pids.size, 4))
subsample_array[:, 0] = pids

# Sort the particle IDs so only one loop of the particles is needed
sort_map = np.argsort(pids)
sorted_pids = pids[sort_map]
last_sort_pid_index = 0
tot_particles = 0

ic_dir = ("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
          "AbacusCosmos_1100box_{}_products/ic_{}_z{}".format(boxid, ic_type,
                                                              z_ic))

print("Assembling initial conditions, elapsed time",
      dt.datetime.now() - start_time)

for i in range(375):
    print("Starting IC file {}".format(i))
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    # Total size of objects is 30, but becomes 32 because of padding
    N_particles = tot_ic_size / 32

    if i <= 31 or i >= 368:
        with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
            bdata = file.read()

            for j in range(last_sort_pid_index, pids.size):
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
    print("last_sort_pid_index", last_sort_pid_index)
    print("tot_particles", tot_particles)

print("Setting up halo array, elapsed time", dt.datetime.now() - start_time)

if halo_type == "FoF":
    halos = halos[halos["subsamp_len"] > 0]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["N"] * 4.e10 / (cosmo_params[0] / 100.)
    halo_array[:, 2] = halos["r90"]
    halo_array[:, 3] = halos["vcirc_max"]

elif halo_type == "rockstar":
    halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
    halo_array = np.empty((halos.shape[0], 13))
    halo_array[:, 1] = halos["m"]
    halo_array[:, 2] = halos["r"]
    halo_array[:, 3] = halos["vmax"]

print("Populating halo array, elapsed time", dt.datetime.now() - start_time)

for i in range(halos.shape[0]):
    # Abacus halo IDs are too long for the emulator code
    halo_array[:, 0] = 1
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

    if ((i + 1) % 10000) == 0:
        print("Finished {} halos, elapsed time".format(i+1),
              dt.datetime.now() - start_time)

print("Saving results, elapsed time", dt.datetime.now() - start_time)

# Columns of the resulting catalogue are shown in this header, which has been
# removed for consistency with the existing emulator code
# header = ("#ID M200b R200b Vmax X Y Z VX VY VZ VX_LIN VY_LIN VZ_LIN")
np.savetxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
           "AbacusCosmos_1100box_{}_products/"
           "AbacusCosmos_1100box_{}_{}_halos/"
           "z0.700/halo_lin_vel_FoF_test.dat".format(boxid, boxid, halo_type),
           halo_array, fmt=["%i", "%.4e", "%.4f", "%.4f", "%.3f", "%.3f",
                            "%.3f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                            "%.6f"])

print("Finished and everthying went great! Elapsed time",
      dt.datetime.now() - start_time)

# Change Log
# v0.2.3, 2021-11-05 - Changed ID and position fields to match emulator code
# v0.2.2, 2021-11-05 - Changed save format to reduce file size and removed
#                      header
# v0.2.1, 2021-10-25 - Updated to be PEP8 compliant and more streamlined
# v0.1.1, 2021-08-12 - Updated the velocity scaling
# v0.1.0, 2021-07-20 - Code started
