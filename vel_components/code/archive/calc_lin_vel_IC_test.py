# Calculate the linear velocity component of halos
# v0.3.3, 2021-11-19 - Testing IC run problems

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
boxid = "02"
# Here for testing, but should be fixed for convenience in the future
halo_type = "rockstar"

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
                                          "IC_test/".format(boxid),
                                  load_subsamples=True, load_pids=True,
                                  halo_type="Rockstar")
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

N_ic_files = len(os.listdir(ic_dir))


i = 291
tot_particles = 2979763200.0
last_sort_pid_index = 2943814
print("pids Size:", pids.size)
print("Starting IC file {}".format(i))
tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
# Total size of objects is 30, but becomes 32 because of padding
N_particles = tot_ic_size / 32

with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
    bdata = file.read()

    for j in range(last_sort_pid_index, pids.size):
        if j == 2950480:
            print("j: {}, sorted_pids[j]: {},".format(j, sorted_pids[j]),
                  "tot_particles: {}, N_particles: {},". format(tot_particles,
                                                                N_particles),
                  "Index: {},".format((sorted_pids[j] - tot_particles)*32),
                  "Condition:", int(sorted_pids[j] - tot_particles) <
                  N_particles)
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

# for i in range(N_ic_files):
#     print("Starting IC file {}".format(i))
#     tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
#     # Total size of objects is 30, but becomes 32 because of padding
#     N_particles = tot_ic_size / 32

#     # if i <= 31 or i >= 368:
#     with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
#         bdata = file.read()

#         for j in range(last_sort_pid_index, pids.size):
#             if int(sorted_pids[j] - tot_particles) < N_particles:
#                 particle_data = unpack_from("3h6f", bdata,
#                                             offset=int((sorted_pids[j] -
#                                                         tot_particles)*32))

#                 subsample_array[sort_map[j], 1:4] = particle_data[-3:]

#             else:
#                 print("Switched IC at {}".format(j))
#                 last_sort_pid_index = j
#                 break

#     tot_particles += N_particles
#     print("last_sort_pid_index", last_sort_pid_index)
#     print("tot_particles", tot_particles)

print("Finished and everthying went great! Elapsed time",
      dt.datetime.now() - start_time)

# Change Log
# v0.2.3, 2021-11-05 - Changed ID and position fields to match emulator code
# v0.2.2, 2021-11-05 - Changed save format to reduce file size and removed
#                      header
# v0.2.1, 2021-10-25 - Updated to be PEP8 compliant and more streamlined
# v0.1.1, 2021-08-12 - Updated the velocity scaling
# v0.1.0, 2021-07-20 - Code started
