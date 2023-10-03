import datetime as dt
import numpy as np
import os
from struct import unpack_from

start_time = dt.datetime.now()
print("Starting at", start_time)

# Prepare array to hold the particle velocities
# EDIT FROM ORIG
subsample_array = np.loadtxt("test_cats/mpi_test_subsamples.dat")

pids = subsample_array[:, 0]
N_pids = pids.size

# Sort the particle IDs so only one loop of the particles is needed
subsample_array = subsample_array[np.argsort(pids)]
sorted_pids = subsample_array[:, 0]

sim_name = "AbacusCosmos_1100box_planck"
boxid = "00-0"
z_ic = 49.0
ic_type = "ngc_nplt"

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))
ic_dir = ("{}/ic_{}_z{}".format(path_root, ic_type, z_ic))
N_ic_files = 10
# N_ic_files = len(os.listdir(ic_dir)) # EDIT_FROM_ORIG

last_sort_pid_index = 0
tot_particles = 0

for i in range(N_ic_files):
    print("Starting IC file {}, elapsed time".format(i),
          dt.datetime.now() - start_time)
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

                subsample_array[j, 4:] = particle_data[-3:]

            else:
                print("Switched IC at {}".format(j))
                last_sort_pid_index = j
                break

    tot_particles += N_particles

print("Particles caclulated, elapsed time", dt.datetime.now() - start_time)
np.savetxt("test_cats/no_mpi_subsample_array.dat", subsample_array)
print("Code finished, elapsed time", dt.datetime.now() - start_time)
