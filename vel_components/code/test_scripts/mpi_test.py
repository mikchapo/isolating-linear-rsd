import datetime as dt
from mpi4py import MPI
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

Np_files = np.empty(N_ic_files)
for i in range(N_ic_files):
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    Np_files[i] = tot_ic_size / 32

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

N_files_per_rank = N_ic_files // size
lower_bound = rank * N_files_per_rank
upper_bound = (rank + 1) * N_files_per_rank
if rank == (size - 1):
    upper_bound = N_ic_files

lower_pid = np.sum(Np_files[:lower_bound])
upper_pid = np.sum(Np_files[:upper_bound])

lower_pid_index = np.sum((lower_pid > pids))
N_subsample = np.sum((lower_pid <= pids) & (upper_pid > pids))
lin_vel_array = np.empty((N_subsample, 3))

last_sort_pid_index = lower_pid_index
tot_particles = lower_pid

for i in range(lower_bound, upper_bound):
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

        for j in range(last_sort_pid_index, lower_pid_index + N_subsample):
            if int(sorted_pids[j] - tot_particles) < N_particles:
                particle_data = unpack_from("3h6f", bdata,
                                            offset=int((sorted_pids[j] -
                                                        tot_particles)*32))

                lin_vel_array[j - lower_pid_index, :] = particle_data[-3:]

            else:
                print("Switched IC at {}".format(j))
                last_sort_pid_index = j
                break

    tot_particles += N_particles

comm.Barrier()
lin_vel_arr_list = comm.gather(lin_vel_array, root=0)

print("\n", "Printing first 10 before broadcast, rank {}".format(rank))
print(subsample_array[:10], "\n")

if rank == 0:
    subsample_array[:, 4:] = np.concatenate(lin_vel_arr_list, axis=0)
    print("Particles caclulated, elapsed time", dt.datetime.now() - start_time)
    np.savetxt("test_cats/mpi_subsample_array.dat", subsample_array)
    print("Code finished, elapsed time", dt.datetime.now() - start_time)

comm.Barrier()
comm.Bcast(subsample_array, root=0)

print("\n", "Printing first 10 after broadcast, rank {}".format(rank))
print(subsample_array[:10], "\n")
