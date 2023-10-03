import numpy as np
import os

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

print(Np_files)

N_subsample = 10
subsample_array = np.zeros((N_ic_files * N_subsample, 7))

for i in range(N_ic_files):
    lower_pid = np.sum(Np_files[:i])
    upper_pid = np.sum(Np_files[:i+1])

    subsample_array[i*N_subsample:(i+1)*N_subsample,
                    0] = np.random.randint(lower_pid, upper_pid,
                                           size=N_subsample)
    subsample_array[i*N_subsample:(i+1)*N_subsample,
                    1:4] = np.random.rand(N_subsample, 3) * 1100.

np.savetxt("test_cats/mpi_test_subsamples.dat", subsample_array)
