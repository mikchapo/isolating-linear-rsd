# Calculate the linear velocity component of halos
# v0.4.0, 2021-12-03 - Added other simulation types and changed FoF back to
#                      M_sun/h mass units

# Imports
import datetime as dt
from mpi4py import MPI
import numpy as np
import os
from struct import unpack_from
import sys

from AbacusCosmos import Halos
from lin_vel_funcs import calc_vel_scaling


def linked_list(positions, M, box_length):
    N = positions.shape[0]
    ll_list = np.zeros(N, dtype=int)
    ll_label = np.zeros((M, M, M), dtype=int)

    xmin = np.min(positions[:, 0])
    ymin = np.min(positions[:, 1])
    zmin = np.min(positions[:, 2])

    for i in range(N):
        x_index = int(np.floor((positions[i, 0] - xmin) / box_length))
        y_index = int(np.floor((positions[i, 1] - ymin) / box_length))
        z_index = int(np.floor((positions[i, 2] - zmin) / box_length))
        ll_list[i] = ll_label[x_index, y_index, z_index]
        ll_label[x_index, y_index, z_index] = i + 1

    return ll_list, ll_label


def calc_sep(x1, y1, z1, x2, y2, z2):
    return np.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))


start_time = dt.datetime.now()
print("Starting at", start_time)

# Get user input parameters
# boxid - ID of the simulation box being used
# halo_type - rockstar or FoF
boxid = sys.argv[1]
# Here for testing, but should be fixed for convenience in the future
halo_type = sys.argv[2]
sim_name = sys.argv[3]

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
                                  load_subsamples=False, load_pids=False)
halos = cat.halos

print("Halos loaded, elapsed time", dt.datetime.now() - start_time)

# Prepare array to hold the particle velocities
subsamples = Halos.read_uniform_subsample_FoF("{}/{}_{}_FoF_halos/"
                                              "z0.700".format(path_root,
                                                              sim_name,
                                                              boxid))

pids = subsamples["pid"]
N_pids = pids.size
subsample_array = np.empty((N_pids, 7))
subsample_array[:, 0] = pids
subsample_array[:, 1] = subsamples["pos"][:, 0]
subsample_array[:, 2] = subsamples["pos"][:, 1]
subsample_array[:, 3] = subsamples["pos"][:, 2]

del subsamples

print("Particles loaded, elapsed time", dt.datetime.now() - start_time)

# Sort the particle IDs so only one loop of the particles is needed
subsample_array = subsample_array[np.argsort(pids)]
sorted_pids = subsample_array[:, 0]

ic_dir = ("{}/ic_{}_z{}".format(path_root, ic_type, z_ic))
N_ic_files = len(os.listdir(ic_dir))

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
    if rank == 0:
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
                # print("Switched IC at {}".format(j))
                last_sort_pid_index = j
                break

    tot_particles += N_particles

comm.Barrier()
lin_vel_arr_list = comm.gather(lin_vel_array, root=0)

if rank == 0:
    subsample_array[4:] = np.concatenate(lin_vel_arr_list, axis=0)

comm.Barrier()
comm.Bcast(subsample_array, root=0)

xmin = np.min(subsample_array[:, 1])
ymin = np.min(subsample_array[:, 2])
zmin = np.min(subsample_array[:, 3])

Lbox = 1100.
M = 100
box_length = Lbox / M
ll_list, ll_label = linked_list(subsample_array[:, 1:4], M, box_length)

if halo_type == "FoF":
    halos = halos[halos["subsamp_len"] > 0]
    halo_array = np.empty((halos.shape[0], 37))
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
    halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
    halo_array = np.empty((halos.shape[0], 37))
    halo_array[:, 1] = halos["m"]
    halo_array[:, 2] = halos["r"]
    halo_array[:, 3] = halos["vmax"]

else:
    raise ValueError("Unrecognized halo type!")

# Abacus halo IDs are too long for the emulator code
halo_array[:, 0] = np.arange(halos.shape[0])
# Default box limits may be [-550, +550], so shift to always be positive
halo_array[:, 4] = halos["pos"][:, 0] - xmin
halo_array[:, 5] = halos["pos"][:, 1] - ymin
halo_array[:, 6] = halos["pos"][:, 2] - zmin
halo_array[:, 7] = halos["vel"][:, 0]
halo_array[:, 8] = halos["vel"][:, 1]
halo_array[:, 9] = halos["vel"][:, 2]

if rank == 0:
    print("Halo array created, elapsed time", dt.datetime.now() - start_time)

smoothing_scales = np.logspace(np.log10(0.1), np.log10(60.), 10)
logDsep = ((np.log10(60.) - np.log10(0.1)) / 9.)

if rank == 0:
    print("Number of halos:", halos.shape[0])

N_halos_per_rank = halos.shape[0] // size
lower_bound = rank * N_halos_per_rank
upper_bound = (rank + 1) * N_halos_per_rank
if rank == (size - 1):
    upper_bound = halos.shape[0]

N_subsample = upper_bound - lower_bound
smooth_vel_array = np.empty((N_subsample, 3 * len(smoothing_scales)-1))

for i in range(lower_bound, upper_bound):
    smoothed_velocities = np.zeros((len(smoothing_scales)-1, 4))

    i_1 = int(np.floor(((halos[i]["pos"][0] - xmin) - smoothing_scales[-1]) /
                       box_length))
    i_2 = int(np.floor(((halos[i]["pos"][0] - xmin) + smoothing_scales[-1]) /
                       box_length))
    j_1 = int(np.floor(((halos[i]["pos"][1] - ymin) - smoothing_scales[-1]) /
                       box_length))
    j_2 = int(np.floor(((halos[i]["pos"][1] - ymin) + smoothing_scales[-1]) /
                       box_length))
    k_1 = int(np.floor(((halos[i]["pos"][2] - zmin) - smoothing_scales[-1]) /
                       box_length))
    k_2 = int(np.floor(((halos[i]["pos"][2] - zmin) + smoothing_scales[-1]) /
                       box_length))

    for il in range(i_1, i_2+1):
        if il < 0:
            im = M + il
            xoff = -Lbox
        elif il > (M - 1):
            im = il - M
            xoff = Lbox
        else:
            im = il
            xoff = 0.

        for jl in range(j_1, j_2+1):
            if jl < 0:
                jm = M + jl
                yoff = -Lbox
            elif jl > (M - 1):
                jm = jl - M
                yoff = Lbox
            else:
                jm = jl
                yoff = 0.

            for kl in range(k_1, k_2+1):
                if kl < 0:
                    km = M + kl
                    zoff = -Lbox
                elif kl > (M - 1):
                    km = kl - M
                    zoff = Lbox
                else:
                    km = kl
                    zoff = 0.

                j = ll_label[im, jm, km]
                while j != 0:
                    sep = calc_sep(halos[i]["pos"][0],
                                   halos[i]["pos"][1],
                                   halos[i]["pos"][2],
                                   subsample_array[j-1, 1] + xoff,
                                   subsample_array[j-1, 2] + yoff,
                                   subsample_array[j-1, 3] + zoff)
                    if sep < smoothing_scales[0]:
                        smoothed_velocities[:, 0] += subsample_array[j-1, 4]
                        smoothed_velocities[:, 1] += subsample_array[j-1, 5]
                        smoothed_velocities[:, 2] += subsample_array[j-1, 6]
                        smoothed_velocities[:, 3] += 1

                    elif sep < smoothing_scales[-2]:
                        minbinsep = int(np.floor((np.log10(smoothing_scales[-1]) -
                                                  np.log10(sep)) / logDsep))
                        smoothed_velocities[-minbinsep:, 0] += subsample_array[j-1, 4]
                        smoothed_velocities[-minbinsep:, 1] += subsample_array[j-1, 5]
                        smoothed_velocities[-minbinsep:, 2] += subsample_array[j-1, 6]
                        smoothed_velocities[-minbinsep:, 3] += 1

                    j = ll_list[j-1]

    for k in range(len(smoothing_scales)-1):
        smooth_vel_array[i - lower_bound, 0+3*k] = (smoothed_velocities[k, 0] /
                                                    smoothed_velocities[k, 3] *
                                                    vel_scaling)
        smooth_vel_array[i - lower_bound, 1+3*k] = (smoothed_velocities[k, 1] /
                                                    smoothed_velocities[k, 3] *
                                                    vel_scaling)
        smooth_vel_array[i - lower_bound, 2+3*k] = (smoothed_velocities[k, 2] /
                                                    smoothed_velocities[k, 3] *
                                                    vel_scaling)

    if rank == 0 and ((i + 1) % 100000) == 0:
        print("Completed {} halos, elapsed time".format(i + 1),
              dt.datetime.now() - start_time)

comm.Barrier()
smooth_vel_arr_list = comm.gather(smooth_vel_array, root=0)

if rank == 0:
    halo_array[10:] = np.concatenate(smooth_vel_arr_list, axis=0)

    # Columns of the resulting catalogue are shown in this header, which has
    # been removed for consistency with the existing emulator code
    # header = ("#ID M200b R200b Vmax X Y Z VX VY VZ VX_LIN_0.1 VY_LIN_0.1
    #            VZ_LIN_0.1 VX_LIN_0.2 VY_LIN_0.2 VZ_LIN_0.2 ...")
    np.savetxt("{}/{}_{}_{}_halos/z0.700/"
               "halo_lin_vel_smoothed.dat".format(path_root, sim_name, boxid,
                                                  halo_type),
               halo_array, fmt=["%i", "%.4e", "%.4f", "%.4f", "%.3f", "%.3f",
                                "%.3f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                                "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                                "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                                "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                                "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                                "%.6f"])

    print("Catalogue saved, and code complete! Total elapsed time",
          dt.datetime.now() - start_time)

# Change Log
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
