# Imports
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
        ll_label[x_index, y_index, z_index] = i

    return ll_list, ll_label


def calc_sep(x1, y1, z1, x2, y2, z2):
    return np.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))


# Get user input parameters
# boxid - ID of the simulation box being used
# halo_type - rockstar or FoF
boxid = "00-0"
# Here for testing, but should be fixed for convenience in the future
halo_type = "rockstar"
sim_name = "AbacusCosmos_1100box_planck"

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

# # Load the halo catalogue
# cat = Halos.make_catalog_from_dir(dirname="{}/{}_{}_{}_halos/"
#                                           "z0.700".format(path_root, sim_name,
#                                                           boxid, halo_type),
#                                   load_subsamples=False, load_pids=False)
# halos = cat.halos

# N_halos = 2
# halos = halos[:N_halos]
# print("Halos to be investigated:\n", halos)

# halo_array = np.zeros((halos.shape[0], 37))
# halo_array[:, 1] = halos["m"]
# halo_array[:, 2] = halos["r"]
# halo_array[:, 3] = halos["vmax"]
# for i in range(N_halos):
#     halo_array[i, 0] = i
#     halo_array[i, 4] = halos[i]["pos"][0]
#     halo_array[i, 5] = halos[i]["pos"][1]
#     halo_array[i, 6] = halos[i]["pos"][2]
#     halo_array[i, 7] = halos[i]["vel"][0]
#     halo_array[i, 8] = halos[i]["vel"][1]
#     halo_array[i, 9] = halos[i]["vel"][2]
# np.savetxt("test_cats/halos.dat", halo_array, fmt=["%i", "%.4e", "%.4f",
#                                                    "%.4f", "%.3f", "%.3f",
#                                                    "%.3f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f", "%.6f", "%.6f",
#                                                    "%.6f"])

halo_array = np.loadtxt("test_cats/halos.dat")
print("Halos to be investigated:\n", halo_array)

# Prepare array to hold the particle velocities
subsamples = Halos.read_uniform_subsample_FoF("{}/{}_{}_FoF_halos/"
                                              "z0.700".format(path_root,
                                                              sim_name,
                                                              boxid))

# window_range = 1.
# subsample_mask = np.zeros(subsamples.shape[0], dtype=bool)
# for i in range(subsamples.shape[0]):
#     subsample_mask[i] = ((abs(subsamples[i]["pos"][0] - halo_array[0, 4]) < window_range) *
#                          (abs(subsamples[i]["pos"][1] - halo_array[0, 5]) < window_range) *
#                          (abs(subsamples[i]["pos"][2] - halo_array[0, 6]) < window_range) *
#                          (abs(subsamples[i]["pos"][0] - halo_array[1, 4]) < window_range) *
#                          (abs(subsamples[i]["pos"][1] - halo_array[1, 5]) < window_range) *
#                          (abs(subsamples[i]["pos"][2] - halo_array[1, 6]) < window_range))
# subsamples[subsample_mask]

# print("Particle subsample to be investigated:\n", subsamples)

pids = subsamples["pid"]
N_pids = pids.size
subsample_array = np.empty((N_pids, 4))
subsample_array[:, 0] = pids
# for i in range(subsamples.shape[0]):
#     subsample_array[i, 1] = subsamples[i]["pos"][0]
#     subsample_array[i, 2] = subsamples[i]["pos"][1]
#     subsample_array[i, 3] = subsamples[i]["pos"][2]

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

print("First 10 particle subsample:\n", subsamples[:10])
print("First 10 IC subsample:\n", subsample_array[:10])
# np.savetxt("test_cats/subsample_array.dat", subsample_array)

xmin = np.min(subsamples["pos"][:, 0])
ymin = np.min(subsamples["pos"][:, 1])
zmin = np.min(subsamples["pos"][:, 2])

Lbox = 1100.
M = 100
box_length = Lbox / M
ll_list, ll_label = linked_list(subsamples["pos"], M, box_length)

# if halo_type == "FoF":
#     halos = halos[halos["subsamp_len"] > 0]
#     halo_array = np.empty((halos.shape[0], 37))
#     halo_array[:, 1] = halos["N"] * 4.e10
#     halo_array[:, 2] = halos["r90"]
#     halo_array[:, 3] = halos["vcirc_max"]

# elif (sim_name == "AbacusCosmos_1100box" and boxid == "planck" and
#       halo_type == "rockstar"):
#     halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
#     halo_array = np.empty((halos.shape[0], 13))
#     halo_array[:, 1] = halos["N"] * 4.e10
#     halo_array[:, 2] = halos["r"]
#     halo_array[:, 3] = halos["vmax"]

# elif halo_type == "rockstar":
#     halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
#     halo_array = np.empty((halos.shape[0], 37))
#     halo_array[:, 1] = halos["m"]
#     halo_array[:, 2] = halos["r"]
#     halo_array[:, 3] = halos["vmax"]

# else:
#     raise ValueError("Unrecognized halo type!")

smoothing_scales = np.logspace(np.log10(0.1), np.log10(60.), 10)
logDsep = ((np.log10(60.) - np.log10(0.1)) / 10.)

for i in range(halo_array.shape[0]):
    # # Abacus halo IDs are too long for the emulator code
    # halo_array[i, 0] = i
    # # Default box limits may be [-550, +550], so shift to always be positive
    # halo_array[i, 4] = halos[i]["pos"][0] - xmin
    # halo_array[i, 5] = halos[i]["pos"][1] - ymin
    # halo_array[i, 6] = halos[i]["pos"][2] - zmin
    # halo_array[i, 7] = halos[i]["vel"][0]
    # halo_array[i, 8] = halos[i]["vel"][1]
    # halo_array[i, 9] = halos[i]["vel"][2]

    smoothed_velocities = np.zeros((len(smoothing_scales)-1, 4))

    i_1 = int(np.floor(((halo_array[i, 4] - xmin) - smoothing_scales[-1]) /
                       box_length))
    i_2 = int(np.floor(((halo_array[i, 4] - xmin) + smoothing_scales[-1]) /
                       box_length))
    j_1 = int(np.floor(((halo_array[i, 5] - ymin) - smoothing_scales[-1]) /
                       box_length))
    j_2 = int(np.floor(((halo_array[i, 5] - ymin) + smoothing_scales[-1]) /
                       box_length))
    k_1 = int(np.floor(((halo_array[i, 6] - zmin) - smoothing_scales[-1]) /
                       box_length))
    k_2 = int(np.floor(((halo_array[i, 6] - zmin) + smoothing_scales[-1]) /
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
                    sep = calc_sep(halo_array[i, 4],
                                   halo_array[i, 5],
                                   halo_array[i, 6],
                                   subsamples[j]["pos"][0] + xoff,
                                   subsamples[j]["pos"][1] + yoff,
                                   subsamples[j]["pos"][2] + zoff)
                    if sep < smoothing_scales[0]:
                        smoothed_velocities[:, 0] += subsample_array[j, 1]
                        smoothed_velocities[:, 1] += subsample_array[j, 2]
                        smoothed_velocities[:, 2] += subsample_array[j, 3]
                        smoothed_velocities[:, 3] += 1

                    elif sep < smoothing_scales[-2]:
                        minbinsep = int(np.floor((np.log10(smoothing_scales[-1]) -
                                                  np.log10(sep)) / logDsep))
                        smoothed_velocities[-minbinsep:, 0] += subsample_array[j, 1]
                        smoothed_velocities[-minbinsep:, 1] += subsample_array[j, 2]
                        smoothed_velocities[-minbinsep:, 2] += subsample_array[j, 3]
                        smoothed_velocities[-minbinsep:, 3] += 1

                    j = ll_list[j]

    for k in range(len(smoothing_scales)-1):
        halo_array[i, 10+3*k] = (smoothed_velocities[k, 0] /
                                 smoothed_velocities[k, 3] * vel_scaling)
        halo_array[i, 11+3*k] = (smoothed_velocities[k, 1] /
                                 smoothed_velocities[k, 3] * vel_scaling)
        halo_array[i, 12+3*k] = (smoothed_velocities[k, 2] /
                                 smoothed_velocities[k, 3] * vel_scaling)

    np.savetxt("test_cats/smoothing_velocities_{}.dat".format(i),
               smoothed_velocities)


# Columns of the resulting catalogue are shown in this header, which has been
# removed for consistency with the existing emulator code
# header = ("#ID M200b R200b Vmax X Y Z VX VY VZ VX_LIN_0.1 VY_LIN_0.1
#            VZ_LIN_0.1 VX_LIN_0.2 VY_LIN_0.2 VZ_LIN_0.2 ...")
np.savetxt("test_cats/halo_lin_vel_smoothed_test.dat",
           halo_array, fmt=["%i", "%.4e", "%.4f", "%.4f", "%.3f", "%.3f",
                            "%.3f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                            "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                            "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                            "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                            "%.6f", "%.6f", "%.6f", "%.6f", "%.6f", "%.6f",
                            "%.6f"])
