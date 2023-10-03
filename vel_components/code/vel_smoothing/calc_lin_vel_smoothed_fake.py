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

    xmin = 0.
    ymin = 0.
    zmin = 0.

    for i in range(N):
        x_index = int(np.floor((positions[i, 0] - xmin) / box_length))
        y_index = int(np.floor((positions[i, 1] - ymin) / box_length))
        z_index = int(np.floor((positions[i, 2] - zmin) / box_length))
        print("x={}, y={}, z={}, x_index={}, y_index={}, z_index="
              "{}".format(positions[i, 0], positions[i, 1], positions[i, 2],
                          x_index, y_index, z_index))
        ll_list[i] = ll_label[x_index, y_index, z_index]
        ll_label[x_index, y_index, z_index] = i + 1

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

vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)

halo_array = np.loadtxt("test_cats/fake_halos.dat")
# print("Halos to be investigated:\n", halo_array)

subsample_array = np.loadtxt("test_cats/fake_subsample.dat")

# print("First 10 IC subsample:\n", subsample_array[:10])

xmin = 0.
ymin = 0.
zmin = 0.

Lbox = 1100.
M = 100
box_length = Lbox / M
ll_list, ll_label = linked_list(subsample_array[:, 4:7], M, box_length)

# print(ll_list)
# print(ll_label)

smoothing_scales = np.logspace(np.log10(0.1), np.log10(60.), 10)
logDsep = ((np.log10(60.) - np.log10(0.1)) / 9.)

smoothed_velocities = np.zeros((len(smoothing_scales)-1, 4))

i_1 = int(np.floor(((halo_array[4] - xmin) - smoothing_scales[-1]) /
                   box_length))
i_2 = int(np.floor(((halo_array[4] - xmin) + smoothing_scales[-1]) /
                   box_length))
j_1 = int(np.floor(((halo_array[5] - ymin) - smoothing_scales[-1]) /
                   box_length))
j_2 = int(np.floor(((halo_array[5] - ymin) + smoothing_scales[-1]) /
                   box_length))
k_1 = int(np.floor(((halo_array[6] - zmin) - smoothing_scales[-1]) /
                   box_length))
k_2 = int(np.floor(((halo_array[6] - zmin) + smoothing_scales[-1]) /
                   box_length))

# print("i_1={}, i_2={}, j_1={}, j_2={}, k_1={}, k_2="
#       "{}".format(i_1, i_2, j_1, j_2, k_1, k_2))

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
            # print("j=", j)
            while j != 0:
                sep = calc_sep(halo_array[4],
                               halo_array[5],
                               halo_array[6],
                               subsample_array[j-1, 4] + xoff,
                               subsample_array[j-1, 5] + yoff,
                               subsample_array[j-1, 6] + zoff)
                print("j-1={}, sep={}".format(j-1, sep))
                if sep < smoothing_scales[0]:
                    print("All bins")
                    smoothed_velocities[:, 0] += subsample_array[j-1, 1]
                    smoothed_velocities[:, 1] += subsample_array[j-1, 2]
                    smoothed_velocities[:, 2] += subsample_array[j-1, 3]
                    smoothed_velocities[:, 3] += 1

                elif sep < smoothing_scales[-2]:
                    minbinsep = int(np.floor((np.log10(smoothing_scales[-1]) -
                                              np.log10(sep)) / logDsep))
                    print("Min bin:", minbinsep)
                    smoothed_velocities[-minbinsep:, 0] += subsample_array[j-1, 1]
                    smoothed_velocities[-minbinsep:, 1] += subsample_array[j-1, 2]
                    smoothed_velocities[-minbinsep:, 2] += subsample_array[j-1, 3]
                    smoothed_velocities[-minbinsep:, 3] += 1

                j = ll_list[j-1]

for k in range(len(smoothing_scales)-1):
    halo_array[10+3*k] = (smoothed_velocities[k, 0] /
                             smoothed_velocities[k, 3] * vel_scaling)
    halo_array[11+3*k] = (smoothed_velocities[k, 1] /
                             smoothed_velocities[k, 3] * vel_scaling)
    halo_array[12+3*k] = (smoothed_velocities[k, 2] /
                             smoothed_velocities[k, 3] * vel_scaling)

correct_x_vel = True
correct_y_vel = True
correct_z_vel = True
correct_N_parts = True
for i in range(len(smoothing_scales)-1):
    if ((smoothed_velocities[i, 0] / smoothed_velocities[i, 3]) !=
        np.mean(subsample_array[:i+1, 1])):
        print("Bin {} has incorrect x velocity: {} when it should be "
              "{}".format(i, smoothed_velocities[i, 0] / smoothed_velocities[i, 3],
                          np.mean(subsample_array[:i+1, 1])))
        correct_x_vel = False

    if ((smoothed_velocities[i, 1] / smoothed_velocities[i, 3]) !=
        np.mean(subsample_array[:i+1, 2])):
        print("Bin {} has incorrect x velocity: {} when it should be "
              "{}".format(i, smoothed_velocities[i, 1] / smoothed_velocities[i, 3],
                          np.mean(subsample_array[:i+1, 2])))
        correct_y_vel = False

    if ((smoothed_velocities[i, 2] / smoothed_velocities[i, 3]) !=
        np.mean(subsample_array[:i+1, 3])):
        print("Bin {} has incorrect x velocity: {} when it should be "
              "{}".format(i, smoothed_velocities[i, 2] / smoothed_velocities[i, 3],
                          np.mean(subsample_array[:i+1, 3])))
        correct_z_vel = False

    if smoothed_velocities[i, 3] != i+1:
        print("Bin {} has incorrect number of parts: "
              "{}".format(i, smoothed_velocities[i, 3]))
        correct_N_parts = False

if correct_N_parts:
    print("All bins have correct number of particles :)")

if correct_x_vel:
    print("All bins have correct x velocity :)")

if correct_y_vel:
    print("All bins have correct y velocity :)")

if correct_z_vel:
    print("All bins have correct z velocity :)")

print(smoothed_velocities)
