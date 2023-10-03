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


def linked_list(positions, M, box_length, mins):
    N = positions.shape[0]
    ll_list = np.zeros(N, dtype=int)
    ll_label = np.zeros((M, M, M), dtype=int)

    xmin, ymin, zmin = mins

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

vel_scaling = 1.25

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    full_subsample_array = np.loadtxt("test_cats/fake_subsample.dat")

    print("Particles loaded, elapsed time", dt.datetime.now() - start_time)

    # Sort the particle IDs so only one loop of the particles is needed
    full_subsample_array = full_subsample_array[np.argsort(full_subsample_array[:, 0])]

    N_parts_per_rank = full_subsample_array.shape[0] // size

    lower_bounds = np.empty(size, dtype=int)
    upper_bounds = np.empty(size, dtype=int)
    for i in range(size):
        lower_bounds[i] = i * N_parts_per_rank
        if i != (size - 1):
            upper_bounds[i] = (i + 1) * N_parts_per_rank
        else:
            upper_bounds[i] = full_subsample_array.shape[0]

    sendcounts = (upper_bounds - lower_bounds) * 7

else:
    full_subsample_array = None
    lower_bounds = np.empty(size, dtype=int)
    upper_bounds = None
    sendcounts = np.empty(size, dtype=int)


comm.Barrier()
comm.Bcast(lower_bounds, root=0)
lower_bound = lower_bounds[rank]
comm.Bcast(sendcounts, root=0)
upper_bound = comm.scatter(upper_bounds, root=0)
subsample_array = np.empty((upper_bound - lower_bound, 7))
print("Sendcounts:", sendcounts)
print("Lower bounds:", lower_bounds)
comm.Scatterv([full_subsample_array, sendcounts, lower_bounds*7, MPI.DOUBLE],
              subsample_array, root=0)

del full_subsample_array

ic_array = np.loadtxt("test_cats/fake_ic.dat")

for i in range(lower_bound, upper_bound):
    subsample_array[i - lower_bound, 4:] = ic_array[i, 1:]

del ic_array

print("Rank {}, subsample first 10:".format(rank), subsample_array[:10, :])

halo_array = np.loadtxt("test_cats/fake_halos.dat")

xmin = 0.
ymin = 0.
zmin = 0.

Lbox = 100.
M = 10
box_length = Lbox / M
ll_list, ll_label = linked_list(subsample_array[:, 1:4], M, box_length,
                                (xmin, ymin, zmin))

smoothing_scales = np.logspace(np.log10(0.1), np.log10(60.), 10)
logDsep = ((np.log10(60.) - np.log10(0.1)) / 9.)

smooth_vel_array = np.empty((halo_array.shape[0],
                             4 * (len(smoothing_scales)-1)))

for i in range(halo_array.shape[0]):
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
        smooth_vel_array[i, 4*k] = (smoothed_velocities[k, 0])
        smooth_vel_array[i, 1+4*k] = (smoothed_velocities[k, 1])
        smooth_vel_array[i, 2+4*k] = (smoothed_velocities[k, 2])
        smooth_vel_array[i, 3+4*k] = (smoothed_velocities[k, 3])

    if rank == 0 and ((i + 1) % 100000) == 0:
        print("Completed {} halos, elapsed time".format(i + 1),
              dt.datetime.now() - start_time)

del subsample_array

np.savetxt("test_cats/smooth_vel_array_mpi_{}.txt".format(rank),
           smooth_vel_array)

comm.Barrier()
smooth_vel_arr_list = comm.gather(smooth_vel_array, root=0)

if rank == 0:
    smooth_vel_array = np.sum(smooth_vel_arr_list, axis=0)

    del smooth_vel_arr_list

    for k in range(len(smoothing_scales)-1):
        halo_array[:, 10+3*k] = (smooth_vel_array[:, 4*k] /
                                 smooth_vel_array[i, 3+4*k] * vel_scaling)
        halo_array[:, 11+3*k] = (smooth_vel_array[:, 1+4*k] /
                                 smooth_vel_array[i, 3+4*k] * vel_scaling)
        halo_array[:, 12+3*k] = (smooth_vel_array[:, 2+4*k] /
                                 smooth_vel_array[i, 3+4*k] * vel_scaling)

    del smooth_vel_array

    halo_array[:, 1] -= xmin
    halo_array[:, 2] -= ymin
    halo_array[:, 3] -= zmin

    # Columns of the resulting catalogue are shown in this header, which has
    # been removed for consistency with the existing emulator code
    # header = ("#ID M200b R200b Vmax X Y Z VX VY VZ VX_LIN_0.1 VY_LIN_0.1
    #            VZ_LIN_0.1 VX_LIN_0.2 VY_LIN_0.2 VZ_LIN_0.2 ...")
    np.savetxt("test_cats/halo_lin_vel_smoothed_mpi_sv.dat",
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
