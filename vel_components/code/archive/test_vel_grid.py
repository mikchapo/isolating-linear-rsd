# Put particle linear velocities on a grid
# v0.1.1, 2022-06-02 - Added minimum position check

# Imports
import datetime as dt
import numpy as np
import os
from struct import unpack_from
import sys

from AbacusCosmos import Halos
from lin_vel_funcs import calc_vel_scaling


'''
start_time = dt.datetime.now()
print("Starting at", start_time)

# Get user input parameters
boxid = sys.argv[1]
sim_name = sys.argv[2]
N_grid = int(sys.argv[3])

# Important global constants, that could be changed to parameters in the future
z_cat = 0.7
z_ic = 49.0
ic_type = "ngc_nplt"
boxsize = 1100.

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

# Load the cosmological parameters of the simulation box
info_path = ("{}/{}_{}_rockstar_halos/"
             "info".format(path_root, sim_name, boxid))

# Calculate the scaling that needs to be applied to the initial condition
# particle velocities
vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)

# Prepare array to hold the particle velocities
# subsample = Halos.read_uniform_subsample_FoF("{}/{}_{}_FoF_halos/"
#                                               "z0.700".format(path_root,
#                                                               sim_name,
#                                                               boxid),
#                                               boxsize=boxsize)

subsample = np.array([([-500.38144, 541.2986, 470.83716],
                       [1640.629, 3.5131242, 379.41745], 101593336)],
                     dtype=[('pos', np.float32, 3), ('vel', np.float32, 3),
                            ('pid', np.uint64)])

l_grid = boxsize / N_grid
vel_grid = np.zeros((N_grid, N_grid, N_grid, 4))

print("Particles loaded, elapsed time", dt.datetime.now() - start_time)

# Sort the particle IDs so only one loop of the particles is needed
subsample = subsample[np.argsort(subsample["pid"])]
N_subsample = subsample["pid"].size
last_index = 0
tot_particles = 0

# xmin = np.floor(np.min(subsample["pos"][:, 0]))
# ymin = np.floor(np.min(subsample["pos"][:, 1]))
# zmin = np.floor(np.min(subsample["pos"][:, 2]))

xmin = -550
ymin = -550
zmin = -550

print("Minimum positions are:", xmin, ymin, zmin)


ic_dir = ("{}/ic_{}_z{}".format(path_root, ic_type, z_ic))

particle_data = [6.51212575e+002, -2.59398532e+001, 3.94510248e+001]

j = 0
vel_grid[int((subsample["pos"][j][0] - xmin) // l_grid),
         int((subsample["pos"][j][1] - ymin) // l_grid),
         int((subsample["pos"][j][2] - zmin) // l_grid), :-1] += particle_data[-3:]
vel_grid[int((subsample["pos"][j][0] - xmin) // l_grid),
         int((subsample["pos"][j][1] - ymin) // l_grid),
         int((subsample["pos"][j][2] - zmin) // l_grid), -1] += 1.
# N_ic_files = len(os.listdir(ic_dir))
# print("Number of IC Files: {}".format(N_ic_files))
# for i in range(N_ic_files):
#     print("Starting IC file {}, elapsed time".format(i),
#           dt.datetime.now() - start_time)
#     with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
#         bdata = file.read()

#         particle_data_first = unpack_from("3h6f", bdata, offset=0)
#         particle_data_last = unpack_from("3h6f", bdata, offset=-32)
#         N_particles = ((particle_data_last[0] -
#                         particle_data_first[0]) * 1440**2 +
#                        (particle_data_last[1] -
#                         particle_data_first[1]) * 1440 +
#                        (particle_data_last[2] -
#                         particle_data_first[2]) + 1)

#         for j in range(last_index, N_subsample):
#             if int(subsample["pid"][j] - tot_particles) < N_particles:
#                 particle_data = unpack_from("3h6f", bdata,
#                                             offset=int((subsample["pid"][j] -
#                                                         tot_particles)*32))

#                 vel_grid[int((subsample["pos"][j][0] - xmin) // l_grid),
#                          int((subsample["pos"][j][1] - ymin) // l_grid),
#                          int((subsample["pos"][j][2] - zmin) // l_grid), :-1] += particle_data[-3:]
#                 vel_grid[int((subsample["pos"][j][0] - xmin) // l_grid),
#                          int((subsample["pos"][j][1] - ymin) // l_grid),
#                          int((subsample["pos"][j][2] - zmin) // l_grid), -1] += 1.

#             else:
#                 print("Switched IC at {}".format(j))
#                 last_index = j
#                 break

#     tot_particles += N_particles

print("Finished populating grid, elapsed time", dt.datetime.now() - start_time)

# vel_grid[:, :, :, 0] *= vel_scaling / vel_grid[:, :, :, -1]
# vel_grid[:, :, :, 1] *= vel_scaling / vel_grid[:, :, :, -1]
# vel_grid[:, :, :, 2] *= vel_scaling / vel_grid[:, :, :, -1]

grid_inds = [49, 1091, 1020]
print("Grid values at expected index:")
print(vel_grid[grid_inds[0], grid_inds[1], grid_inds[2], :])
# np.save("{}/vel_grid_{}_v2".format(path_root, N_grid), vel_grid)

print("Grid saved, elapsed time", dt.datetime.now() - start_time)
'''


'''
xmin = -550.
ymin = -550.
zmin = -550.
positions = [-475.5, 4.6, 253.8]
velocities = [50., -25., 75]

boxsize = 1100.
N_grid = int(1100)
l_grid = boxsize / N_grid
print("l_grid={}, expected 1.".format(l_grid))

vel_grid = np.zeros((N_grid, N_grid, N_grid, 4))
vel_grid[int((positions[0] - xmin) // l_grid),
         int((positions[1] - ymin) // l_grid),
         int((positions[2] - zmin) // l_grid), :-1] += velocities
vel_grid[int((positions[0] - xmin) // l_grid),
         int((positions[1] - ymin) // l_grid),
         int((positions[2] - zmin) // l_grid), -1] += 1.

grid_inds = [int((positions[0] - xmin) // l_grid),
             int((positions[1] - ymin) // l_grid),
             int((positions[2] - zmin) // l_grid)]

print("Expected grid indices:", [74, 554, 803])
print("Resulting grid indices:", grid_inds)
print("vel_grid:", vel_grid[grid_inds[0], grid_inds[1], grid_inds[2]])
'''

boxid = "00-0"
sim_name = "AbacusCosmos_1100box_planck"
N_grid = 1100

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

# TEST CHANGE
vel_grid = np.load("{}/vel_grid_{}_v2_test.npy".format(path_root, N_grid))

# TEST CHANGE
vel_grid_orig = np.load("{}/vel_grid_{}_v2.npy".format(path_root, N_grid))

grid_inds = [49, 1091, 1020]

# TEST CHANGE
print("Current velocty grid test cell")
print(vel_grid[grid_inds[0], grid_inds[1], grid_inds[2], :])
print("Old velocty grid test cell")
print(vel_grid_orig[grid_inds[0], grid_inds[1], grid_inds[2], :])

print("Shifted current velocty grid test cell")
print(vel_grid[grid_inds[0] - 550, grid_inds[1] - 550, grid_inds[2] - 550, :])
print("Shifted old velocty grid test cell")
print(vel_grid_orig[grid_inds[0] - 550, grid_inds[1] - 550, grid_inds[2] - 550,
      :])


print("Grid saved, elapsed time", dt.datetime.now() - start_time)