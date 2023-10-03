# Imports
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import os
from struct import unpack_from

from AbacusCosmos import Halos
from lin_vel_funcs import calc_vel_scaling


def load_uniform_subsample(cat_dir, field_len=100., field_height=10.,
                           x0=0., y0=0., z0=0., load_pids=False):
    # Load the halo catalogue
    cat = Halos.make_catalog_from_dir(dirname=cat_dir,
                                      load_uniform_subsample=True,
                                      load_halos=False, load_pids=load_pids)
    subsample = cat.uniform_subsample

    # These filters could be applied in a single step by multiplying them, but
    # since each removes about half of the remaining objects it will require
    # many fewer calculations to apply them one after another
    subsample = subsample[subsample['pos'][:, 0] >= x0]
    subsample = subsample[subsample['pos'][:, 0] < (x0 + field_len)]
    subsample = subsample[subsample['pos'][:, 1] >= y0]
    subsample = subsample[subsample['pos'][:, 1] < (y0 + field_len)]
    subsample = subsample[subsample['pos'][:, 2] >= (z0 - field_height / 2)]
    subsample = subsample[subsample['pos'][:, 2] < (z0 + field_height / 2)]

    return subsample


sim_name = "AbacusCosmos_1100box_planck"
boxid = "00-0"
z_cat = 0.7
smooth_type = "tophat"
R_smooth = 5.
boxsize = 1100.
N_grid = 1100
cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
           "{}_{}_products/{}_{}_FoF_halos/"
           "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
                            z_cat))

'''
# Testing whether the old and new catalogue loading methods agree, if the
# correct smoothed velocity grid is used, and if the smoothed linear velocity
# looks like the total velocity.
# - Both methods give the same nubmer of objects, the same range, and the same
#   position and velocity for the first particle
# - The correct linear velocity grid is used when the offset for position range
#   is included
# - The smoothed linear velocity does not look like the regular velocity, which
# may be a fluke or may mean the velocity grid is off

new_cat = Halos.make_catalog_from_dir(dirname=cat_dir,
                                      load_uniform_subsample=True,
                                      load_halos=False, load_pids=True)
new_cat = new_cat.uniform_subsample

old_cat = Halos.read_uniform_subsample_FoF(cat_dir, boxsize=1100.)

print("Number of objects in new cat:", new_cat["pos"].shape[0])
print("Number of objects in old cat:", old_cat["pos"].shape[0])

print("New cat position range [{:.3f}, {:.3f}"
      "]".format(np.min(new_cat["pos"][:, 0]), np.max(new_cat["pos"][:, 0])))
print("Old cat position range [{:.3f}, {:.3f}"
      "]".format(np.min(old_cat["pos"][:, 0]), np.max(old_cat["pos"][:, 0])))

print("New cat particle 1")
print("Position:", new_cat["pos"][0, :])
print("Velocity:", new_cat["vel"][0, :])

print("Old cat particle 1")
print("Position:", old_cat["pos"][0, :])
print("Velocity:", old_cat["vel"][0, :])

l_grid = boxsize / N_grid
grid_inds = [int((old_cat["pos"][0, 0] + boxsize / 2.) // l_grid),
             int((old_cat["pos"][0, 1] + boxsize / 2.) // l_grid),
             int((old_cat["pos"][0, 2] + boxsize / 2.) // l_grid)]

print("OC-P1 Smoothed vel grid indices:", grid_inds)

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))
vel_grid = np.load("{}/vel_grid_{}_{}_{}_smoothed_v2"
                   ".npy".format(path_root, N_grid, smooth_type, R_smooth))

print("OC-P1 smoothed velocity:", vel_grid[grid_inds[0], grid_inds[1],
                                           grid_inds[2], :])
'''


# Testing whether the unsmoothed velocity grid matches the raw calculation
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))
# subsample = load_uniform_subsample(cat_dir, field_len=1., field_height=1.,
#                                    x0=-501., y0=541., z0=470.5, load_pids=True)

subsample = np.array([([-500.38144, 541.2986, 470.83716],
                       [1640.629, 3.5131242, 379.41745], 101593336)],
                     dtype=[('pos', np.float32, 3), ('vel', np.float32, 3),
                            ('pid', np.uint64)])

print("Cell positions:")
print(subsample["pos"])

print("Cell velocities:")
print(subsample["vel"])

print("Cell pids:")
print(subsample["pid"])

# Convert the subsample object to a numpy array
boxsize = 1100.
pids = subsample["pid"]
N_sample = pids.size
subsample_array = np.empty((N_sample, 13))
subsample_array[:, 0] = pids
# Note: The AbacusCosmos particle loading function arranges the particles
# between [-boxsize, boxsize] by default. Here I've hard coded a shift to
# change the range to [0, boxsize]
subsample_array[:, 1:4] = subsample["pos"] + boxsize / 2.
subsample_array[:, 4:7] = subsample["vel"]
del subsample

# subsample_array = np.empty((1, 13))
# subsample_array[0, 0] = 

# Initialize some variables for loading the initial conditions
# Sort the particle IDs so only one loop of the particles is needed
sort_map = np.argsort(pids)
sorted_pids = pids[sort_map]
last_sort_pid_index = 0
tot_particles = 0
# Boolean to make sure the final value isn't overwritten
first_final = True
ic_dir = ("{}/ic_ngc_nplt_z49.0".format(path_root))
N_ic_files = len(os.listdir(ic_dir))

# Loop through the initial conditions files and find the linear velocity
# for each particle
# for i in range(N_ic_files):
#     print("Starting IC file {}".format(i))
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

#         for j in range(last_sort_pid_index, N_sample):
#             if int(sorted_pids[j] - tot_particles) < N_particles:
#                 if j == (N_sample-1):
#                     print("Final particle, j={}".format(j))
#                     last_sort_pid_index = j

#                     if first_final:
#                         particle_data = unpack_from("3h6f", bdata,
#                                                     offset=int((sorted_pids[j] -
#                                                                 tot_particles)*32))

#                         subsample_array[sort_map[j], 7:10] = particle_data[-3:]
#                         first_final = False

#                 else:
#                     particle_data = unpack_from("3h6f", bdata,
#                                                 offset=int((sorted_pids[j] -
#                                                             tot_particles)*32))

#                     subsample_array[sort_map[j], 7:10] = particle_data[-3:]

#             else:
#                 print("Switched IC at {}".format(j))
#                 last_sort_pid_index = j
#                 break

#     tot_particles += N_particles

j = 0
particle_data = [6.51212575e+002, -2.59398532e+001, 3.94510248e+001]
subsample_array[sort_map[j], 7:10] = particle_data[-3:]

# Scale the initial condition linear velocities to the catalogue redshift
# Load the cosmological parameters of the simulation box
info_path = ("{}/{}_{}_rockstar_halos/"
             "info".format(path_root, sim_name, boxid))
# Calculate the scaling that needs to be applied to the initial condition
# particle velocities
z_ic = 49.0
# vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)
# subsample_array[:, 7:10] *= vel_scaling

print("In array form")
print(subsample_array)

# Load the smoothed velocity grid
vel_grid = np.load("{}/vel_grid_{}_{}_{}_smoothed_v2"
                   ".npy".format(path_root, N_grid, smooth_type, R_smooth))
l_grid = boxsize / N_grid

for i in range(subsample_array.shape[0]):
    subsample_array[i, 10] = vel_grid[int(subsample_array[i, 1] // l_grid),
                                      int(subsample_array[i, 2] // l_grid),
                                      int(subsample_array[i, 3] // l_grid), 0]
    subsample_array[i, 11] = vel_grid[int(subsample_array[i, 1] // l_grid),
                                      int(subsample_array[i, 2] // l_grid),
                                      int(subsample_array[i, 3] // l_grid), 1]
    subsample_array[i, 12] = vel_grid[int(subsample_array[i, 1] // l_grid),
                                      int(subsample_array[i, 2] // l_grid),
                                      int(subsample_array[i, 3] // l_grid), 2]
    grid_inds = [int(subsample_array[i, 1] // l_grid),
                 int(subsample_array[i, 2] // l_grid),
                 int(subsample_array[i, 3] // l_grid)]
    print("Smoothed vel grid indices:", grid_inds)
    print("Smoothed vel grid value:", vel_grid[grid_inds[0], grid_inds[1],
                                               grid_inds[2], :])
    print("Shifted smoothed vel grid value:", vel_grid[grid_inds[0] - 550,
                                                       grid_inds[1] - 550,
                                                       grid_inds[2] - 550, :])

raw_vel_grid = np.load("{}/vel_grid_{}_v2"
                       ".npy".format(path_root, N_grid))
print("Unsmoothed vel grid value:", raw_vel_grid[grid_inds[0], grid_inds[1],
                                                 grid_inds[2], :])
print("Shifted unsmoothed vel grid value:", raw_vel_grid[grid_inds[0] - 550,
                                                         grid_inds[1] - 550,
                                                         grid_inds[2] - 550,
                                                         :])

print("Total number of particles in smoothed velocity grid:",
      np.sum(vel_grid[:, :, :, 3]))
print("Total number of particles in unsmoothed velocity grid:",
      np.sum(raw_vel_grid[:, :, :, 3]))
