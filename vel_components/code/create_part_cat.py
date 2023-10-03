# Downsample the particle catalogue for calculating the pairwise velocity
# v0.2.0, 2021-12-10 - Updated with final velocity scaling and for gener

# Imports
import camb
import numpy as np
import os
from struct import iter_unpack, unpack_from
import sys

from lin_vel_funcs import calc_vel_scaling


boxid = sys.argv[1]
sample_size = int(sys.argv[2])

z_ic = 49.
z_cat = 0.7
info_path = ("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
             "AbacusCosmos_1100box_{}_products/"
             "AbacusCosmos_1100box_{}_rockstar_halos/"
             "info".format(boxid, boxid))
ic_dir = ("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
          "AbacusCosmos_1100box_{}_products/ic_ngc_nplt_z49.0".format(boxid))

vel_scaling, disp_scaling = calc_vel_scaling(z_ic, z_cat, info_path,
                                             calc_disp_scaling=True)

seed = 1337
rng = np.random.default_rng(seed)
sample_indices = np.empty(0)
while sample_indices.size != sample_size:
    random_indices = rng.integers(1440**3, size=(sample_size -
                                                 sample_indices.size))
    sample_indices = np.unique(np.concatenate((sample_indices,
                                               random_indices)))
    print("Try {}, Seed is {}, # of Indices is {}".format(seed - 1337 + 1,
                                                          seed,
                                                          sample_indices.size))
    seed += 1

last_sort_index = 0
tot_particles = 0

particle_sample = np.empty((sample_size, 10))
particle_index = 0

for i in range(375):
    print("Starting IC file {}".format(i))
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    # Total size of objects is 30, but becomes 32 because of padding
    N_particles = tot_ic_size / 32
    with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
        bdata = file.read()

        for j in range(last_sort_index, sample_size):
            if (sample_indices[j] - tot_particles) < N_particles:
                particle_data = unpack_from("3h6f", bdata,
                                            offset=int((sample_indices[j] -
                                                        tot_particles)*32))

                particle_sample[particle_index, 0] = sample_indices[j]
                particle_sample[particle_index, 1:4] = particle_data[:3]
                for k in range(3):
                    particle_sample[particle_index, 4+k] = ((particle_data[k] /
                                                            1440 * 1050.) +
                                                            (particle_data[3+k]
                                                             * disp_scaling))
                    particle_sample[particle_index, 7+k] = (particle_data[6+k]
                                                            * vel_scaling)
                particle_index += 1

            else:
                print("Switched IC at {}".format(j))
                last_sort_index = j
                break

    tot_particles += N_particles

header = "#ID\tI\tJ\tK\tX\tY\tZ\tVX\tVY\tVZ"
np.savetxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
           "AbacusCosmos_1100box_{}_products/"
           "AbacusCosmos_1100box_{}_rockstar_halos/z0.700/"
           "ic_part_{}.dat".format(boxid, boxid, sample_size),
           particle_sample, header=header)


# Change Log
# v0.1.1, 2021-10-25 - Updated for PEP8 linting
# v0.1.0 - v0.1.1 - Many bug fixes and changes
# v0.1.0, 2021-09-13 - Code started with snippets from calc_lin_vel.py
