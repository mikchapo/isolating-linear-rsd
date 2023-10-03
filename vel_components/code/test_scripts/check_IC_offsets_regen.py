import os
from struct import unpack_from
import sys

boxid = sys.argv[1]
ic_dir = ("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
          "AbacusCosmos_1100box_{}_products/ic_ngc_nplt_z49.0".format(boxid))

print("Starting box {}".format(boxid))
double_files = 0

tot_particles = 0
N_ic_files = len(os.listdir(ic_dir))
for i in range(N_ic_files):
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    N_particles = int(tot_ic_size / 32)

    with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
        bdata = file.read()
        particle_data_0 = unpack_from("3h6f", bdata, offset=0)
        init_index = [int(particle_data_0[0]), int(particle_data_0[1]),
                      int(particle_data_0[2])]

        particle_data_end = unpack_from("3h6f", bdata, offset=-32)
        end_index = [int(particle_data_end[0]), int(particle_data_end[1]),
                     int(particle_data_end[2])]

    N_index = ((end_index[0] - init_index[0]) * 1440**2 +
               (end_index[1] - init_index[1]) * 1440 +
               (end_index[2] - init_index[2]) + 1)

    tot_particles += N_index

    if (N_particles / N_index) != 1:
        if (N_particles / N_index) != 2:
            print("File {}, N_len: {},".format(i, N_particles),
                  "N_ind: {},".format(N_index),
                  "Ratio: {}".format(N_particles / N_index))

        else:
            double_files += 1

print("Tot: {}, Double Files: {}".format(tot_particles, double_files), "\n")
