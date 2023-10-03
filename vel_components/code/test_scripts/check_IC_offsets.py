import os
from struct import unpack_from
import sys

boxid = sys.argv[1]
ic_dir = ("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
          "AbacusCosmos_1100box_{}_products/ic_ngc_nplt_z49.0".format(boxid))

print("Starting box {}".format(boxid))
double_files = 0

offsets_to_check = [1, 2, 1439, 1440, (1440**2 - 1), (1440**2)]

tot_particles = 0
N_ic_files = len(os.listdir(ic_dir))
expected_init_index = [0, 0, 0]
for i in range(N_ic_files):
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    N_particles = int(tot_ic_size / 32)

    with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
        bdata = file.read()
        particle_data_0 = unpack_from("3h6f", bdata, offset=0)
        init_index = [int(particle_data_0[0]), int(particle_data_0[1]),
                      int(particle_data_0[2])]
        # print("Expected Init Index:", expected_init_index)
        # print("Init Index:", init_index)
        for j in range(3):
            if init_index[j] != expected_init_index[j]:
                print("Oh no, we got a problem!",
                      "File {}, Offset 0, dimension {}".format(i, j),
                      "Expected {}, Got {}".format(expected_init_index[j],
                                                   init_index[j]))
        for offset in offsets_to_check:
            particle_data = unpack_from("3h6f", bdata, offset=offset*32)
            expected_index = [init_index[0] + (offset // 1440**2), 0,
                              0]
            expected_index[2] = init_index[2] + (offset % 1440)
            if expected_index[2] >= 1440:
                expected_index[2] -= 1440
                expected_index[1] += 1
            expected_index[1] += init_index[1] + ((offset % 1440**2) // 1440)
            if expected_index[1] >= 1440:
                expected_index[1] -= 1440
                expected_index[0] += 1
            # print("Offset:", offset)
            # print("Expected Index:", expected_index)
            # print("Particle data:", particle_data)
            for j in range(3):
                if particle_data[j] != expected_index[j]:
                    print("Oh no, we got a problem!",
                          "File {}, Offset {}, dimension {}".format(i, offset,
                                                                    j),
                          "Expected {}, Got {}".format(expected_index[j],
                                                       particle_data[j]))

        particle_data_end = unpack_from("3h6f", bdata, offset=-32)
        end_index = [int(particle_data_end[0]), int(particle_data_end[1]),
                     int(particle_data_end[2])]
        # print("End Index:", end_index)

        particle_data_aend = unpack_from("3h6f", bdata, offset=-64)
        expected_index = [end_index[0], end_index[1], end_index[2] - 1]
        if expected_index[2] < 0:
            expected_index[2] += 1440
            expected_index[1] -= 1
        if expected_index[1] < 0:
            expected_index[1] += 1440
            expected_index[0] -= 1
        # print("Expected A. End index:", expected_index)
        # print("Prticle Data A. End:", particle_data_aend)
        for j in range(3):
            if particle_data_aend[j] != expected_index[j]:
                print("Oh no, we got a problem!",
                      "File {}, Offset -2, dimension {}".format(i, j),
                      "Expected {}, Got {}".format(expected_index[j],
                                                   particle_data_aend[j]))

        expected_init_index = [end_index[0], end_index[1], end_index[2] + 1]
        if expected_init_index[2] >= 1440:
            expected_init_index[2] -= 1440
            expected_init_index[1] += 1
        if expected_init_index[1] >= 1440:
            expected_init_index[1] -= 1440
            expected_init_index[0] += 1

        # print("Next Expectected init_index:", expected_init_index)

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
