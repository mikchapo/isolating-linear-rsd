# Check IC offsets follow up. What did you think it meant?

import os
from struct import unpack_from
import sys

boxid = sys.argv[1]
file_id = sys.argv[2]
print("Starting box {}, file {}".format(boxid, file_id))

ic_dir = ("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
          "AbacusCosmos_1100box_{}_products/ic_ngc_nplt_z49.0".format(boxid))

offsets_to_check = [1, 2, 1439, 1440]

tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, file_id))
N_particles = int(tot_ic_size / 32)

with open("{}/ic_{}".format(ic_dir, file_id), "rb") as file:
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

    for i in range(1, 100):
        if (1440**2) * i + 1 < N_particles:
            offsets_to_check.append(int((1440**2) * i - 1))
            offsets_to_check.append(int((1440**2) * i))
            offsets_to_check.append(int((1440**2) * i + 1))

        else:
            break

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

        for j in range(3):
            if particle_data[j] != expected_index[j]:
                print("Oh no, we got a problem!",
                      "File {}, Offset {}, dimension {}".format(file_id,
                                                                offset,
                                                                j),
                      "Expected {}, Got {}".format(expected_index[j],
                                                   particle_data[j]))

        print("Offset {}, particle data:".format(offset), particle_data, "\n")

    particle_data_aend = unpack_from("3h6f", bdata, offset=-64)
    expected_index = [end_index[0], end_index[1], end_index[2] - 1]
    if expected_index[2] < 0:
        expected_index[2] += 1440
        expected_index[1] -= 1
    if expected_index[1] < 0:
        expected_index[1] += 1440
        expected_index[0] -= 1

    for j in range(3):
        if particle_data_aend[j] != expected_index[j]:
            print("Oh no, we got a problem!",
                  "File {}, Offset -2, dimension {}".format(file_id, j),
                  "Expected {}, Got {}".format(expected_index[j],
                                               particle_data_aend[j]))

particle_data_N_index = unpack_from("3h6f", bdata, offset=(N_index - 1) * 32)

for j in range(len(particle_data_N_index)):
    if particle_data_N_index[j] != particle_data_end[j]:
        print("Oh no, we got a problem!",
              "File {}, N index, dimension {}".format(file_id, j),
              "Expected {}, Got {}".format(particle_data_end[j],
                                           particle_data_N_index[j]))
