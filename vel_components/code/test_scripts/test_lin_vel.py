import numpy as np
import os
from struct import iter_unpack

from AbacusCosmos import Halos


cat = Halos.make_catalog_from_dir(dirname='/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700', load_subsamples=True, load_pids=True)
halos = cat.halos
subsamples = cat.subsamples
pids = subsamples["pid"]

i = 0
print("Subsamples - Start:", halos[i]["subsamp_start"], "\tLen:", halos[i]["subsamp_len"])

subsample_array = np.empty((pids.size, 4))
subsample_array[:, 0] = pids

particle_ids = subsample_array[halos[i]["subsamp_start"]:halos[i]["subsamp_start"]+halos[i]["subsamp_len"], 0]
for pid in particle_ids:
    print(pid)

# sort_map = np.argsort(pids)
# sorted_pids = pids[sort_map]

# last_sort_pid_index = 0
# tot_particles = 0

# for i in range(375):
# # for i in range(2):
#     print("Starting IC file {}".format(i))
#     tot_ic_size = os.path.getsize("/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/{}/ic_{}".format(ic_dir, i))
#     N_particles = tot_ic_size / 32 # Total size of objects is 30, but becomes 32 because of padding
#     ic_array = np.empty(int(N_particles), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
#     with open("/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/{}/ic_{}".format(ic_dir, i), "rb") as file:
#         bdata = file.read()
#         data = iter_unpack("3h6f", bdata)
#         for j, line in enumerate(data):
#             ic_array[j] = line[-3:]

#     for j in range(last_sort_pid_index, pids.size):
#         try:
#             subsample_array[sort_map[j], 1] = ic_array["vx"][int(sorted_pids[j] - tot_particles)]
#             subsample_array[sort_map[j], 2] = ic_array["vy"][int(sorted_pids[j] - tot_particles)]
#             subsample_array[sort_map[j], 3] = ic_array["vz"][int(sorted_pids[j] - tot_particles)]

#         except IndexError:
#             print("Hit IndexError at {}".format(j))
#             last_sort_pid_index = j
#             break

#     tot_particles += N_particles