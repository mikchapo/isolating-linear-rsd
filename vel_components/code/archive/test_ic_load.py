import numpy as np
import os
from struct import iter_unpack


# tot_ic_size = 0
# # for i in range(375):
# for i in range(1):
# 	tot_ic_size += os.path.getsize("/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/ic_32c/ic_{}".format(i))

# tot_particles = tot_ic_size / 30
# print("Total Particles: {}, 1440^3: {}".format(tot_particles, 1440**3.))

# ic_array = np.empty((int(1440**3/375), 10), dtype=([('i', '<H'), ('j', '<H'), ('k', '<H'), ('x', 'f'), ('y', 'f'), ('z', 'f'), ('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
# ic_array = np.empty((int(1440**3), 3), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))

for i in range(375):
# for i in range(1):
	print("Starting IC file {}".format(i))
	tot_ic_size = os.path.getsize("/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/ic_32c/ic_{}".format(i))
	tot_particles = tot_ic_size / 30
	ic_array = np.empty((int(tot_particles), 3), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
	with open("/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/ic_32c/ic_{}".format(i), "rb") as file:
		bdata = file.read()
		data = iter_unpack("3h6f", bdata)
		for j, line in enumerate(data):
			ic_array[j, :] = line[-3:]

print("Finished")
