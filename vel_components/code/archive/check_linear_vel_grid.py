import numpy as np

# vg550 = np.load("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_00_products/vel_grid_550.npy")
vg1100 = np.load("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_00_products/vel_grid_1100.npy")
vg1100_v2 = np.load("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_00_products/vel_grid_1100_v2.npy")

# print("Shape 550 is", vg550.shape)
print("Shape 1100 is", vg1100.shape)
print("Shape 1100 v2 is", vg1100_v2.shape)

# print("Tot. Part. 550 is {}".format(np.sum(vg550[:, :, :, 3])))
print("Tot. Part. 1100 is {}".format(np.sum(vg1100[:, :, :, 3])))
print("Tot. Part. 1100 v2 is {}".format(np.sum(vg1100_v2[:, :, :, 3])))

print("Excerpt from grid")
#print(vg550[:10, 0, 0, :])
print(vg1100[:10, 0, 0, :])
print(vg1100_v2[:10, 0, 0, :])
