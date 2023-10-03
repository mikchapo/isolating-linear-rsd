# Check that the difference between old and new halo smoothed lienar velocity
# catalogues are caused by machine precision
# v0.1.0, 2022-09-06 - Code Started

# Imports
import numpy as np

sim_name = "AbacusCosmos_1100box_planck"
boxid = "00-0"

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products/{}_{}_rockstar_halos/"
             "z0.700".format(sim_name, sim_name, boxid, sim_name, boxid))

old_cat = np.loadtxt("{}/halo_lin_vel_1100_tophat_5.0_smoothed_v2"
                     ".dat".format(path_root))
new_cat = np.loadtxt("{}/halo_lin_vel_1100_tophat_5.0_smoothed_v3"
                     ".dat".format(path_root))

vel_diffs = new_cat[:, -3:] - old_cat[:, -3:]
sorted_vel_diffs = np.sort(np.abs(vel_diffs.flatten()))
print("Largest 10 velocity differences:",
      sorted_vel_diffs[-10:])
print("Number of velocity differences greater than 0.1:",
      np.sum(np.abs(vel_diffs.flatten()) > 0.1))
print("Number of velocity differences greater than 1:",
      np.sum(np.abs(vel_diffs.flatten()) > 1))
print("Number of velocity differences greater than 10:",
      np.sum(np.abs(vel_diffs.flatten()) > 10))

edge_halo_mask = (((new_cat[:, 4] % 1) == 0) | ((new_cat[:, 5] % 1) == 0) |
                  ((new_cat[:, 6] % 1) == 0) | ((old_cat[:, 4] % 1) == 0) |
                  ((old_cat[:, 5] % 1) == 0) | ((old_cat[:, 6] % 1) == 0))

print("Number of edge halos:", np.sum(edge_halo_mask))

edge_vel_diffs = new_cat[edge_halo_mask, -3:] - old_cat[edge_halo_mask, -3:]
N_edge_vel_diffs = np.sum(np.abs(edge_vel_diffs.flatten()) > 0)
print("Number of edge velocity differences:", N_edge_vel_diffs)

no_edge_vel_diffs = (new_cat[np.logical_not(edge_halo_mask), -3:] -
                     old_cat[np.logical_not(edge_halo_mask), -3:])
N_no_edge_vel_diffs = np.sum(np.abs(no_edge_vel_diffs.flatten()) > 0)
print("Number of non-edge velocity differences:", N_no_edge_vel_diffs)
sorted_no_edge_vel_diffs = np.sort(np.abs(no_edge_vel_diffs.flatten()))
print("Largest 10 non-edge velocity differences:",
      sorted_no_edge_vel_diffs[-10:])
print("Number of non-edge velocity differences greater than 0.1:",
      np.sum(np.abs(no_edge_vel_diffs.flatten()) > 0.1))
print("Number of non-edge velocity differences greater than 1:",
      np.sum(np.abs(no_edge_vel_diffs.flatten()) > 1))
print("Number of non-edge velocity differences greater than 10:",
      np.sum(np.abs(no_edge_vel_diffs.flatten()) > 10))



'''
new_boundary_x_cat = new_cat[(new_cat[:, 4] % 1) == 0]
new_boundary_y_cat = new_cat[(new_cat[:, 5] % 1) == 0]
new_boundary_z_cat = new_cat[(new_cat[:, 6] % 1) == 0]

print("Number of x new boundary halos:", new_boundary_x_cat.shape[0])
print("Number of y new boundary halos:", new_boundary_y_cat.shape[0])
print("Number of z new boundary halos:", new_boundary_z_cat.shape[0])

print("New boundary halos x:")
print(new_boundary_x_cat[:10, 4:7])
print("New boundary halos y:")
print(new_boundary_y_cat[:10, 4:7])
print("New boundary halos z:")
print(new_boundary_z_cat[:10, 4:7])

new_boundary_halos = np.unique(np.concatenate((np.array(new_boundary_x_cat),
                                               np.array(new_boundary_y_cat),
                                               np.array(new_boundary_z_cat)),
                                              axis=0), axis=0)

print("Number of new boundary halos:", new_boundary_halos.shape)

old_boundary_x_cat = old_cat[(old_cat[:, 4] % 1) == 0]
old_boundary_y_cat = old_cat[(old_cat[:, 5] % 1) == 0]
old_boundary_z_cat = old_cat[(old_cat[:, 6] % 1) == 0]

print("Number of x old boundary halos:", old_boundary_x_cat.shape[0])
print("Number of y old boundary halos:", old_boundary_y_cat.shape[0])
print("Number of z old boundary halos:", old_boundary_z_cat.shape[0])

print("Old boundary halos x:")
print(old_boundary_x_cat[:10, 4:7])
print("Old boundary halos y:")
print(old_boundary_y_cat[:10, 4:7])
print("Old boundary halos z:")
print(old_boundary_z_cat[:10, 4:7])

old_boundary_halos = np.unique(np.concatenate((np.array(old_boundary_x_cat),
                                               np.array(old_boundary_y_cat),
                                               np.array(old_boundary_z_cat)),
                                              axis=0), axis=0)

print("Number of old boundary halos:", old_boundary_halos.shape)
'''


# same = 0
# diff = 0
# for index in boundary_halos:
#     if old_cat[index, -3:] == new_cat[index, -3:]:
#         same += 1
#     else:
#         diff += 1

# print("# Same = {}, # Diff = {}".format(same, diff))


# Change Log
