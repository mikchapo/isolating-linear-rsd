import matplotlib.pyplot as plt
import numpy as np


# ic_locations = np.empty((1440**3, 3))
# locs = np.arange(1440)
# ic_locations[:, 0] = np.repeat(locs, 1440**2)
# print("Finished x")
# ic_locations[:, 1] = np.tile(np.repeat(locs, 1440), 1440)
# print("Finished y")
# ic_locations[:, 2] = np.tile(locs, 1440**2)
# print("Finished z")

# print("First 10")
# print(ic_locations[:10, :])

# print("First transition")
# print(ic_locations[1435:1445, :])

# print("Second transition")
# print(ic_locations[1440**2-5:1440**2+5, :])

# print("End 10")
# print(ic_locations[-10:, :])

# for i in range(1440):
#     for j in range(1440):
#     	for k in range(1440):
#     		ic_locations[i*(1440**2) + j*1440 + k, :] = [i, j, k]

#     if ((i + 1) % 10) == 0:
#     	print("Finished ", (i+1))

pids = np.loadtxt("../../output/data_products/halo_0_pids.dat")
xvals = np.empty(pids.size)
yvals = np.empty(pids.size)
zvals = np.empty(pids.size)

for i in range(pids.size):
	xvals[i] = pids[i] / 1440**2
	yvals[i] = (pids[i] % 1440**2) / 1440
	zvals[i] = pids[i] % 1440

plt.figure(0, figsize=(8., 6.))
plt.hist(xvals)
plt.xlabel("x [Mpc/h]")
plt.title("Halo 0 subsample x values")
plt.savefig("../../output/plots/halo_0_xvals.jpg")

plt.figure(1, figsize=(8., 6.))
plt.hist(yvals)
plt.xlabel("y [Mpc/h]")
plt.title("Halo 0 subsample y values")
plt.savefig("../../output/plots/halo_0_yvals.jpg")

plt.figure(2, figsize=(8., 6.))
plt.hist(zvals)
plt.xlabel("z [Mpc/h]")
plt.title("Halo 0 subsample z values")
plt.savefig("../../output/plots/halo_0_zvals.jpg")