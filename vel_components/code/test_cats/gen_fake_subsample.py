import numpy as np

lbox = 100.
N_halos = 1000
N_parts = 10000

halo_array = np.zeros((N_halos, 37))
halo_array[:, 0] = np.arange(N_halos)
halo_array[:, 1] = 55
halo_array[:, 2] = 1000
halo_array[:, 3] = 100
halo_array[:, 4] = np.random.rand(N_halos) * lbox
halo_array[:, 5] = np.random.rand(N_halos) * lbox
halo_array[:, 6] = np.random.rand(N_halos) * lbox
halo_array[:, 7] = 20.
halo_array[:, 8] = 30.
halo_array[:, 9] = 40.

np.savetxt("fake_halos.dat", halo_array)

part_array = np.zeros((N_parts, 7))
part_array[:, 0] = np.arange(N_parts)
part_array[:, 1] = np.random.rand(N_parts) * lbox
part_array[:, 2] = np.random.rand(N_parts) * lbox
part_array[:, 3] = np.random.rand(N_parts) * lbox

np.savetxt("fake_subsample.dat", part_array)

ic_array = np.zeros((N_parts, 4))
ic_array[:, 0] = np.arange(N_parts)
ic_array[:, 1] = (np.random.rand(N_parts) - 0.5) * 200.
ic_array[:, 2] = (np.random.rand(N_parts) - 0.5) * 200.
ic_array[:, 3] = (np.random.rand(N_parts) - 0.5) * 200.

np.savetxt("fake_ic.dat", ic_array)
