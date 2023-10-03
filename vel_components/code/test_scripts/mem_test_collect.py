# Imports
import datetime as dt
import gc
import numpy as np
import os

from AbacusCosmos import Halos


start_time = dt.datetime.now()
print("Starting at", start_time)

# Get user input parameters
# boxid - ID of the simulation box being used
# halo_type - rockstar or FoF
boxid = "00-0"
# Here for testing, but should be fixed for convenience in the future
halo_type = "rockstar"
sim_name = "AbacusCosmos_1100box_planck"

# Important global constants, that could be changed to parameters in the future
z_ic = 49.0
ic_type = "ngc_nplt"

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

# Load the halo catalogue
cat = Halos.make_catalog_from_dir(dirname="{}/{}_{}_{}_halos/"
                                          "z0.700".format(path_root, sim_name,
                                                          boxid, halo_type),
                                  load_subsamples=False, load_pids=False)
halos = cat.halos
del cat
gc.collect()

xmin = np.min(halos["pos"][:, 0])
ymin = np.min(halos["pos"][:, 1])
zmin = np.min(halos["pos"][:, 2])

if halo_type == "FoF":
    halos = halos[halos["subsamp_len"] > 0]
    halo_array = np.empty((halos.shape[0], 37))
    halo_array[:, 1] = halos["N"] * 4.e10
    halo_array[:, 2] = halos["r90"]
    halo_array[:, 3] = halos["vcirc_max"]

# elif (sim_name == "AbacusCosmos_1100box" and boxid == "planck" and
#       halo_type == "rockstar"):
#     halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
#     halo_array = np.empty((halos.shape[0], 13))
#     halo_array[:, 1] = halos["N"] * 4.e10
#     halo_array[:, 2] = halos["r"]
#     halo_array[:, 3] = halos["vmax"]

elif halo_type == "rockstar":
    halos = halos[(halos['parent_id'] == -1) * (halos["subsamp_len"] > 0)]
    halo_array = np.empty((halos.shape[0], 37))
    halo_array[:, 1] = halos["m"]
    halo_array[:, 2] = halos["r"]
    halo_array[:, 3] = halos["vmax"]

else:
    raise ValueError("Unrecognized halo type!")

# Abacus halo IDs are too long for the emulator code
halo_array[:, 0] = np.arange(halos.shape[0])
# Default box limits may be [-550, +550], so shift to always be positive
halo_array[:, 4] = halos["pos"][:, 0] - xmin
halo_array[:, 5] = halos["pos"][:, 1] - ymin
halo_array[:, 6] = halos["pos"][:, 2] - zmin
halo_array[:, 7] = halos["vel"][:, 0]
halo_array[:, 8] = halos["vel"][:, 1]
halo_array[:, 9] = halos["vel"][:, 2]

del halos
gc.collect()

# Prepare array to hold the particle velocities
subsamples = Halos.read_uniform_subsample_FoF("{}/{}_{}_FoF_halos/"
                                              "z0.700".format(path_root,
                                                              sim_name,
                                                              boxid))

pids = subsamples["pid"]
N_pids = pids.size
subsample_array = np.empty((N_pids, 7))
subsample_array[:, 0] = pids
subsample_array[:, 1] = subsamples["pos"][:, 0]
subsample_array[:, 2] = subsamples["pos"][:, 1]
subsample_array[:, 3] = subsamples["pos"][:, 2]

del subsamples
gc.collect()
