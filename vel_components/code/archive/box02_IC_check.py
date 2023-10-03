# Check box 02 for IC error
# v0.1.0, 2021-11-19 - Copied from calc_lin_vel_IC_test.py

# Imports
import numpy as np
import os
from struct import unpack_from

from AbacusCosmos import Halos


# Run parameters
boxid = "02"
halo_type = "rockstar"
z_cat = 0.7
z_ic = 49.0
ic_type = "ngc_nplt"

# Load the cosmological parameters of the simulation box
info_path = ("/home/mj3chapm/scratch/abacus/"
             "AbacusCosmos_1100box_products/"
             "AbacusCosmos_1100box_{}_products/"
             "AbacusCosmos_1100box_{}_rockstar_halos/"
             "info".format(boxid, boxid))

# Load the cosmological parameters of the simulation box
cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_{}_products/"
                          "AbacusCosmos_1100box_{}_rockstar_halos/"
                          "info/cosmo_params.dat".format(boxid, boxid))

# Load the halo catalogue
cat = Halos.make_catalog_from_dir(dirname="/home/mj3chapm/scratch/abacus/"
                                          "AbacusCosmos_1100box_products/"
                                          "AbacusCosmos_1100box_{}_products/"
                                          "IC_test/".format(boxid),
                                  load_subsamples=True, load_pids=True,
                                  halo_type="Rockstar")
halos = cat.halos

# Prepare array to hold the particle velocities
subsamples = cat.subsamples
pids = np.array([part[0] for part in subsamples])
print("pids size:", pids.size)
print("Min particle ID:", np.min(pids))
print("Max particle ID:", np.max(pids))

# # Sort the particle IDs so only one loop of the particles is needed
# sort_map = np.argsort(pids)
# sorted_pids = pids[sort_map]

# Change Log
