# Check the FoF vs Rockstar catalogues to determine the source of the
# v0.1.0, 2022-01-27 - Code started

# Imports
import matplotlib.pyplot as plt
import numpy as np

from AbacusCosmos import Halos

sim_name = "AbacusCosmos_1100box"
boxid = "planck"
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

rockstar_cat = Halos.make_catalog_from_dir(dirname="{}/{}_{}_rockstar_halos/"
                                                   "z0.700".format(path_root,
                                                                   sim_name,
                                                                   boxid),
                                           load_subsamples=False,
                                           load_pids=False)
rockstar_halos = rockstar_cat.halos

rockstar_host_halos = rockstar_halos[rockstar_halos['parent_id'] == -1]["N"]

FoF_cat = Halos.make_catalog_from_dir(dirname="{}/{}_{}_FoF_halos/"
                                              "z0.700".format(path_root,
                                                              sim_name,
                                                              boxid),
                                      load_subsamples=False,
                                      load_pids=False)
FoF_halos = FoF_cat.halos

FoF_mass_hist = np.histogram(FoF_halos["N"] * 4.e10,
                             bins=np.logspace(12., 16., 21))[0]
rockstar_mass_hist = np.histogram(rockstar_halos["N"] * 4.e10,
                                  bins=np.logspace(12., 16., 21))[0]
rockstar_host_mass_hist = np.histogram(rockstar_host_halos * 4.e10,
                                       bins=np.logspace(12., 16., 21))[0]

dlogm = 0.2

bin_centres = 10.**(np.linspace(12., 16., 20) + 0.2)


fig, axes = plt.subplots(2, 1, figsize=(8., 8.), sharex=True,
                         gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
axes[0].plot(bin_centres, FoF_mass_hist, label="Friends-of-friends",
             marker='o')
axes[0].plot(bin_centres, rockstar_mass_hist, label="Rockstar", marker='o')
axes[0].plot(bin_centres, rockstar_host_mass_hist,
             label="Rockstar no subhalos", marker='o')
axes[1].plot(bin_centres[:16], rockstar_mass_hist[:16] / FoF_mass_hist[:16],
             color="C1", marker='o')
axes[1].plot(bin_centres[:16],
             rockstar_host_mass_hist[:16] / FoF_mass_hist[:16],
             color="C2", marker='o')
axes[1].axhline(y=1., color='k', linestyle='-')
axes[1].set_xlabel("Halo Mass [M_sun/h]")
axes[0].set_ylabel("Count")
axes[1].set_xscale("log")
axes[0].set_yscale("log")
axes[0].legend()
plt.savefig("../output/plots/FoF_rockstar_mass_comp.jpg")
