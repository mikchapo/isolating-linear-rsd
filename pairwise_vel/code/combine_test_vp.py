# Combine testbox pairwise velocities
# v0.1.0, 2022-08-15 - Code started, copied from plot_vp_check_halo_type.py

# Imports
import numpy as np

simname = "AbacusCosmos_1100box_planck"
boxids = ["00-0", "00-1", "00-2", "00-3", "00-4",
          "00-5", "00-6", "00-7", "00-8", "00-9",
          "00-10", "00-11", "00-12", "00-13", "00-14",
          "00-15", "00-16", "00-17", "00-18", "00-19"]
# boxids = ["00-0", "00-1", "00-3", "00-4",
#           "00-5", "00-8",
#           "00-11", "00-12", "00-13", "00-14",
#           "00-15", "00-17", "00-18", "00-19"]
mass_bin_lims = np.arange(12., 15.5, 0.5)
Nseps = 80
seps = np.logspace(-2, 2, Nseps, endpoint=False)
vel_components = ["v-tot", "v-lin", "v-nl"]

combined_vp = np.zeros((5, mass_bin_lims.size-1, 3, Nseps))

for boxid in boxids:
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
                 "{}_{}_rockstar_halos".format(simname, simname, boxid,
                                               simname, boxid))

    for i in range(mass_bin_lims.size-1):
        for j, component in enumerate(vel_components):
            vp_path = ("{}/z0.700/pairwise_vel/"
                       "halo_{}_{:.2f}-{:.2f}_vp_v2_all-halos"
                       ".dat".format(path_root, component,
                                     mass_bin_lims[i],
                                     mass_bin_lims[i+1]))
            vps = np.loadtxt(vp_path)
            vps[:, 1] = np.nan_to_num(vps[:, 1]) * vps[:, 3]
            for k in range(vps.shape[0]):
                if vps[k, 3] == 1.:
                    vps[k, 2] = vps[k, 1]**2.
                else:
                    vps[k, 2] = np.nan_to_num(vps[k, 2] * (vps[k, 3] - 1) +
                                              (vps[k, 1] * vps[k, 1]) /
                                              vps[k, 3])

            combined_vp[j, i, 0, :] += vps[:, 1]
            combined_vp[j, i, 1, :] += vps[:, 2]
            combined_vp[j, i, 2, :] += vps[:, 3]

        for j, component in enumerate(vel_components[1:]):
            vp_path = ("{}/z0.700/pairwise_vel/"
                       "halo_{}_{:.2f}-{:.2f}_vp_v2_vel-sm_all-halos"
                       ".dat".format(path_root, component,
                                     mass_bin_lims[i],
                                     mass_bin_lims[i+1]))
            vps = np.loadtxt(vp_path)
            vps[:, 1] = np.nan_to_num(vps[:, 1]) * vps[:, 3]
            for k in range(vps.shape[0]):
                if vps[k, 3] == 1.:
                    vps[k, 2] = vps[k, 1]**2.
                else:
                    vps[k, 2] = np.nan_to_num(vps[k, 2] * (vps[k, 3] - 1) +
                                              (vps[k, 1] * vps[k, 1]) /
                                              vps[k, 3])
            combined_vp[3+j, i, 0, :] += vps[:, 1]
            combined_vp[3+j, i, 1, :] += vps[:, 2]
            combined_vp[3+j, i, 2, :] += vps[:, 3]

comb_path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                  "combined_pairwise_vel".format(simname))

for i in range(mass_bin_lims.size-1):
    for j, component in enumerate(vel_components):
        combined_vp[j, i, 1, :] = ((combined_vp[j, i, 1, :] -
                                   combined_vp[j, i, 0, :]**2. /
                                   combined_vp[j, i, 2, :]) /
                                   (combined_vp[j, i, 2, :] - 1))
        combined_vp[j, i, 0, :] = (combined_vp[j, i, 0, :] /
                                   combined_vp[j, i, 2, :])

        vp_path = ("{}/halo_{}_{:.2f}-{:.2f}_vp_v2_all-halos"
                   ".dat".format(comb_path_root, component,
                                 mass_bin_lims[i],
                                 mass_bin_lims[i+1]))
        vp_array = np.empty((Nseps, 4))
        vp_array[:, 0] = seps
        vp_array[:, 1] = combined_vp[j, i, 0, :]
        vp_array[:, 2] = combined_vp[j, i, 1, :]
        vp_array[:, 3] = combined_vp[j, i, 2, :]
        np.savetxt(vp_path, vp_array)

    for j, component in enumerate(vel_components[1:]):
        combined_vp[3+j, i, 1, :] = ((combined_vp[3+j, i, 1, :] -
                                      combined_vp[3+j, i, 0, :]**2. /
                                      combined_vp[3+j, i, 2, :]) /
                                     (combined_vp[3+j, i, 2, :] - 1))
        combined_vp[3+j, i, 0, :] = (combined_vp[3+j, i, 0, :] /
                                     combined_vp[3+j, i, 2, :])

        vp_path = ("{}/halo_{}_{:.2f}-{:.2f}_vp_v2_vel-sm_all-halos"
                   ".dat".format(comb_path_root, component,
                                 mass_bin_lims[i],
                                 mass_bin_lims[i+1]))
        vp_array = np.empty((Nseps, 4))
        vp_array[:, 0] = seps
        vp_array[:, 1] = combined_vp[3+j, i, 0, :]
        vp_array[:, 2] = combined_vp[3+j, i, 1, :]
        vp_array[:, 3] = combined_vp[3+j, i, 2, :]
        np.savetxt(vp_path, vp_array)
