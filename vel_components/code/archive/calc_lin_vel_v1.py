# Calculate the linear velocity component of halos
# v0.2.0, 2021-10-22 - Updated the velocity sclaing to match the final version
#                      of the particle tests

# Imports
import camb
import numpy as np
import os
from struct import iter_unpack
import sys

from AbacusCosmos import Halos


def H_z(z, H0, omegam0):
    return H0 * np.sqrt(omegam0 * (1. + z) * (1. + z) * (1. + z) + (1. - omegam0))


def omegam(z, omegam0=0.315):
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
    l = 0.83 + 0.19 / (omegam0**0.8)
    beta = 0.98 + 0.01 / (omegam0**1.2)
    gamma = 0.64 + 0.09 / (omegam0**0.4)
    # print("lambda:", l, "beta:", beta, "gamma:", gamma)
    return l * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


def calc_vel_scaling(z, H0, omch2, ombh2, ns, sigma8, low_redshift=False):
    # Default CAMB power spectrum normalization
    As_default = 2e-09

    # List of redshifts required for CAMB calulcations without duplicates, to
    # avoid CAMB error
    redshifts = list(dict.fromkeys([49., z, 0.]))

    # Initialize CAMB parameters for power spectrum calculation
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
    pars.InitPower.set_params(ns=ns)
    pars.set_matter_power(redshifts=redshifts, kmax=2.0)
    pars.NonLinear = camb.model.NonLinear_none

    # Get the value of sigma8(z=0) with the default normalization in order to
    # calculate a new normalization
    results = camb.get_results(pars)
    sigma8_0_default = np.array(results.get_sigma8_0())
    As_rescale = (sigma8 / sigma8_0_default) ** 2.

    # Re-initalize parameters and get results with new normalization
    pars.InitPower.set_params(ns=ns, As=As_default*As_rescale)
    pars.set_matter_power(redshifts=redshifts, kmax=2.0)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)

    sigma8_camb = np.array(results.get_sigma8_0())
    print("CAMB sigma8(z=0)=", sigma8_camb)

    H_z2 = results.hubble_parameter(z)
    H_49 = results.hubble_parameter(49.)
    z_index = redshifts.index(z)
    fsigma8_z2 = np.array(results.get_fsigma8())[z_index]
    sigma8_z2 = np.array(results.get_sigma8())[z_index]
    fsigma8_49 = np.array(results.get_fsigma8())[0]
    sigma8_49 = np.array(results.get_sigma8())[0]

    if low_redshift:
        vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / (H0 / 100.)

    else:
        vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / fsigma8_49 / (H0 / 100.)

    return vel_scaling

cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/cosmo_params.dat")

omh2_planck = 0.14212
h_planck = 0.6726
omegam0_planck = omh2_planck / h_planck**2.
H_0p7 = H_z(0.7, h_planck*100., omegam0_planck)
H_49 = H_z(49., h_planck*100., omegam0_planck)
sigma8_planck = 0.830

fsigma8_0p7 = fsigma8_approximate(0.7, sigma8=sigma8_planck, omegam0=omegam0_planck)
fsigma8_49 = fsigma8_approximate(49., sigma8=sigma8_planck, omegam0=omegam0_planck)


ic_type = sys.argv[1]
print("IC Type:", ic_type)

if ic_type=="low_redshift":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr_z0.7"
    z = 0.7
    vel_scaling = calc_vel_scaling(z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4], low_redshift=True)

elif ic_type=="ic_no_growth_corr":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr"
    z = 0.7
    vel_scaling = calc_vel_scaling(z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4])

else:
    ic_dir = "/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/ic_32c"
    z = 0.7
    vel_scaling = calc_vel_scaling(z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4])


cat = Halos.make_catalog_from_dir(dirname='/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700',
                                  load_subsamples=True, load_pids=True)
halos = cat.halos

# if type(halos) is np.ndarray:
#   print("Halos stored in numpy array :)")
# else:
#   print("Halos are not stored as a numpy array :'(")

subsamples = cat.subsamples
pids = subsamples["pid"]
# pids = cat.subsample_pids

# pids_0 = pids[halos[0]['subsamp_start']:halos[0]['subsamp_start']+halos[0]['subsamp_len']]
# print("IDs of halo 0:")
# print(pids_0)
# print("Total number of particles: {}".format(halos[0]['N']))
# print("Total number of subsampled particles: {}".format(pids_0.size))

subsample_array = np.empty((pids.size, 4))
subsample_array[:, 0] = pids

sort_map = np.argsort(pids)
sorted_pids = pids[sort_map]

last_sort_pid_index = 0
tot_particles = 0

for i in range(375):
# for i in range(2):
    print("Starting IC file {}".format(i))
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    N_particles = tot_ic_size / 32 # Total size of objects is 30, but becomes 32 because of padding
    ic_array = np.empty(int(N_particles), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
    with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
        bdata = file.read()
        data = iter_unpack("3h6f", bdata)
        for j, line in enumerate(data):
            ic_array[j] = line[-3:]

    for j in range(last_sort_pid_index, pids.size):
        try:
            subsample_array[sort_map[j], 1] = ic_array["vx"][int(sorted_pids[j] - tot_particles)]
            subsample_array[sort_map[j], 2] = ic_array["vy"][int(sorted_pids[j] - tot_particles)]
            subsample_array[sort_map[j], 3] = ic_array["vz"][int(sorted_pids[j] - tot_particles)]

        except IndexError:
            print("Hit IndexError at {}".format(j))
            last_sort_pid_index = j
            break

    tot_particles += N_particles


# print("Halos Shape:", halos.shape)

# for field in sorted(halos[0].dtype.fields):
#     print(field, ':', halos[0][field])

halo_array = np.empty((halos.shape[0], 21))
halo_array[:, 0] = halos["id"]
halo_array[:, 1] = -1
halo_array[:, 2] = halos["N"] * 4.e10
halo_array[:, 3] = halos["vcirc_max"]
# halo_array[:, 4] = halos["sigma_v"]
halo_array[:, 4] = 0
halo_array[:, 5] = halos["r90"]
halo_array[:, 6] = 0
halo_array[:, 7] = halos["N"]
halo_array[:, 14] = -1

mask = np.ones(halos.shape[0], dtype=bool)
for i in range(halos.shape[0]):
    if halos[i]["subsamp_len"] > 0:
        halo_array[i, 8] = halos[i]["pos"][0]
        halo_array[i, 9] = halos[i]["pos"][1]
        halo_array[i, 10] = halos[i]["pos"][2]
        halo_array[i, 11] = halos[i]["vel"][0]
        halo_array[i, 12] = halos[i]["vel"][1]
        halo_array[i, 13] = halos[i]["vel"][2]
        halo_array[i, 15] = np.mean(subsample_array[halos[i]["subsamp_start"]:halos[i]["subsamp_start"]+halos[i]["subsamp_len"], 1]) * vel_scaling
        halo_array[i, 16] = np.mean(subsample_array[halos[i]["subsamp_start"]:halos[i]["subsamp_start"]+halos[i]["subsamp_len"], 2]) * vel_scaling
        halo_array[i, 17] = np.mean(subsample_array[halos[i]["subsamp_start"]:halos[i]["subsamp_start"]+halos[i]["subsamp_len"], 3]) * vel_scaling
        halo_array[i, 18] = halo_array[i, 11] - halo_array[i, 15]
        halo_array[i, 19] = halo_array[i, 12] - halo_array[i, 16]
        halo_array[i, 20] = halo_array[i, 13] - halo_array[i, 17]
    else:
        mask[i] = False

print("N_skipped:", halos.shape[0] - np.sum(mask))
halo_array = halo_array[mask, :]

header = "#ID DescID M200b Vmax Vrms R200b Rs Np X Y Z VX VY VZ Parent_ID VX_LIN VY_LIN VZ_LIN VX_NL VY_NL VZ_NL"
np.savetxt('/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/halos_lin_vel_{}.dat'.format(ic_type), halo_array, header=header)


# Change Log
# v0.1.1, 2021-08-12 - Updated the velocity scaling
# v0.1.0, 2021-07-20 - Code started



