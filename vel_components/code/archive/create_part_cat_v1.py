# Downsample the particle catalogue for calculating the pairwise velocity
# v0.1.0, 2021-09-13 - Code started with snippets from calc_lin_vel.py

# Imports
import camb
import numpy as np
import os
from struct import iter_unpack, unpack_from
import sys


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


# omegam0_planck = omh2_planck / h_planck**2.
# H_z2 = H_z(0.7, h_planck*100., omegam0_planck)
# H_49 = H_z(49., h_planck*100., omegam0_planck)
# sigma8_planck = 0.830

# fsigma8_z2 = fsigma8_approximate(0.7, sigma8=sigma8_planck, omegam0=omegam0_planck)
# fsigma8_49 = fsigma8_approximate(49., sigma8=sigma8_planck, omegam0=omegam0_planck)


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

ic_type = sys.argv[1]
sample_size = int(sys.argv[2])
z = float(sys.argv[3])
print("IC Type:", ic_type)
print("Sample Size:", sample_size)

omch2_planck = 0.1199
ombh2_planck = 0.02222
omh2_planck = 0.14212
h_planck = 0.6726
ns_planck = 0.9652
sigma8_planck = 0.830
As_default = 2e-09
sigma8_default = 0.7917199
As_rescale = (sigma8_planck / sigma8_default) ** 2.

pars = camb.CAMBparams()
pars.set_cosmology(H0=h_planck*100., ombh2=ombh2_planck, omch2=omch2_planck)
pars.InitPower.set_params(ns=ns_planck, As=As_default*As_rescale)
if z != 0. and z != 49.:
    pars.set_matter_power(redshifts=[0., z, 49.], kmax=2.0)
elif z == 0.:
    pars.set_matter_power(redshifts=[z, 49.], kmax=2.0)
else:
    pars.set_matter_power(redshifts=[0., z], kmax=2.0)

pars.NonLinear = camb.model.NonLinear_none
results = camb.get_results(pars)
kh, zcamb, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints=200)
sigma8 = np.array(results.get_sigma8())
if z != 0. and z != 49.:
    print("Camb sigma8(z=0)=", sigma8[2])
else:
    print("Camb sigma8(z=0)=", sigma8[1])

H_z2 = results.hubble_parameter(z)
H_49 = results.hubble_parameter(49.)
if z != 49.:
    fsigma8_z2 = np.array(results.get_fsigma8())[1]
    sigma8_z2 = sigma8[1]
else:
    fsigma8_z2 = np.array(results.get_fsigma8())[0]
    sigma8_z2 = sigma8[0]
fsigma8_49 = np.array(results.get_fsigma8())[0]
sigma8_49 = sigma8[0]

del pars
del results

print("Camb H(z={}) = {}".format(z, H_z2))
print("Camb H(z=49) = {}".format(H_49))
print("Camb fsig8(z={}) = {}".format(z, fsigma8_z2))
print("Camb fsig8(z=49) = {}".format(fsigma8_49))
print()


if ic_type=="ic_no_growth_corr":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr"
    # disp_scaling = fsigma8_z2 / fsigma8_49
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / fsigma8_49 / h_planck
    disp_scaling = sigma8_z2 / sigma8_49
    # vel_scaling = H_z2 / (1. + z) * sigma8_z2 / sigma8_49

elif ic_type=="ic_ngc_nplt":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_ngc_nplt"
    # disp_scaling = fsigma8_z2 / fsigma8_49
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / fsigma8_49 / h_planck
    disp_scaling = sigma8_z2 / sigma8_49
    # vel_scaling = H_z2 / (1. + z) * sigma8_z2 / sigma8_49

elif ic_type=="ic_32c":
    ic_dir = "/home/mj3chapm/scratch/abacus_test/AbacusCosmos_1100box_products/planck_products/ic_32c"
    # disp_scaling = fsigma8_z2 / fsigma8_49
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / fsigma8_49 / h_planck
    disp_scaling = sigma8_z2 / sigma8_49
    # vel_scaling = H_z2 / (1. + z) * sigma8_z2 / sigma8_49

elif ic_type=="ic_no_growth_corr_z0":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr_z0"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

elif ic_type=="ic_no_growth_corr_z0.7":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr_z0.7"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

elif ic_type=="ic_ngc_nplt_z0.7":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_ngc_nplt_z0.7"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

elif ic_type=="ic_z0.7":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_z0.7"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

elif ic_type=="ic_no_growth_corr_z1.5":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr_z1.5"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

elif ic_type=="ic_no_growth_corr_z5":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr_z5"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

elif ic_type=="ic_no_growth_corr_z20":
    ic_dir = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/ic_no_growth_corr_z20"
    disp_scaling = 1.
    vel_scaling = H_z2 / (1. + z) * fsigma8_z2 / sigma8_z2 / h_planck

else:
    sys.exit("Unrecognized initial condition type!!")

cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/cosmo_params.dat")
if "z" in ic_type:
    vel_scaling_func = calc_vel_scaling(z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4], low_redshift=True)
else:
    vel_scaling_func =  calc_vel_scaling(z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4])

print("Velocity scaling:", vel_scaling)
print("Velocity scaling func:", vel_scaling_func)

seed = 1337
rng = np.random.default_rng(seed)
sample_indices = np.empty(0)
while sample_indices.size != sample_size:
    random_indices = rng.integers(1440**3, size=(sample_size - sample_indices.size))
    sample_indices = np.unique(np.concatenate((sample_indices, random_indices)))
    print("Try {}, Seed is {}, # of Indices is {}".format(seed - 1337 + 1, seed, sample_indices.size))
    seed += 1

last_sort_index = 0
tot_particles = 0

particle_sample = np.empty((sample_size, 10))
particle_index = 0

for i in range(375):
# for i in range(2):
    print("Starting IC file {}".format(i))
    tot_ic_size = os.path.getsize("{}/ic_{}".format(ic_dir, i))
    N_particles = tot_ic_size / 32 # Total size of objects is 30, but becomes 32 because of padding
    # ic_array = np.empty((int(N_particles), 9))
    with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
        bdata = file.read()
        # data = iter_unpack("3h6f", bdata)
        # for j, line in enumerate(data):
        #     ic_array[j] = line[-3:]

        for j in range(last_sort_index, sample_size):
            if (sample_indices[j] - tot_particles) < N_particles:
                particle_data = unpack_from("3h6f", bdata, offset=int((sample_indices[j] - tot_particles)*32))
                # print("Particle Data:")
                # print(particle_data)
                # print("Sample Index:", sample_indices[j], "Particle Index:", particle_data[0]*1440**2 + particle_data[1]*1440 + particle_data[2])
                particle_sample[particle_index, 0] = sample_indices[j]
                particle_sample[particle_index, 1:4] = particle_data[:3]
                for k in range(3):
                    particle_sample[particle_index, 4+k] = (particle_data[k] / 1440 * 1050.) + (particle_data[3+k] * disp_scaling)
                    particle_sample[particle_index, 7+k] = particle_data[6+k] * vel_scaling
                # print("Particle Sample:")
                # print(particle_sample[particle_index, :])
                particle_index += 1

            else:
                print("Switched IC at {}".format(j))
                last_sort_index = j
                break

    tot_particles += N_particles

header = "#ID\tI\tJ\tK\tX\tY\tZ\tVX\tVY\tVZ"
np.savetxt('/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_sample_z{}_{}_{}_v7.dat'.format(z, sample_size, ic_type), particle_sample, header=header)

