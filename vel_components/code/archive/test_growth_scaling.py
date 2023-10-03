import numpy as np
import os
from struct import iter_unpack

z1 = 49.
z2 = 20.
ic_z2_options = ["ic_no_growth_corr_z0", "ic_no_growth_corr_z0.7", "ic_no_growth_corr_z1.5", "ic_no_growth_corr_z5", "ic_no_growth_corr_z20"]
# ic_z2_options = ["ic_ngc_nplt_z0.7", "ic_ngc_nplt_z20"]
ic_z2_dir = ic_z2_options[-1]

ic_z1_dir = "ic_no_growth_corr"
# ic_z1_dir = "ic_ngc_nplt"
ic_z1_size = os.path.getsize("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/{}/ic_0".format(ic_z1_dir))
N_particles_z1 = ic_z1_size / 32 # Total size of objects is 30, but becomes 32 because of padding
ic_z1_array = np.empty(int(N_particles_z1), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
with open("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/{}/ic_0".format(ic_z1_dir), "rb") as file:
    bdata = file.read()
    data = iter_unpack("3h6f", bdata)
    for j, line in enumerate(data):
        ic_z1_array[j] = line[-3:]
        if j > 20:
        	break

# ic_z1_alt_dir = "ic_no_growth_corr"
# ic_z1_alt_size = os.path.getsize("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/{}/ic_0".format(ic_z1_alt_dir))
# N_particles_z1_alt = ic_z1_alt_size / 32 # Total size of objects is 30, but becomes 32 because of padding
# ic_z1_alt_array = np.empty(int(N_particles_z1_alt), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
# with open("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/{}/ic_0".format(ic_z1_alt_dir), "rb") as file:
#     bdata = file.read()
#     data = iter_unpack("3h6f", bdata)
#     for j, line in enumerate(data):
#         ic_z1_alt_array[j] = line[-3:]
#         if j > 20:
#         	break

ic_z2_size = os.path.getsize("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/{}/ic_0".format(ic_z2_dir))
N_particles_z2 = ic_z2_size / 32 # Total size of objects is 30, but becomes 32 because of padding
ic_z2_array = np.empty(int(N_particles_z2), dtype=([('vx', 'f'), ('vy', 'f'), ('vz', 'f')]))
with open("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/{}/ic_0".format(ic_z2_dir), "rb") as file:
    bdata = file.read()
    data = iter_unpack("3h6f", bdata)
    for j, line in enumerate(data):
        ic_z2_array[j] = line[-3:]
        if j > 20:
        	break


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

omh2_planck = 0.14212
h_planck = 0.6726
omegam0_planck = omh2_planck / h_planck**2.
print("\nIntermediate steps:")
print("omegam0:", omegam0_planck)
H_z2 = H_z(z2, h_planck*100., omegam0_planck)
H_z1 = H_z(z1, h_planck*100., omegam0_planck)
print("H(z2):", H_z2)
print("H(z1):", H_z1)
sigma8_planck = 0.830

fsigma8_z2 = fsigma8_approximate(z2, sigma8=sigma8_planck, omegam0=omegam0_planck)
fsigma8_z1 = fsigma8_approximate(z1, sigma8=sigma8_planck, omegam0=omegam0_planck)

print("fsigma8_z2:", fsigma8_z2)
print("fsigma8_z1:", fsigma8_z1)
print("IC z1 first entry:", ic_z1_array[0])
# print("IC z1_alt first entry:", ic_z1_alt_array[0])
print("IC z2 first entry:", ic_z2_array[0], "\n")

vel_scaling = H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * fsigma8_z2 / fsigma8_z1

print("My results")
print("Ratio x, particle 0:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vx"] / ic_z1_array[0]["vx"])
print("Ratio y, particle 0:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vy"] / ic_z1_array[0]["vy"])
print("Ratio z, particle 0:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vz"] / ic_z1_array[0]["vz"])
print("Ratio x, particle 13:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[13]["vx"] / ic_z1_array[13]["vx"])
print("Ratio y, particle 13:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[13]["vy"] / ic_z1_array[13]["vy"])
print("Ratio z, particle 13:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[13]["vz"] / ic_z1_array[13]["vz"])
print("My scaling:", vel_scaling, "\n")

# Cosmology calculation with Omega_R=0, z2=0.7
H_z2 = 100.411
H_z1 = 13325.4
fsigma8_z2 = 0.472727
fsigma8_z1 = 0.0210862
vel_scaling = H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * fsigma8_z2 / fsigma8_z1

print("Cosmology calculator results:")
print("Ratio x:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vx"] / ic_z1_array[0]["vx"])
print("Ratio y:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vy"] / ic_z1_array[0]["vy"])
print("Ratio z:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vz"] / ic_z1_array[0]["vz"])
print("My scaling:", vel_scaling, "\n")

# Cosmology calculation with Omega_R != 0, z2=0.7
H_z2 = 100.436
H_z1 = 13426.8
fsigma8_z2 = 0.473348
fsigma8_z1 = 0.0213691
vel_scaling = H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * fsigma8_z2 / fsigma8_z1

print("Cosmology calculator results with Omega_R != 0:")
print("Ratio x:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vx"] / ic_z1_array[0]["vx"])
print("Ratio y:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vy"] / ic_z1_array[0]["vy"])
print("Ratio z:", H_z2 / (1. + z2) / (H_z1 / (1. + z1)) * ic_z2_array[0]["vz"] / ic_z1_array[0]["vz"])
print("My scaling:", vel_scaling, "\n")

# Calculation with some Omegam_R, z2=20
# H_z2 = 3638.68
# H_z1 = 13423.6
# fsigma8_z2 = 0.0505205
# fsigma8_z1 = 0.0213716