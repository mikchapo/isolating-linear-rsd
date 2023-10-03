# Check the z=0.7 sigma8 value used to normalize the initial conditions
# for each simulation box
# v0.1.0, 2021-10-10 - Code started with snippets from check_ic_norms.py

# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


def D_integrand(a, omegam, omegal):
    return ((omegam / (a**3.) / (omegam / (a**3.) + omegal))**0.55 - 1.) / a


def calc_D(z, omegam0, omegal):
    a = 1. / (1. + z)
    integral = quad(D_integrand, 0., a, args=(omegam0, omegal))
    return a*np.exp(integral[0])


z1 = 0.
z2 = 0.7

load_ratios = False

if load_ratios:
    s8_ratios = np.loadtxt("../output/data_products/"
                           "check_ic_norms_z0.7_s8_ratios.dat")
    fs8_ratios = np.loadtxt("../output/data_products/"
                            "check_ic_norms_z0.7_fs8_ratios.dat")
    H_ratios = np.loadtxt("../output/data_products/"
                          "check_ic_norms_z0.7_H_ratios.dat")

else:
    s8_ratios = np.empty((41, 2))
    fs8_ratios = np.empty(41)
    H_ratios = np.empty(41)

    for i in range(41):
        print("\nStarting Box {}".format(i), "\n")

        # Load the cosmological parameters of the simulation box
        if i == 40:
            cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                      "AbacusCosmos_1100box_products/"
                                      "AbacusCosmos_1100box_planck_products/"
                                      "AbacusCosmos_1100box_planck_rockstar_"
                                      "halos/info/cosmo_params.dat")
            ini_file = ("/home/mj3chapm/scratch/abacus/"
                        "AbacusCosmos_1100box_products/"
                        "AbacusCosmos_1100box_planck_products/"
                        "AbacusCosmos_1100box_planck_rockstar_halos/"
                        "info/camb_params.ini")

        elif i > 9:
            cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                      "AbacusCosmos_1100box_products/"
                                      "AbacusCosmos_1100box_{}_products/"
                                      "AbacusCosmos_1100box_{}_rockstar_halos/"
                                      "info/cosmo_params.dat".format(i, i))
            ini_file = ("/home/mj3chapm/scratch/abacus/"
                        "AbacusCosmos_1100box_products/"
                        "AbacusCosmos_1100box_{}_products/"
                        "AbacusCosmos_1100box_{}_rockstar_halos/"
                        "info/camb_params.ini".format(i, i))

        else:
            cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                                      "AbacusCosmos_1100box_products/"
                                      "AbacusCosmos_1100box_0{}_products/"
                                      "AbacusCosmos_1100box_0{}_rockstar_"
                                      "halos/info/cosmo_params.dat".format(i,
                                                                           i))
            ini_file = ("/home/mj3chapm/scratch/abacus/"
                        "AbacusCosmos_1100box_products/"
                        "AbacusCosmos_1100box_0{}_products/"
                        "AbacusCosmos_1100box_0{}_rockstar_halos/"
                        "info/camb_params.ini".format(i, i))

        omegam0 = ((cosmo_params[1] + cosmo_params[2]) /
                   (cosmo_params[0] / 100.)**2)
        sigma8_z2 = cosmo_params[4] * (calc_D(z2, omegam0, 1 - omegam0) /
                                       calc_D(z1, omegam0, 1 - omegam0))

        # Default CAMB power spectrum normalization
        As_default = 2e-09

        # List of redshifts required for CAMB calulcations without duplicates,
        # to avoid CAMB error
        redshifts = list(dict.fromkeys([z1, z2, 0.]))
        redshifts.sort()
        redshifts.reverse()

        # Initialize CAMB parameters for power spectrum calculation
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=cosmo_params[0], omch2=cosmo_params[1],
                           ombh2=cosmo_params[2], num_massive_neutrinos=0,
                           mnu=0)
        pars.DarkEnergy = camb.dark_energy.DarkEnergyFluid(w=cosmo_params[6])
        pars.InitPower.set_params(ns=cosmo_params[3])
        pars.set_matter_power(redshifts=redshifts, kmax=32.0)
        pars.NonLinear = camb.model.NonLinear_none

        # Get the value of sigma8(z=0) with the default normalization in order
        # to calculate a new normalization
        results = camb.get_results(pars)
        sigma8_0_default = np.array(results.get_sigma8_0())
        As_rescale = (cosmo_params[4] / sigma8_0_default) ** 2.

        # Re-initalize parameters and get results with new normalization
        pars.InitPower.set_params(ns=cosmo_params[3], As=As_default*As_rescale)
        pars.set_matter_power(redshifts=redshifts, kmax=32.0)
        pars.NonLinear = camb.model.NonLinear_none
        results = camb.get_results(pars)

        z1_index = redshifts.index(z1)
        z2_index = redshifts.index(z2)

        sigma8 = np.array(results.get_sigma8())
        my_camb_fsigma8 = np.array(results.get_fsigma8())[z2_index]
        my_camb_H = results.hubble_parameter(z2)

        IC_pars = camb.read_ini(ini_file)
        IC_pars.set_matter_power(redshifts=redshifts, kmax=32)
        IC_results = camb.get_results(IC_pars)
        IC_sigma8_0 = np.array(IC_results.get_sigma8_0())
        IC_camb_sigma8_z2 = (np.array(IC_results.get_sigma8())[z2_index] *
                             cosmo_params[4] / IC_sigma8_0)
        IC_camb_fsigma8 = (np.array(IC_results.get_fsigma8())[z2_index] *
                           cosmo_params[4] / IC_sigma8_0)
        IC_camb_H = IC_results.hubble_parameter(z2)

        s8_ratios[i, 0] = sigma8_z2 / sigma8[z2_index]
        s8_ratios[i, 1] = IC_camb_sigma8_z2 / sigma8[z2_index]
        fs8_ratios[i] = IC_camb_fsigma8 / my_camb_fsigma8
        H_ratios[i] = IC_camb_H / my_camb_H

        print("Cosmology sigma8(z=0) =", cosmo_params[4])
        print("My Camb sigma8(z={}) =".format(z1), sigma8[z1_index])
        print("IC Camb sigma8(z=0) =", IC_sigma8_0, "\n")

        print("Growth scaled sigma8(z={}) =".format(z2), sigma8_z2)
        print("My CAMB sigma8(z={}) =".format(z2), sigma8[z2_index])
        print("IC CAMB sigma8(z={}) =".format(z2), IC_camb_sigma8_z2)
        print("Growth to My CAMB:", sigma8_z2 / sigma8[z2_index])
        print("IC CAMB to My CAMB:", IC_camb_sigma8_z2 / sigma8[z2_index],
              "\n")

        print("My CAMB fsigma8(z={}) =".format(z2), my_camb_fsigma8)
        print("IC CAMB fsigma8(z={}) =".format(z2), IC_camb_fsigma8)
        print("IC CAMB to My CAMB:", IC_camb_fsigma8 / my_camb_fsigma8, "\n")

        print("My CAMB H(z={}) =".format(z2), my_camb_H)
        print("IC CAMB H(z={}) =".format(z2), IC_camb_H)
        print("IC CAMB to My CAMB:", IC_camb_H / my_camb_H, "\n")

    np.savetxt("../output/data_products/check_ic_norms_z0.7_s8_ratios.dat",
               s8_ratios)
    np.savetxt("../output/data_products/check_ic_norms_z0.7_fs8_ratios.dat",
               fs8_ratios)
    np.savetxt("../output/data_products/check_ic_norms_z0.7_H_ratios.dat",
               H_ratios)

plt.figure(figsize=(8., 6.), dpi=300)
plt.axvline(x=1., linestyle="-", color="k")
# plt.hist(s8_ratios[:, 0], bins=10, color="C0", alpha=1,
#          label="Growth scaled sigma8(z=0.7)")
plt.hist(s8_ratios[:, 1], bins=10, color="C2", alpha=1,
         label="IC CAMB sigma8(z=0.7)")
plt.xlabel(r"$\sigma_8(z=0.7) / \sigma_8^{CAMB}(z=0.7)$")
plt.ylabel("Count")
plt.legend()
plt.savefig("../output/plots/check_ic_norms_z0.7_s8_hist.jpg")

plt.figure(figsize=(8., 6.), dpi=300)
plt.axvline(x=1., linestyle="-", color="k")
plt.hist(fs8_ratios[:], bins=10, color="C2", alpha=1,
         label="IC CAMB fsigma8(z=0.7)")
plt.xlabel(r"$f\sigma_8(z=0.7) / f\sigma_8^{CAMB}(z=0.7)$")
plt.ylabel("Count")
plt.legend()
plt.savefig("../output/plots/check_ic_norms_z0.7_fs8_hist.jpg")

plt.figure(figsize=(8., 6.), dpi=300)
plt.axvline(x=1., linestyle="-", color="k")
plt.hist(H_ratios[:], bins=10, color="C2", alpha=1,
         label="IC CAMB sigma8(z=0.7)")
plt.xlabel(r"$H(z=0.7) / H^{CAMB}(z=0.7)$")
plt.ylabel("Count")
plt.legend()
plt.savefig("../output/plots/check_ic_norms_z0.7_H_hist.jpg")
