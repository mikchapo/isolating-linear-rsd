# Check the high redshift sigma8 value used to normalize the initial conditions
# for each simulation box
# v0.1.0, 2021-10-09 - Code started with snippets from test_simga8_rescaling.py

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
z2 = 49.

load_ratios = True

if load_ratios:
    ratios = np.loadtxt("../output/data_products/check_ic_norms_ratios.dat")

else:
    ratios = np.empty((41, 3))

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
                           ombh2=cosmo_params[2])
        pars.DarkEnergy = camb.dark_energy.DarkEnergyFluid(w=cosmo_params[6])
        pars.InitPower.set_params(ns=cosmo_params[3])
        pars.set_matter_power(redshifts=redshifts, kmax=2.0)
        pars.NonLinear = camb.model.NonLinear_none

        # Get the value of sigma8(z=0) with the default normalization in order
        # to calculate a new normalization
        results = camb.get_results(pars)
        sigma8_0_default = np.array(results.get_sigma8_0())
        As_rescale = (cosmo_params[4] / sigma8_0_default) ** 2.

        # Re-initalize parameters and get results with new normalization
        pars.InitPower.set_params(ns=cosmo_params[3], As=As_default*As_rescale)
        pars.set_matter_power(redshifts=redshifts, kmax=2.0)
        pars.NonLinear = camb.model.NonLinear_none
        results = camb.get_results(pars)

        sigma8 = np.array(results.get_sigma8())
        z1_index = redshifts.index(z1)
        z2_index = redshifts.index(z2)

        IC_pars = camb.read_ini(ini_file)
        IC_pars.set_matter_power(redshifts=redshifts, kmax=32)
        IC_results = camb.get_results(IC_pars)
        IC_sigma8_0 = np.array(IC_results.get_sigma8_0())
        IC_camb_sigma8_z2 = (np.array(IC_results.get_sigma8())[z2_index] *
                             cosmo_params[4] / IC_sigma8_0)

        ratios[i, 0] = sigma8_z2 / cosmo_params[-1]
        ratios[i, 1] = sigma8[z2_index] / cosmo_params[-1]
        ratios[i, 2] = IC_camb_sigma8_z2 / cosmo_params[-1]

        print("Cosmology sigma8(z=0) =", cosmo_params[4])
        print("My Camb sigma8(z={}) =".format(z1), sigma8[z1_index])
        print("IC Camb sigma8(z=0) =", IC_sigma8_0, "\n")

        print("IC sigma8(z=49.0) =", cosmo_params[-1])
        print("Growth scaled sigma8(z={}) =".format(z2), sigma8_z2)
        print("My CAMB sigma8(z={}) =".format(z2), sigma8[z2_index])
        print("IC CAMB sigma8(z={}) =".format(z2), IC_camb_sigma8_z2)
        print("Growth to IC ratio:", sigma8_z2 / cosmo_params[-1])
        print("My CAMB to IC ratio:", sigma8[z2_index] / cosmo_params[-1])
        print("IC CAMB to IC ratio:", IC_camb_sigma8_z2 / cosmo_params[-1])

    np.savetxt("../output/data_products/check_ic_norms_ratios.dat", ratios)

plt.figure(figsize=(8., 6.), dpi=300)
plt.axvline(x=1., linestyle="-", color="k")
plt.hist(ratios[:, 0], bins=10, color="C0",
         label="Growth scaled sigma8(z=49)")
plt.hist(ratios[:, 1], bins=10, color="C1",
         label="My CAMB sigma8(z=49)")
plt.hist(ratios[:, 2], bins=10, color="C2",
         label="IC CAMB sigma8(z=49)")
plt.xlabel(r"$\sigma_8(z=49) / \sigma_8^{IC}(z=49)$")
plt.ylabel("Count")
plt.legend()
plt.savefig("../output/plots/check_ic_norms_hist.jpg")
