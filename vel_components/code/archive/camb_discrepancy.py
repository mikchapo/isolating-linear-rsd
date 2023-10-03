# Explore the discrepancy between my CAMB results and the Abacus CAMB results
# v0.1.0, 2021-10-10 - Code started with snippets from check_ic_norms.py

# Imports
import camb
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad


z1 = 0.
z2 = 49.
z3 = 0.7

# Load the cosmological parameters of the simulation box
cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                          "AbacusCosmos_1100box_products/"
                          "AbacusCosmos_1100box_00_products/"
                          "AbacusCosmos_1100box_00_rockstar_halos/"
                          "info/cosmo_params.dat")
ini_file = ("/home/mj3chapm/scratch/abacus/"
            "AbacusCosmos_1100box_products/"
            "AbacusCosmos_1100box_00_products/"
            "AbacusCosmos_1100box_00_rockstar_halos/"
            "info/camb_params.ini")

omegam0 = (cosmo_params[1] + cosmo_params[2]) / (cosmo_params[0] / 100.)**2

# Default CAMB power spectrum normalization
As_default = 2e-09

# List of redshifts required for CAMB calulcations without duplicates, to
# avoid CAMB error
redshifts = list(dict.fromkeys([z1, z2, z3, 0.]))
redshifts.sort()
redshifts.reverse()

# Initialize CAMB parameters for power spectrum calculation
pars = camb.CAMBparams()
pars.set_cosmology(H0=cosmo_params[0], omch2=cosmo_params[1],
                   ombh2=cosmo_params[2])
pars.DarkEnergy = camb.dark_energy.DarkEnergyFluid(w=cosmo_params[6])
pars.InitPower.set_params(ns=cosmo_params[3])
# pars.set_accuracy(AccuracyBoost=3.0, lSampleBoost=3.0, lAccuracyBoost=3.0,
#                   DoLateRadTruncation=True)
# pars.WantCls = False
# pars.WantScalars = False
# pars.DoLensing = False
# pars.MassiveNuMethod = "Nu_trunc"
# pars.DoLateRadTruncation = False
pars.set_matter_power(redshifts=redshifts, kmax=32)
pars.NonLinear = camb.model.NonLinear_none

# Get the value of sigma8(z=0) with the default normalization in order to
# calculate a new normalization
results = camb.get_results(pars)
sigma8_0_default = np.array(results.get_sigma8_0())
As_rescale = (cosmo_params[4] / sigma8_0_default) ** 2.

# Re-initalize parameters and get results with new normalization
pars.InitPower.set_params(ns=cosmo_params[3], As=As_default*As_rescale)
pars.set_matter_power(redshifts=redshifts, kmax=32)
# pars.set_matter_power(redshifts=redshifts, kmax=32)
pars.NonLinear = camb.model.NonLinear_none
results = camb.get_results(pars)

z1_index = redshifts.index(z1)
z2_index = redshifts.index(z2)
z3_index = redshifts.index(z3)

sigma8 = np.array(results.get_sigma8())
my_camb_fsigma8_z3 = np.array(results.get_fsigma8())[z3_index]
my_camb_H_z3 = results.hubble_parameter(z3)

IC_pars = camb.read_ini(ini_file)
IC_pars.set_matter_power(redshifts=redshifts, kmax=32)
IC_results = camb.get_results(IC_pars)
IC_sigma8_0 = np.array(IC_results.get_sigma8_0())
IC_camb_sigma8_z2 = (np.array(IC_results.get_sigma8())[z2_index] *
                     cosmo_params[4] / IC_sigma8_0)
IC_camb_fsigma8_z3 = (np.array(IC_results.get_fsigma8())[z3_index] *
                      cosmo_params[4] / IC_sigma8_0)
IC_camb_H_z3 = IC_results.hubble_parameter(z3)

print("My parameters:")
print(pars)

# print("Abacus parameters:")
# print(IC_pars)

print("Cosmology sigma8(z=0) =", cosmo_params[4])
print("My Camb sigma8(z={}) =".format(z1), sigma8[z1_index])
print("IC Camb sigma8(z=0) =", IC_sigma8_0, "\n")

print("IC sigma8(z=49.0) =", cosmo_params[-1])
print("My CAMB sigma8(z={}) =".format(z2), sigma8[z2_index])
print("IC CAMB sigma8(z={}) =".format(z2), IC_camb_sigma8_z2)
print("Ratio:", sigma8[z2_index] / IC_camb_sigma8_z2)

print("My CAMB fsigma8(z={}) =".format(z3), my_camb_fsigma8_z3)
print("IC CAMB fsigma8(z={}) =".format(z3), IC_camb_fsigma8_z3)
print("IC CAMB to My CAMB:", IC_camb_fsigma8_z3 / my_camb_fsigma8_z3, "\n")

print("My CAMB H(z={}) =".format(z3), my_camb_H_z3)
print("IC CAMB H(z={}) =".format(z3), IC_camb_H_z3)
print("IC CAMB to My CAMB:", IC_camb_H_z3 / my_camb_H_z3, "\n")
