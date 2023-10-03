# Compare high and low redshift output of the initial conditions code to
# determine why there is a velocity offset
# v0.1.0, 2021-11-09 - Code started, copied from check_plt.py

# Imports
import camb
import numpy as np
from scipy.integrate import quad
from struct import unpack_from


def D_integrand(a, omegam, omegal):
    return ((omegam / (a**3.) / (omegam / (a**3.) + omegal))**0.55 - 1.) / a


def calc_D(z, omegam0, omegal):
    a = 1. / (1. + z)
    integral = quad(D_integrand, 0., a, args=(omegam0, omegal))
    return a*np.exp(integral[0])


def calc_vel_scaling(z_ic, z_cat, H0, omch2, ombh2, ns, sigma8, camb_s8=True):
    '''
    Caclulates the scaling required for Abacus initial condition particle
    velocities to match the give redshift. Taken directly from
    lin_vel_funcs.py. ALL CHANGES MARKED WITH TAG [CHANGE]!!

    Input:
    ------
    z_ic - Redshift of the initial conditions
    z_cat - Redshift of the catalogue
    H0 - Hubble constant at z=0
    omch2 - Omega cold dark matter multiplied by h**2
    ombh2 - Omega baryon multiplied by h**2
    ns - The initial power spectrum slope
    sigma8 - RMS fluctuations in a sphere of 8 Mpc/h
    camb_s8 - [CHANGE] Boolean to determine whether to use CAMB or growth
               scaling for sigma8(z_ic)

    Note: This code assumes that z_ic > z_cat when indexing CAMB results.

    Output:
    -------
    vel_scaling - Factor used to scale the initial condition velocities
    '''

    # Default CAMB power spectrum normalization
    As_default = 2e-09

    # List of redshifts required for CAMB calulcations without duplicates, to
    # avoid CAMB error
    redshifts = list(dict.fromkeys([z_ic, z_cat, 0.]))

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

    # Print sigma8(z=0) from CAMB to make sure the As scaling was done properly
    sigma8_camb = np.array(results.get_sigma8_0())
    print("CAMB sigma8(z=0)=", sigma8_camb)

    # Get all the cosmological parameters needed for the scaling calculation
    z_ic_index = redshifts.index(z_ic)
    z_cat_index = redshifts.index(z_cat)
    H_z_cat = results.hubble_parameter(z_cat)
    # [CHANGE] Allow either CAMB sigma8(z_ic) or growth scaled
    omegam0 = (ombh2 + omch2) / (H0 / 100.)**2
    if camb_s8:
        sigma8_z_ic = np.array(results.get_sigma8())[z_ic_index]
    else:
        sigma8_z_ic = sigma8 * (calc_D(z_ic, omegam0, 1 - omegam0) /
                                calc_D(0., omegam0, 1 - omegam0))

    fsigma8_z_cat = np.array(results.get_fsigma8())[z_cat_index]

    '''
    The velocity scaling from z1 to z2 should be the ratio of
    aHfsig8(z2)/aHfsig8(z1) (see Eq.4.77 of Mo, van den Bosch, and White).
    However, the particle 'velocities' output by the zeldovich-PLT code are
    actually co-moving displacements with units of Mpc/h (using the code
    defaults). To convert these co-moving displacements to velocities we need
    to multiply by aHf(z1)/h, which cancel with several of the factors from the
    linear theory redshift scaling. The result is given below.

    Note: If z1=z2 this works out to the same expression.
    '''

    vel_scaling = (H_z_cat / (1. + z_cat) * fsigma8_z_cat / sigma8_z_ic /
                   (H0 / 100.))
    return vel_scaling


# Load the cosmological parameters of the simulation box
cosmo_params = np.loadtxt("AbacusCosmos_1100box_planck_FoF_halos/info/"
                          "cosmo_params.dat")

with open("ic_ngc_nplt_z49.0/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)

    z_49_vel_scaling_camb_s8 = calc_vel_scaling(49., 0.7, cosmo_params[0],
                                                cosmo_params[1],
                                                cosmo_params[2],
                                                cosmo_params[3],
                                                cosmo_params[4])
    z_49_growth_ic_v_camb_scale = (np.array(particle_data[-3:]) *
                                   z_49_vel_scaling_camb_s8)

    z_49_vel_scaling_growth_s8 = calc_vel_scaling(49., 0.7, cosmo_params[0],
                                                  cosmo_params[1],
                                                  cosmo_params[2],
                                                  cosmo_params[3],
                                                  cosmo_params[4],
                                                  camb_s8=False)
    z_49_growth_ic_v_growth_scale = (np.array(particle_data[-3:]) *
                                     z_49_vel_scaling_growth_s8)

with open("ic_ngc_nplt_z0.7/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)

    z_0p7_vel_scaling_camb_s8 = calc_vel_scaling(0.7, 0.7, cosmo_params[0],
                                                 cosmo_params[1],
                                                 cosmo_params[2],
                                                 cosmo_params[3],
                                                 cosmo_params[4])
    z_0p7_growth_ic_v_camb_scale = (np.array(particle_data[-3:]) *
                                    z_0p7_vel_scaling_camb_s8)

    z_0p7_vel_scaling_growth_s8 = calc_vel_scaling(0.7, 0.7, cosmo_params[0],
                                                   cosmo_params[1],
                                                   cosmo_params[2],
                                                   cosmo_params[3],
                                                   cosmo_params[4],
                                                   camb_s8=False)
    z_0p7_growth_ic_v_growth_scale = (np.array(particle_data[-3:]) *
                                      z_0p7_vel_scaling_growth_s8)

with open("ic_ngc_nplt_z49.0_camb_s8/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)

    z_49_vel_scaling_camb_s8 = calc_vel_scaling(49., 0.7, cosmo_params[0],
                                                cosmo_params[1],
                                                cosmo_params[2],
                                                cosmo_params[3],
                                                cosmo_params[4])
    z_49_camb_ic_v_camb_scale = (np.array(particle_data[-3:]) *
                                 z_49_vel_scaling_camb_s8)

    z_49_vel_scaling_growth_s8 = calc_vel_scaling(49., 0.7, cosmo_params[0],
                                                  cosmo_params[1],
                                                  cosmo_params[2],
                                                  cosmo_params[3],
                                                  cosmo_params[4],
                                                  camb_s8=False)
    z_49_camb_ic_v_growth_scale = (np.array(particle_data[-3:]) *
                                   z_49_vel_scaling_growth_s8)

with open("ic_ngc_nplt_z0.7_camb_s8/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)

    z_0p7_vel_scaling_camb_s8 = calc_vel_scaling(0.7, 0.7, cosmo_params[0],
                                                 cosmo_params[1],
                                                 cosmo_params[2],
                                                 cosmo_params[3],
                                                 cosmo_params[4])
    z_0p7_camb_ic_v_camb_scale = (np.array(particle_data[-3:]) *
                                  z_0p7_vel_scaling_camb_s8)

    z_0p7_vel_scaling_growth_s8 = calc_vel_scaling(0.7, 0.7, cosmo_params[0],
                                                   cosmo_params[1],
                                                   cosmo_params[2],
                                                   cosmo_params[3],
                                                   cosmo_params[4],
                                                   camb_s8=False)
    z_0p7_camb_ic_v_growth_scale = (np.array(particle_data[-3:]) *
                                    z_0p7_vel_scaling_growth_s8)

print("\n", "High vs low redshift comparisons with the same method:", "\n")

print("v(z=49), Growth IC CAMB Scale / v(z=0.7), Growth IC CAMB Scale:")
print(z_49_growth_ic_v_camb_scale / z_0p7_growth_ic_v_camb_scale, "\n")

print("v(z=49), Growth IC Growth Scale / v(z=0.7), Growth IC Growth Scale:")
print(z_49_growth_ic_v_growth_scale / z_0p7_growth_ic_v_growth_scale, "\n")

print("v(z=49), CAMB IC CAMB Scale / v(z=0.7), CAMB IC CAMB Scale:")
print(z_49_camb_ic_v_camb_scale / z_0p7_camb_ic_v_camb_scale, "\n")

print("v(z=49), CAMB IC Growth Scale / v(z=0.7), CAMB IC Growth Scale:")
print(z_49_camb_ic_v_growth_scale / z_0p7_camb_ic_v_growth_scale, "\n")

print("\n", "Same redshift comparisons between different methods:", "\n")

print("v(z=49), Growth IC CAMB Scale / v(z=49), Growth IC Growth Scale:")
print(z_49_growth_ic_v_camb_scale / z_49_growth_ic_v_growth_scale)
print("v(z=0.7), Growth IC CAMB Scale / v(z=0.7), Growth IC Growth Scale:")
print(z_0p7_growth_ic_v_camb_scale / z_0p7_growth_ic_v_growth_scale, "\n")

print("v(z=49), CAMB IC CAMB Scale / v(z=49), Growth IC Growth Scale:")
print(z_49_camb_ic_v_camb_scale / z_49_growth_ic_v_growth_scale)
print("v(z=0.7), CAMB IC CAMB Scale / v(z=0.7), Growth IC Growth Scale:")
print(z_0p7_camb_ic_v_camb_scale / z_0p7_growth_ic_v_growth_scale, "\n")

print("v(z=49), CAMB IC Growth Scale / v(z=49), Growth IC Growth Scale:")
print(z_49_camb_ic_v_growth_scale / z_49_growth_ic_v_growth_scale)
print("v(z=0.7), CAMB IC Growth Scale / v(z=0.7), Growth IC Growth Scale:")
print(z_0p7_camb_ic_v_growth_scale / z_0p7_growth_ic_v_growth_scale, "\n")

# Change Log
