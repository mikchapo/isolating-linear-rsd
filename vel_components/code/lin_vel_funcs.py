# Functions for calculating the linear velocity
# v0.2.0, 2021-11-11 - Updated calc_vel_scaling, saving old versions as v2, for
#                      more accurate cosmological parameters

# Imports
import camb
import numpy as np


def calc_vel_scaling(z_ic, z_cat, info_path, calc_disp_scaling=False):
    '''
    Caclulates the scaling required for Abacus initial condition particle
    velocities to match the give redshift.

    Input:
    ------
    z_ic - Redshift of the initial conditions
    z_cat - Redshift of the catalogue
    info_path - Path to the directory containing the cosmo_params.dat and
                camb_params.ini files

    Note: This code assumes that z_ic > z_cat when indexing CAMB results.

    Output:
    -------
    vel_scaling - Factor used to scale the initial condition velocities
    '''

    cosmo_params = np.loadtxt("{}/cosmo_params.dat".format(info_path))
    ini_file = ("{}/camb_params.ini".format(info_path))

    # Default CAMB power spectrum normalization
    As_default = 2e-09

    # List of redshifts required for CAMB calulcations without duplicates, to
    # avoid CAMB error and put them in the correct order for CAMB
    # (highest to lowest)
    redshifts = list(dict.fromkeys([z_ic, z_cat, 0.]))
    redshifts.sort()
    redshifts.reverse()

    # Changed to reading from Abacus CAMB .ini file for consistency with
    # parameters, in particular switching to massless neutrinos which give
    # 0.7% offset
    pars = camb.read_ini(ini_file)
    pars.set_matter_power(redshifts=redshifts, kmax=32)
    results = camb.get_results(pars)
    sigma8_0 = np.array(results.get_sigma8_0())
    sigma8_scaling = cosmo_params[4] / sigma8_0

    # Get all the cosmological parameters needed for the scaling calculation
    z_ic_index = redshifts.index(z_ic)
    z_cat_index = redshifts.index(z_cat)
    H_z_cat = results.hubble_parameter(z_cat)
    fsigma8_z_cat = (np.array(results.get_fsigma8())[z_cat_index] *
                     sigma8_scaling)
    # For sigma8 at the redshift of the inititial conditions use the actual
    # value from the abacus.par file, since it is biased low by 1% compared to
    # the CAMB output
    sigma8_ic = cosmo_params[-1]

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

    vel_scaling = (H_z_cat / (1. + z_cat) * fsigma8_z_cat / sigma8_ic /
                   (cosmo_params[0] / 100.))

    if calc_disp_scaling:
        sigma8_z_cat = (np.array(results.get_sigma8())[z_cat_index] *
                        sigma8_scaling)
        disp_scaling = sigma8_z_cat / sigma8_ic
        return vel_scaling, disp_scaling

    else:
        return vel_scaling


def calc_vel_scaling_v2(z_ic, z_cat, H0, omch2, ombh2, ns, sigma8):
    '''
    Caclulates the scaling required for Abacus initial condition particle
    velocities to match the give redshift.

    Input:
    ------
    z_ic - Redshift of the initial conditions
    z_cat - Redshift of the catalogue
    H0 - Hubble constant at z=0
    omch2 - Omega cold dark matter multiplied by h**2
    ombh2 - Omega baryon multiplied by h**2
    ns - The initial power spectrum slope
    sigma8 - RMS fluctuations in a sphere of 8 Mpc/h

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
    sigma8_z_ic = np.array(results.get_sigma8())[z_ic_index]
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


def calc_vel_scaling_v1(z, H0, omch2, ombh2, ns, sigma8, low_redshift=False):
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

# Change Log
# v0.1.0, 2021-10-25 - Code started with function from calc_lin_vel.py
