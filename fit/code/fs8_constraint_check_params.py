"""
Calculate the parameter values needed for fsigma8 constraints check.

v0.1.0, 2023-06-21 - Code started
"""

# Imports
import camb
import datetime as dt
import numpy as np
from scipy.optimize import minimize
import sys

import fit_funcs as ff

sys.path.insert(0, "/home/mj3chapm/P2/fit/smoothed_emulator")
from tools import GET_PriorND


def calc_fs8_meas(gamma_l, omegam, sigma8, z_emulator=0.7):
    """Calculate fsigma8 using the same method as the mcmc fit."""
    fsigma8_wcdm = ff.fsigma8_approximate(z_emulator, sigma8=sigma8,
                                          omegam0=omegam)
    return gamma_l * fsigma8_wcdm


def find_omegam(target_fs8, omegam, gamma_l, sigma8, threshold=0.0005):
    """
    Find the correct value of Omega_m for a given fsigma8.

    Did not use because also need to make sure the parameters are within the
    training prior.
    """
    omegam_min = 0.1
    omegam_max = 1.0
    omegam_current = omegam

    current_fs8 = calc_fs8_meas(gamma_l, omegam_current, sigma8)
    print("Target fs8:", target_fs8, "Current fs8:", current_fs8)
    n = 0
    while abs(current_fs8 - target_fs8) > threshold and n < 20:
        if current_fs8 > target_fs8:
            omegam_max = omegam_current
        else:
            omegam_min = omegam_current

        omegam_current = (omegam_min + omegam_max) / 2.
        current_fs8 = calc_fs8_meas(gamma_l, omegam_current, sigma8)
        print("Current omegam:", omegam_current, "Current fs8:", current_fs8)
        n += 1

    print("Correct fs8 found")
    return omegam_current


def calc_fs8_camb(omegam, sigma8, gamma_l, H0_marg=68.887601,
                  ombh2_marg=0.022678243, ns_marg=0.96669858, As_default=2e-9):
    """Calculate fsigma8 using CAMB."""
    pars1 = camb.CAMBparams()
    pars1.set_cosmology(H0=H0_marg, ombh2=ombh2_marg, TCMB=2.7255,
                        omch2=(omegam * (H0_marg / 100.)**2. - ombh2_marg))
    pars1.InitPower.set_params(ns=ns_marg, As=As_default)
    pars1.set_matter_power(redshifts=[0.7, 0.0], kmax=2.0)

    pars1.NonLinear = camb.model.NonLinear_none
    results1 = camb.get_results(pars1)
    kh1, zcamb1, pk1 = results1.get_matter_power_spectrum(minkh=1e-4, maxkh=1,
                                                          npoints=200)
    sigma8_default = np.array(results1.get_sigma8())

    ratio = sigma8 / sigma8_default[1]
    As_desired = As_default * ratio**2.

    pars2 = camb.CAMBparams()
    pars2.set_cosmology(H0=H0_marg, ombh2=ombh2_marg, TCMB=2.7255,
                        omch2=(omegam * (H0_marg / 100.)**2. - ombh2_marg))
    pars2.InitPower.set_params(ns=ns_marg, As=As_desired)
    pars2.set_matter_power(redshifts=[0.7, 0.0], kmax=2.0)

    pars2.NonLinear = camb.model.NonLinear_none
    results2 = camb.get_results(pars2)
    kh2, zcamb2, pk2 = results2.get_matter_power_spectrum(minkh=1e-4, maxkh=1,
                                                          npoints=200)
    sigma8_final = np.array(results2.get_sigma8())
    fsigma8_final = np.array(results2.get_fsigma8())

    print("CAMB fsigma8 for omegam={}, sigma8={}, gamma_l="
          "{}".format(omegam, sigma8, gamma_l))
    print("Resulting sigma8:", sigma8_final[1])
    print("Resulting fsigma8:", fsigma8_final)
    print("Resulting gamma_l fsigma8:", gamma_l * fsigma8_final, "\n")


def calc_fs8_meas_tp(x, sigma8_marg=8.1526127E-01, gamma_l=7.9454733E-01,
                     z_emulator=0.7, h_marg=0.68887601, ombh2_marg=0.022678243,
                     ns_marg=0.96669858, CUT=12):
    """Tried to do this using scipy.optimize.minimized, but did not work."""
    GG = GET_PriorND()
    omegam = x[0]
    omegab_marg = ombh2_marg / h_marg**2.
    w = -1.
    params = [omegam, omegab_marg, sigma8_marg, h_marg, ns_marg, w]

    if GG.isinornot(params, CUT):
        fsigma8_wcdm = ff.fsigma8_approximate(z_emulator, sigma8=sigma8_marg,
                                              omegam0=omegam)
        return gamma_l * fsigma8_wcdm

    else:
        return 10.


def find_tp_extreme_fixed(find_om=True, find_s8=True, ext_type="min",
                          omegam_marg=2.9657580E-01, sigma8_marg=8.1526127E-01,
                          h_marg=0.68887601, ombh2_marg=0.022678243,
                          ns_marg=0.96669858, CUT=12, increment=0.001):
    """Didn't work because the other parameters affect the training prior."""
    GG = GET_PriorND()
    omegab_marg = ombh2_marg / h_marg**2.
    w = -1.
    params = [omegam_marg, omegab_marg, sigma8_marg, h_marg, ns_marg, w]

    if ext_type == "min":
        increment *= -1

    if find_om and find_s8:
        print("Good luck!")
        return 0

    elif find_om or find_s8:
        if find_om:
            index = 0

        else:
            index = 2

        searching = True
        while searching:
            if GG.isinornot(params, CUT):
                params[index] += increment

            else:
                searching = False
                return (params[index] - increment)

    else:
        print("Well you have to find something!")
        return 0


def find_tp_extreme(N_iter, fixed_omegam=None, fixed_sigma8=None, CUT=12,
                    seed=42):
    start_time = dt.datetime.now()
    np.random.seed(seed)
    GG = GET_PriorND()
    prior_ranges = [{"label": "omegam", "max": 0.3668427, "min": 0.2529504},
                    {"label": "omegab", "max": 0.0601736, "min": 0.0393986},
                    {"label": "sigma8", "max": 0.9986873, "min": 0.6466325},
                    {"label": "h", "max": 0.7479261, "min": 0.6156746},
                    {"label": "ns", "max": 0.9897832, "min": 0.9300325}]

    min_info = {"fsigma8": 0.466, "omegam": 0.297, "ns": 0.967,
                "omegab": 0.0227 / 0.689**2., "sigma8": 0.815, "h": 0.689}
    max_info = {"fsigma8": 0.466, "omegam": 0.297, "ns": 0.967,
                "omegab": 0.0227 / 0.689**2., "sigma8": 0.815, "h": 0.689}

    params = [min_info["omegam"], min_info["omegab"],
              min_info["sigma8"], min_info["h"], min_info["ns"], -1.]

    if fixed_omegam is not None and fixed_sigma8 is not None:
        print("You need to vary either omegam or sigma8")
        return False

    elif fixed_omegam is not None:
        params[0] = fixed_omegam
        param_indices = range(1, 5)

    elif fixed_sigma8 is not None:
        params[2] = fixed_sigma8
        param_indices = [0, 1, 3, 4]

    else:
        param_indices = range(5)

    for i in range(N_iter):
        for j in param_indices:
            params[j] = np.random.uniform(low=prior_ranges[j]["min"],
                                          high=prior_ranges[j]["max"])
        # params = [2.9657580E-01, 0.022678243, 8.1526127E-01, 0.68887601,
        #           0.96669858, -1.]
        # print("i={}, params are".format(i), params)
        # params[1] = params[1] / params[3]**2.
        if GG.isinornot(params, CUT):
            fsigma8 = calc_fs8_meas(1., params[0], params[2])
            # print("Within training prior, i={}, fsigma8={}, params are".format(i, fsigma8), params)
            if fsigma8 < min_info["fsigma8"]:
                min_info = {"fsigma8": fsigma8, "omegam": params[0],
                            "omegab": params[1],
                            "sigma8": params[2], "h": params[3],
                            "ns": params[4]}

            elif fsigma8 > max_info["fsigma8"]:
                max_info = {"fsigma8": fsigma8, "omegam": params[0],
                            "omegab": params[1],
                            "sigma8": params[2], "h": params[3],
                            "ns": params[4]}

        if (i + 1) % (N_iter / 10) == 0:
            print("Completed {} iterations, elapsed time "
                  "{}".format(i+1, dt.datetime.now() - start_time))

    return min_info, max_info


# Marginalized constraints
gamma_l_marg = 7.9454733E-01
omegam_marg = 2.9657580E-01
sigma8_marg = 8.1526127E-01
H0_marg = 68.887601
ombh2_marg = 0.022678243
ns_marg = 0.96669858

# omegam_min = 0.2529504
# omegam_max = 0.3668427

fs8_marg = 3.6753341E-01
fs8_std_marg = 4.0653156E-02
fs8_calc_marg = calc_fs8_meas(gamma_l_marg, omegam_marg, sigma8_marg)
print()
print("Marginalized fsigma8: {}, fsigma8 calculated from marginalized "
      "parameters: {}".format(fs8_marg, fs8_calc_marg), "\n")


"""
# Find fixed values of omegam and sigma8 within the training prior

min_info_100000, max_info_100000 = find_tp_extreme(100000)
print("100000 iterations min fsigma8={}; Om={}, Ob={}, s8={}, h={}, ns="
      "{}".format(min_info_100000["fsigma8"], min_info_100000["omegam"],
                  min_info_100000["omegab"], min_info_100000["sigma8"],
                  min_info_100000["h"], min_info_100000["ns"]))
print("100000 iterations max fsigma8={}; Om={}, Ob={}, s8={}, h={}, ns="
      "{}".format(max_info_100000["fsigma8"], max_info_100000["omegam"],
                  max_info_100000["omegab"], max_info_100000["sigma8"],
                  max_info_100000["h"], max_info_100000["ns"]), "\n")

min_info_1000000, max_info_1000000 = find_tp_extreme(1000000)
print("1000000 iterations min fsigma8={}; Om={}, Ob={}, s8={}, h={}, ns="
      "{}".format(min_info_1000000["fsigma8"], min_info_1000000["omegam"],
                  min_info_1000000["omegab"], min_info_1000000["sigma8"],
                  min_info_1000000["h"], min_info_1000000["ns"]))
print("1000000 iterations max fsigma8={}; Om={}, Ob={}, s8={}, h={}, ns="
      "{}".format(max_info_1000000["fsigma8"], max_info_1000000["omegam"],
                  max_info_1000000["omegab"], max_info_1000000["sigma8"],
                  max_info_1000000["h"], max_info_1000000["ns"]), "\n")

(min_info_om_100000,
 max_info_om_100000) = find_tp_extreme(100000, fixed_omegam=omegam_marg)
print("100000 iterations with fixed omegam, min fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(min_info_om_100000["fsigma8"],
                                  min_info_om_100000["omegam"],
                                  min_info_om_100000["omegab"],
                                  min_info_om_100000["sigma8"],
                                  min_info_om_100000["h"],
                                  min_info_om_100000["ns"]))
print("100000 iterations with fixed omegam, max fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(max_info_om_100000["fsigma8"],
                                  max_info_om_100000["omegam"],
                                  max_info_om_100000["omegab"],
                                  max_info_om_100000["sigma8"],
                                  max_info_om_100000["h"],
                                  max_info_om_100000["ns"]), "\n")

(min_info_om_1000000,
 max_info_om_1000000) = find_tp_extreme(1000000, fixed_omegam=omegam_marg)
print("1000000 iterations with fixed omegam, min fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(min_info_om_1000000["fsigma8"],
                                  min_info_om_1000000["omegam"],
                                  min_info_om_1000000["omegab"],
                                  min_info_om_1000000["sigma8"],
                                  min_info_om_1000000["h"],
                                  min_info_om_1000000["ns"]))
print("1000000 iterations with fixed omegam, max fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(max_info_om_1000000["fsigma8"],
                                  max_info_om_1000000["omegam"],
                                  max_info_om_1000000["omegab"],
                                  max_info_om_1000000["sigma8"],
                                  max_info_om_1000000["h"],
                                  max_info_om_1000000["ns"]), "\n")

(min_info_s8_100000,
 max_info_s8_100000) = find_tp_extreme(100000, fixed_sigma8=sigma8_marg)
print("100000 iterations with fixed sigma8, min fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(min_info_s8_100000["fsigma8"],
                                  min_info_s8_100000["omegam"],
                                  min_info_s8_100000["omegab"],
                                  min_info_s8_100000["sigma8"],
                                  min_info_s8_100000["h"],
                                  min_info_s8_100000["ns"]))
print("100000 iterations with fixed sigma8, max fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(max_info_s8_100000["fsigma8"],
                                  max_info_s8_100000["omegam"],
                                  max_info_s8_100000["omegab"],
                                  max_info_s8_100000["sigma8"],
                                  max_info_s8_100000["h"],
                                  max_info_s8_100000["ns"]))

(min_info_s8_1000000,
 max_info_s8_1000000) = find_tp_extreme(1000000, fixed_sigma8=sigma8_marg)
print("1000000 iterations with fixed sigma8, min fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(min_info_s8_1000000["fsigma8"],
                                  min_info_s8_1000000["omegam"],
                                  min_info_s8_1000000["omegab"],
                                  min_info_s8_1000000["sigma8"],
                                  min_info_s8_1000000["h"],
                                  min_info_s8_1000000["ns"]))
print("1000000 iterations with fixed sigma8, max fsigma8={}; Om={}, Ob={}, "
      "s8={}, h={}, ns={}".format(max_info_s8_1000000["fsigma8"],
                                  max_info_s8_1000000["omegam"],
                                  max_info_s8_1000000["omegab"],
                                  max_info_s8_1000000["sigma8"],
                                  max_info_s8_1000000["h"],
                                  max_info_s8_1000000["ns"]))
"""


# Calculate parameters for varying

omegam_min = 0.2763982054070388
omegam_max = 0.3293361232823127
sigma8_min = 0.736129870588657
sigma8_max = 0.8989857305264024
fs8_both_min = calc_fs8_meas(gamma_l_marg, omegam_min, sigma8_min)
fs8_both_max = calc_fs8_meas(gamma_l_marg, omegam_max, sigma8_max)

omegam_min_fixed_s8 = 0.25299337446824643
# omegam_max_fixed_s8 = 0.36226835384298794
omegam_max_fixed_s8 = 0.36065646710102783 # 0.9505458812459266
fs8_omegam_min = calc_fs8_meas(gamma_l_marg, omegam_min_fixed_s8, sigma8_marg)
fs8_omegam_max = calc_fs8_meas(gamma_l_marg, omegam_max_fixed_s8, sigma8_marg)

sigma8_min_fixed_om = 0.7400843220242559
sigma8_max_fixed_om = 0.893933143111853
fs8_sigma8_min = calc_fs8_meas(gamma_l_marg, omegam_marg, sigma8_min_fixed_om)
fs8_sigma8_max = calc_fs8_meas(gamma_l_marg, omegam_marg, sigma8_max_fixed_om)

# print("Omegam min: {}, Omegam max: {}".format(omegam_min, omegam_max))
# print("Sigma8 min: {}, Sigma8 max: {}".format(sigma8_min, sigma8_max),
#       "\n\n")

# Varied parameters
for std_var in range(1, 4):
    fs8_low = fs8_marg - std_var * fs8_std_marg
    fs8_high = fs8_marg + std_var * fs8_std_marg
    var_ratio_low = fs8_low / fs8_marg
    var_ratio_high = fs8_high / fs8_marg
    print("Std var: {}, fs8 low: {}, fs8 high: "
          "{}".format(std_var, fs8_low, fs8_high))

    gamma_l_low = var_ratio_low * gamma_l_marg
    gamma_l_high = var_ratio_high * gamma_l_marg
    print("gamma low: {}, gamma high: "
          "{}".format(gamma_l_low, gamma_l_high), "\n")

    var_ratio_split_both_low = ((fs8_marg - std_var * fs8_std_marg) /
                                fs8_both_min)
    gamma_l_split_both_low = var_ratio_split_both_low * gamma_l_marg
    fs8_split_both_low = calc_fs8_meas(gamma_l_split_both_low, omegam_min,
                                       sigma8_min)
    print("Split both low, gamma_l={}, omegam={}, sigma8={}, fs8="
          "{}".format(gamma_l_split_both_low, omegam_min, sigma8_min,
                      fs8_split_both_low))

    var_ratio_split_both_high = ((fs8_marg + std_var * fs8_std_marg) /
                                 fs8_both_max)
    gamma_l_split_both_high = var_ratio_split_both_high * gamma_l_marg
    fs8_split_both_high = calc_fs8_meas(gamma_l_split_both_high,
                                        omegam_max, sigma8_max)
    print("Split both high, gamma_l={}, omegam={}, sigma8={}, fs8="
          "{}".format(gamma_l_split_both_high, omegam_max, sigma8_max,
                      fs8_split_both_high), "\n")

    var_ratio_split_omegam_low = ((fs8_marg - std_var * fs8_std_marg) /
                                  fs8_omegam_min)
    gamma_l_split_omegam_low = var_ratio_split_omegam_low * gamma_l_marg
    fs8_split_omegam_low = calc_fs8_meas(gamma_l_split_omegam_low,
                                         omegam_min_fixed_s8, sigma8_marg)
    print("Split Omegam low, gamma_l={}, omegam={}, sigma8={}, fs8="
          "{}".format(gamma_l_split_omegam_low, omegam_min_fixed_s8,
                      sigma8_marg, fs8_split_omegam_low))

    var_ratio_split_omegam_high = ((fs8_marg + std_var * fs8_std_marg) /
                                   fs8_omegam_max)
    gamma_l_split_omegam_high = var_ratio_split_omegam_high * gamma_l_marg
    fs8_split_omegam_high = calc_fs8_meas(gamma_l_split_omegam_high,
                                          omegam_max_fixed_s8, sigma8_marg)
    print("Split Omegam high, gamma_l={}, omegam={}, sigma8={}, fs8="
          "{}".format(gamma_l_split_omegam_high, omegam_max_fixed_s8,
                      sigma8_marg, fs8_split_omegam_high), "\n")

    var_ratio_split_sigma8_low = ((fs8_marg - std_var * fs8_std_marg) /
                                  fs8_sigma8_min)
    gamma_l_split_sigma8_low = var_ratio_split_sigma8_low * gamma_l_marg
    fs8_split_sigma8_low = calc_fs8_meas(gamma_l_split_sigma8_low, omegam_marg,
                                         sigma8_min_fixed_om)
    print("Split Sigma8 low, gamma_l={}, omegam={}, sigma8={}, fs8="
          "{}".format(gamma_l_split_sigma8_low, omegam_marg,
                      sigma8_min_fixed_om, fs8_split_sigma8_low))

    var_ratio_split_sigma8_high = ((fs8_marg + std_var * fs8_std_marg) /
                                   fs8_sigma8_max)
    gamma_l_split_sigma8_high = var_ratio_split_sigma8_high * gamma_l_marg
    fs8_split_sigma8_high = calc_fs8_meas(gamma_l_split_sigma8_high,
                                          omegam_marg, sigma8_max_fixed_om)
    print("Split Sigma8 high, gamma_l={}, omegam={}, sigma8={}, fs8="
          "{}".format(gamma_l_split_sigma8_high, omegam_marg,
                      sigma8_max_fixed_om, fs8_split_sigma8_high), "\n\n")


# calc_fs8_camb(omegam_marg, sigma8_marg, gamma_l_marg)
# calc_fs8_camb(omegam_min, sigma8_marg, gamma_l_marg)
# calc_fs8_camb(omegam_max, sigma8_marg, gamma_l_marg)

# fs8_om_min_fit = calc_fs8_meas(gamma_l_marg, omegam_min, sigma8_marg)
# fs8_om_max_fit = calc_fs8_meas(gamma_l_marg, omegam_max, sigma8_marg)

# print("fs8 om min fit: {}, fs8 om max fit: "
#       "{}".format(fs8_om_min_fit, fs8_om_max_fit), "\n")
