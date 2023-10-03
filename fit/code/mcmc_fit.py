# Perform MCMC fit using the Abacus emulator
# v1.1.0, 2022-01-10 - Updated for smoothed linear velocity emulator, including
#                      moving scripts to the code directory and adding
#                      "emulator" as an input

# Imports
from cobaya.run import run
import numpy as np
import sys
import yaml

import fit_funcs as ff

print("Import Complete")

c = 299792.458
Omega_M_fid = 0.31
z_emulator = 0.7

# Load Params
input_file = sys.argv[1]
with open(input_file) as file:
    input_dict = yaml.load(file, Loader=yaml.FullLoader)

sys.path.insert(0, "/home/mj3chapm/P2/fit/{}".format(input_dict["emulator"]))
# Importing the emulator takes 30-120 seconds
import GP_prediction as emulator
from tools import GET_PriorND

zz_seps = np.loadtxt("/home/mj3chapm/P2/fit/{}/mean_mps_seps"
                     ".dat".format(input_dict["emulator"]))

print("Output Path:", input_dict["chain_dir"])

indices = ff.find_indices(input_dict)

# v0.5 will need to add input_dict["z_data"] as z arg
D_M_fid = ff.D_M(Omega_M_fid, input_dict["z_data"])
D_H_fid = ff.D_H(Omega_M_fid, input_dict["z_data"])

data_corr_funcs = np.loadtxt(input_dict["data_path"])
DATA_VECTOR = data_corr_funcs[indices, -1]
cov_matrix = np.loadtxt(input_dict["cov_mat_path"])
n_cm = 200
INV_COV_MATRIX = ff.invert_cov_matrix(cov_matrix[indices][:, indices], n_cm)

print("Data and Cov Matrix Generated")

GG = GET_PriorND()
CUT = 12


def jpa_emulator_likelihood(ombh2, ns, w, logM_sat, alpha, logM_cut, sigma_logM, v_bc, v_bs, c_vir, gamma_n, fmax, gamma_l, sigma8_free, _self=None):
    '''
    Calculates log probability of model fit to data using inverse covariance matrix, and emulator correlation function predictions.
    '''    
    h = _self.provider.get_param("H0") / 100.
    omegab = ombh2 / h / h
    omegam = _self.provider.get_param("omegam")
    sigma8 = _self.provider.get_param("sigma8")
    fsigma8_wcdm = _self.provider.get_fsigma8(input_dict["z_data"])[0]

    fsigma8_meas = gamma_l*fsigma8_wcdm

    if input_dict["s8_treatment"] == "scale":
        sigma8_emulator = sigma8 * ff.calc_D(input_dict["z_data"], omegam, 1. - omegam) / ff.calc_D(z_emulator, omegam, 1. - omegam)
        fsig8 = ff.fsigma8_approximate(z_emulator, sigma8=sigma8, omegam0=omegam) 
        fsigma8_comp = gamma_l * fsig8 / sigma8 * sigma8_emulator

    elif input_dict["s8_treatment"] == "no_scale":
        sigma8_emulator = sigma8
        fsig8 = ff.fsigma8_approximate(z_emulator, sigma8=sigma8, omegam0=omegam) 
        fsigma8_comp = gamma_l * fsig8

    elif input_dict["s8_treatment"] == "free":
        sigma8_emulator = sigma8_free
        fsig8 = ff.fsigma8_approximate(z_emulator, sigma8=sigma8_free, omegam0=omegam)
        fsigma8_comp = gamma_l * fsig8

    else:
        raise ValueError("Unrecognized s8 treatment type")

    derived = {'sigma8_aem': sigma8_emulator, 'fsigma8_wcdm': fsigma8_wcdm, 'fsigma8_comp': fsigma8_comp, 'fsigma8_meas': fsigma8_meas}

    params = [omegam, omegab, sigma8_emulator, h, ns, w, logM_sat, alpha, logM_cut, sigma_logM, v_bc, v_bs, c_vir, gamma_n, fmax, gamma_l]

    if GG.isinornot(params[0:6], CUT) or not input_dict["training_prior"]:
        mono_pred = emulator.mono_pre(params)
        quad_pred = emulator.quad_pre(params) / (zz_seps**2.)
        wp_pred = emulator.wp_pre(params)

        if input_dict["include_ap"]:
            D_M_camb = _self.provider.get_angular_diameter_distance(input_dict["z_data"]) * (1. + input_dict["z_data"]) * h
            D_H_camb = 1. / _self.provider.get_Hubble(input_dict["z_data"]) * h

            alpha_perp = D_M_camb / D_M_fid
            alpha_para = D_H_camb / D_H_fid

            mono_pred, quad_pred, wp_pred = ff.ap_scaling(mono_pred, quad_pred, wp_pred, alpha_perp, alpha_para)

        model_vector = np.concatenate([mono_pred, quad_pred, wp_pred])[indices]

        return -0.5 * np.asmatrix(DATA_VECTOR - model_vector) * np.asmatrix(INV_COV_MATRIX) * np.asmatrix(DATA_VECTOR - model_vector).T, derived

    else:
        return -np.inf, derived


def ao_emulator_likelihood(omegam, ombh2, sigma8, H0, ns, w, logM_sat, alpha, logM_cut, sigma_logM, v_bc, v_bs, c_vir, gamma_n, fmax, gamma_l, sigma8_free, _self=None):
    '''
    Calculates log probability of model fit to data using inverse covariance matrix, and emulator correlation function predictions.
    '''    
    h = H0 / 100.
    omegab = ombh2 / h / h
    fsigma8_wcdm = ff.fsigma8_approximate(z_emulator, sigma8=sigma8, omegam0=omegam) 

    fsigma8_meas = gamma_l * fsigma8_wcdm

    derived = {'sigma8_aem': sigma8, 'fsigma8_wcdm': fsigma8_wcdm, 'fsigma8_comp': fsigma8_meas, 'fsigma8_meas': fsigma8_meas}

    if input_dict["equal_gammas"]:
        gamma_n = gamma_l

    params = [omegam, omegab, sigma8, h, ns, w, logM_sat, alpha, logM_cut, sigma_logM, v_bc, v_bs, c_vir, gamma_n, fmax, gamma_l]

    if GG.isinornot(params[0:6], CUT) or not input_dict["training_prior"]:
        mono_pred = emulator.mono_pre(params)
        quad_pred = emulator.quad_pre(params) / (zz_seps**2.)
        wp_pred = emulator.wp_pre(params)

        if input_dict["include_ap"]:
            D_M_camb = ff.D_M(omegam, input_dict["z_data"])
            D_H_camb = ff.D_H(omegam, input_dict["z_data"])

            alpha_perp = D_M_camb / D_M_fid
            alpha_para = D_H_camb / D_H_fid

            mono_pred, quad_pred, wp_pred = ff.ap_scaling(mono_pred, quad_pred, wp_pred, alpha_perp, alpha_para)

        model_vector = np.concatenate([mono_pred, quad_pred, wp_pred])[indices]

        return -0.5 * np.asmatrix(DATA_VECTOR - model_vector) * np.asmatrix(INV_COV_MATRIX) * np.asmatrix(DATA_VECTOR - model_vector).T, derived

    else:
        return -np.inf, derived

with open(input_dict["cobaya_input_file"]) as file:
    info = yaml.load(file, Loader=yaml.FullLoader)

# This code needs to be checked
info["sampler"]["mcmc"]["covmat"] = input_dict["fit_cov_mat"]
info["output"] = input_dict["chain_dir"]
info["resume"] = input_dict["resume"]
info["force"] = not input_dict["resume"]

with open(input_dict["ref_input_file"]) as file:
    ref_input = yaml.load(file, Loader=yaml.FullLoader)

with open(input_dict["prior_input_file"]) as file:
    prior_input = yaml.load(file, Loader=yaml.FullLoader)

if "joint_planck_emulator" in input_dict["cobaya_input_file"] or "jpa" in input_dict["cobaya_input_file"]:
    print("Using joint Planck-emulator likelihood")
    PARAM_LABELS = ['ombh2', 'ns', "logM_sat", "alpha", "logM_cut", "sigma_logM", "v_bc", "v_bs", "c_vir", "gamma_n", "fmax", "gamma_l"]
    info["likelihood"]["emulator_likelihood"]["external"] = jpa_emulator_likelihood
    info["likelihood"]["emulator_likelihood"]["requires"]["angular_diameter_distance"]["z"] = [input_dict["z_data"]]
    info["likelihood"]["emulator_likelihood"]["requires"]["Hubble"]["z"] = [input_dict["z_data"]]
    info["likelihood"]["emulator_likelihood"]["requires"]["fsigma8"]["z"] = [input_dict["z_data"]]
    info["theory"]["camb"]["extra_args"]["redshifts"] = [0, input_dict["z_data"]]

if "emulator_only" in input_dict["cobaya_input_file"] or "ao" in input_dict["cobaya_input_file"]:
    print("Using emulator only likelihood")
    PARAM_LABELS = ['omegam', 'ombh2', 'sigma8', 'H0', 'ns', "logM_sat", "alpha", "logM_cut", "sigma_logM", "v_bc", "v_bs", "c_vir", "gamma_n", "fmax", "gamma_l"]
    info["likelihood"]["emulator_likelihood"]["external"] = ao_emulator_likelihood

for param in PARAM_LABELS:
    fixed = False
    for key in ref_input[param]:
        info['params'][param][key] = ref_input[param][key]
        if key == "value":
            fixed = True

    if not fixed:
        info['params'][param]['prior'] = prior_input[param]


for key in ref_input['sigma8']:
    fixed = False
    info['params']['sigma8_free'][key] = ref_input['sigma8'][key]
    if key == "value":
        fixed = True

    if not fixed:
        info['params']['sigma8_free']['prior'] = prior_input['sigma8']

updated_info, products = run(info)

print("Run Complete")

# Change Log
# v1.0.0, 2022-01-10 - Code started from RSD/fit/aemulus_fmax/mcmc_fit_v2.py
