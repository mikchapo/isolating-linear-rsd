import matplotlib.pyplot as plt
import numpy as np

from vp_plot_funcs import *


output_path = "../output/plots/lin_pred_unit_test.jpg"
cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/cosmo_params.dat")
z = 0
vp_lin_seps = np.logspace(np.log10(1.), np.log10(100.), 40)

initialize_vp_plot()

vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4], hubble_units=True, k_hunit=True)
plt.plot(vp_lin_seps, vp_lin, color="C0", linestyle="-", label="CAMB units k[h/Mpc] and P[(Mpc/h)^3]")

vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4], hubble_units=False, k_hunit=True)
plt.plot(vp_lin_seps, vp_lin, color="C1", linestyle="--", label="CAMB units k[h/Mpc] and P[Mpc^3]")

vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4], hubble_units=True, k_hunit=False)
plt.plot(vp_lin_seps, vp_lin, color="C2", linestyle="-.", label="CAMB units k[1/Mpc] and P[(Mpc/h)^3]")

vp_lin = vp_lin_pred(vp_lin_seps, z, cosmo_params[0], cosmo_params[1], cosmo_params[2], cosmo_params[3], cosmo_params[4], hubble_units=False, k_hunit=False)
plt.plot(vp_lin_seps, vp_lin, color="C3", linestyle=":", label="CAMB units k[1/Mpc] and P[Mpc^3]")

finish_vp_plot(output_path)