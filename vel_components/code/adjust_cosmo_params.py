import numpy as np

for i in range(5):
    file_path = ("emulator_1100box_planck_0{}_products/"
                 "emulator_1100box_planck_0{}_rockstar_halos/info/"
                 "cosmo_params.dat".format(i, i))
    cosmo_params = np.loadtxt(file_path)
    cosmo_params[1] = (cosmo_params[1] * (cosmo_params[0]/100)**2 -
                       cosmo_params[2])
    np.savetxt(file_path, cosmo_params, fmt=["%.10f"])

for i in range(16):
    file_path = ("emulator_1100box_planck_00-{}_products/"
                 "emulator_1100box_planck_00-{}_rockstar_halos/info/"
                 "cosmo_params.dat".format(i, i))
    cosmo_params = np.loadtxt(file_path)
    cosmo_params[1] = (cosmo_params[1] * (cosmo_params[0]/100)**2 -
                       cosmo_params[2])
    np.savetxt(file_path, cosmo_params, fmt=["%.10f"])
