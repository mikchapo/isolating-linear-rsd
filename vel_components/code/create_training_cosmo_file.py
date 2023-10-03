# Create training cosmology file for the emulator
# v0.1.0, 2021-11-25 - Code started

# Import
import numpy as np

training_cosmos = np.empty((40, 7))

for i in range(40):
    if i < 10:
        cosmo = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                           "AbacusCosmos_1100box_products/"
                           "AbacusCosmos_1100box_0{}_products/"
                           "AbacusCosmos_1100box_0{}_rockstar_halos/info/"
                           "cosmo_params.dat".format(i, i))

    else:
        cosmo = np.loadtxt("/home/mj3chapm/scratch/abacus/"
                           "AbacusCosmos_1100box_products/"
                           "AbacusCosmos_1100box_{}_products/"
                           "AbacusCosmos_1100box_{}_rockstar_halos/info/"
                           "cosmo_params.dat".format(i, i))

    training_cosmos[i, 0] = (cosmo[1] + cosmo[2]) / (cosmo[0] / 100.)**2
    training_cosmos[i, 1] = cosmo[2] / (cosmo[0] / 100.)**2
    training_cosmos[i, 2] = cosmo[4]
    training_cosmos[i, 3] = cosmo[0] / 100.
    training_cosmos[i, 4] = cosmo[3]
    training_cosmos[i, 5] = 3.046
    training_cosmos[i, 6] = cosmo[6]

np.savetxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/"
           "training_cosmology.dat", training_cosmos,
           fmt=["%.7f", "%.7f", "%.7f", "%.7f", "%.7f", "%.3f", "%.7f"])
