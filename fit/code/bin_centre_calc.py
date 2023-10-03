# Aemulus bin centre calculation
# v0.1.0, 2022-01-11 - Copied from aemulus_fmax

# Imports
import numpy as np


# Store all the training separations
training_seps = np.empty((4000, 9))
clustering_dir = "/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Abacus/eBOSS/training_clustering/Results_Mean"
datatemp = np.loadtxt("{}/CF_cosmo_0_HOD_0_test_0_mean.dat".format(clustering_dir))
for i in range(40):
    for j in range(100):
        training_cf = np.loadtxt("{}/CF_cosmo_{}_HOD_{}_test_0_mean.dat".format(clustering_dir, i, i*100+j))
        if np.isnan(training_cf).any() == True:
            training_seps[i*100+j, :] = datatemp[:,2]
        else:
            training_seps[i*100+j, :] = training_cf[:, 2]

# Calculate mean and std of training separations
mean_seps = np.mean(training_seps, axis=0)
std_seps = np.std(training_seps, axis=0)

bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))

print("Mean separations")
print(mean_seps)
print("Sep. Std.")
print(std_seps)
print("My Seps.")
print(seps)

np.savetxt("/home/mj3chapm/P2/fit/smoothed_emulator/mean_mps_seps.dat", mean_seps)
