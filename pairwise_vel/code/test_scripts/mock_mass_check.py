import matplotlib.pyplot as plt
import numpy as np


cat_names = ["hod_z0p70_m200c1e12_test113_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test371_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test605_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test640_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test890_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test142_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test389_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test607_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test666_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test892_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test210_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test406_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test623_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test703_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test900_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test233_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test41_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test627_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test764_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test910_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test274_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test498_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test630_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test825_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test954_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test326_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test579_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test635_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test873_n1e-4_zerr.dat",
             "hod_z0p70_m200c1e12_test986_n1e-4_zerr.dat"]

# cat_names = ["hod_z0p70_m200c1e12_test113_n1e-4_zerr.dat",
#              "hod_z0p70_m200c1e12_test371_n1e-4_zerr.dat",
#              "hod_z0p70_m200c1e12_test605_n1e-4_zerr.dat",
#              "hod_z0p70_m200c1e12_test640_n1e-4_zerr.dat",
#              "hod_z0p70_m200c1e12_test890_n1e-4_zerr.dat"]

plt.figure(figsize=(8., 6.))
for i, cat_name in enumerate(cat_names):
    print("Starting Mock {}".format(i))
    cat = np.loadtxt(cat_name)
    if i==2:
        plt.hist(cat[:, 1], bins=np.logspace(10.75, 15., 100), histtype="step", color="k")
    else:
        plt.hist(cat[:, 1], bins=np.logspace(10.75, 15., 100), histtype="step")
plt.xscale("log")
plt.savefig("mock_mass_check.jpg")