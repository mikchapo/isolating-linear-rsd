import numpy as np
import george
from george.kernels import *
import scipy.optimize as op
from scipy.linalg import cholesky, cho_solve
import sys
sys.path.insert(0, '/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASS_emulators/tools')
# from gp_training import *
# import gp_training

data = 'RSD_multiple'
x = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Abacus/eBOSS/HOD_parameters_training.dat")

x[:,0] = np.log10(x[:,0])
x[:,2] = np.log10(x[:,2])

xc = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/abacus_share/AbacusCosmos_1100box_rockstar_halos_z0.700_lin_vel/training_cosmology.dat")
xc = xc[:,[0,1,2,3,4,6]]  # remove the column for Neff  

NewKernel = False
HODLR = False

Nsize1 = 0
Nsize2 = 80

N_hod_up = 0 # the parameter to change HODs: e.g. 0~4000 -> 1000~5000 when this number is 1000                                                                       
HH = np.array(range(0,4000))
HH  = HH.reshape(40, 100)
HH = HH + N_hod_up
HH = HH[:,Nsize1:Nsize2]
CC = range(40)
rr = np.empty((HH.shape[1]*len(CC), x.shape[1]+xc.shape[1]))
YY = np.empty((9, HH.shape[1]*len(CC)))

HHmask = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Abacus/eBOSS/GP_MeanAxis/HOD_mask_monopole.dat")
pp = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Abacus/eBOSS/GP_MeanAxis/GP_mo/"+data+"_mono_9bins_pp_Nsize_"+str(Nsize1)+"_"+str(Nsize2)+"+N_hod_up"+str(N_hod_up)+".dat")

##################   find the mean of the data  #############

s2 = 0
for CID in CC:
    HH2 = HHmask[CID]
    HH3 = HH2[np.where(HH2!=-1)]
    for HID in HH3[Nsize1:Nsize2]:
        HID = int(HID)

        d = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Abacus/eBOSS/training_clustering/Results_Mean/CF_cosmo_"+str(CID)+"_HOD_"+str(HID)+"_test_0_mean.dat")
        YY[:,s2] = d[:,3]
        s2 = s2+1

Ymean = np.mean(YY, axis=1)

##################  found the mean of the data ################

GP_error = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Aemulus/CMASS_0.7/training_clustering/MedianError_from_test.dat")

GP_err = GP_error[:,3]

y2 = np.empty((len(rr)*9))
ss2 = 0
for j in range(9):
    DC = j
    Ym = Ymean[DC]
    ss = 0
    yerr = np.zeros((len(rr)))
    y = np.empty((len(rr)))
    for CID in CC:
        HH2 = HHmask[CID]
        HH3 = HH2[np.where(HH2!=-1)]
        for HID in HH3[Nsize1:Nsize2]:
            HID = int(HID)
            
            d = np.loadtxt("/home/mj3chapm/projects/rrg-wperciva/zhai/Emulator/Abacus/eBOSS/training_clustering/Results_Mean/CF_cosmo_"+str(CID)+"_HOD_"+str(HID)+"_test_0_mean.dat")

            rr[ss,0:6]=xc[CID, :]
            rr[ss,6:16]=x[HID, :]
            
            d = d[:,3]
            d1 = d[DC]
            y[ss] = np.log10(d1/Ym)
            yerr[ss] = GP_err[j]/2.303
            y2[ss2] = y[ss]
            ss = ss+1
            ss2 = ss2+1

######
    p0 = np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])

    k1 = ExpSquaredKernel(p0, ndim=len(p0))
    k2 = Matern32Kernel(p0, ndim=len(p0))
    k3 = ConstantKernel(0.1, ndim=len(p0))
    k4 = WhiteKernel(0.1, ndim=len(p0))
    k5 = ConstantKernel(0.1, ndim=len(p0))

    if NewKernel == False:
        kernel = k1*k5+k2+k3+k4
    else:
        kernel = k2+k5

    ppt = pp[j]

    if j == 0:
        if HODLR == True:
            gp0 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp0 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp0.compute(rr, yerr)

        gp0.kernel.vector = ppt
        gp0.compute(rr, yerr)

    if j == 1:
        if HODLR == True:
            gp1 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp1 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp1.compute(rr, yerr)
        
        gp1.kernel.vector = ppt
        gp1.compute(rr, yerr)

    if j == 2:
        if HODLR == True:
            gp2 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp2 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp2.compute(rr, yerr)
        
        gp2.kernel.vector = ppt
        gp2.compute(rr, yerr)

    if j == 3:
        if HODLR == True:
            gp3 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp3 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp3.compute(rr, yerr)

        gp3.kernel.vector = ppt
        gp3.compute(rr, yerr)

    if j == 4:
        if HODLR == True:
            gp4 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp4 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp4.compute(rr, yerr)
        
        gp4.kernel.vector = ppt
        gp4.compute(rr, yerr)

    if j == 5:
        if HODLR == True:
            gp5 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp5 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp5.compute(rr, yerr)
        
        gp5.kernel.vector = ppt
        gp5.compute(rr, yerr)

    if j == 6:
        if HODLR == True:
            gp6 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp6 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp6.compute(rr, yerr)
        
        gp6.kernel.vector = ppt
        gp6.compute(rr, yerr)

    if j == 7:
        if HODLR == True:
            gp7 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp7 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp7.compute(rr, yerr)
        
        gp7.kernel.vector = ppt
        gp7.compute(rr, yerr)

    if j == 8:
        if HODLR == True:
            gp8 = george.GP(kernel, mean=np.mean(y), solver=george.HODLRSolver)
        else:
            gp8 = george.GP(kernel, mean=np.mean(y), solver=george.BasicSolver)
        gp8.compute(rr, yerr)
        
        gp8.kernel.vector = ppt
        gp8.compute(rr, yerr)


