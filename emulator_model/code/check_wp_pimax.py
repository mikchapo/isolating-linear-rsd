# Test if wp calculation is working correclty, and what effect pimax has
# v0.1.0, 2022-05-27 - Code started

# Imports
from Corrfunc.theory import DDrppi, wp
from Corrfunc.utils import convert_rp_pi_counts_to_wp
import matplotlib.pyplot as plt
import numpy as np

from mock_meas import gen_red_cat


def calc_wp(cat, pimax=80., boxsize=1100., N_bins=9):
    N = cat.shape[0]

    nthreads = 1
    bins = np.logspace(np.log10(0.1), np.log10(60.0), N_bins+1)

    wp_results = wp(boxsize, pimax, nthreads, bins, cat[:, 4],
                    cat[:, 5], cat[:, 6])

    wp_meas = np.empty(len(wp_results))
    for i in range(len(wp_results)):
        wp_meas[i] = wp_results[i][3]

    return wp_meas


simname = "AbacusCosmos_1100box"
boxid = "00"
halotype = "rockstar"
N_grid = "1100"
smooth_type = "tophat"
path_root = ("/home/mj3chapm/scratch/abacus/{}_products/{}_{}_products/"
             "{}_{}_{}_halos/z0.700" .format(simname, simname, boxid, simname,
                                             boxid, halotype))

real_cat_paths = []
real_cat_paths.append("{}/halo_lin_vel.dat".format(path_root))
real_cat_paths.append("{}/halo_lin_vel_{}_{}_{}_smoothed"
                      ".dat".format(path_root, N_grid, smooth_type, 5.0))
real_cat_paths.append("{}/halo_lin_vel_{}_{}_{}_smoothed"
                      ".dat".format(path_root, N_grid, smooth_type, 7.0))

output_root = "test/test"
zs_real = 0.7
boxsize = 1100.
save_output = False

# pimaxs = [80., 160., boxsize/4, boxsize/2]
pimaxs = [80., 160., boxsize/4]

bin_edges = np.logspace(np.log10(0.1), np.log10(60.), 10)
log_dsep = np.log10(bin_edges[1]) - np.log10(bin_edges[0])
seps = np.power(10., (np.log10(bin_edges[:9]) + log_dsep / 2.))

linestyles = ["-", "--", "-.", ":"]
labels = ["Unsmoothed", "TH R=5", "TH R=7"]

# wp_reals = []
# plt.figure(0, figsize=(4., 3.), dpi=100)

# for i, real_cat_path in enumerate(real_cat_paths):
#     real_cat, red_cat = gen_red_cat(real_cat_path, output_root, zs_real,
#                                     boxsize, save_output,
#                                     scaling_method="split", gamma_l=1.16)

#     for j, pimax in enumerate(pimaxs):
#         wp_real = calc_wp(real_cat, pimax=pimax)
#         wp_red = calc_wp(real_cat, pimax=pimax)

#         plt.figure(0)
#         if j == 0:
#             plt.plot(seps, (wp_red - wp_real)/wp_real, color="C{}".format(i),
#                      linestyle=linestyles[j], label=labels[i])
#         else:
#             plt.plot(seps, (wp_red - wp_real)/wp_real, color="C{}".format(i),
#                      linestyle=linestyles[j])

#         if i == 0:
#             wp_reals.append(wp_real)

# plt.axhline(y=0., color='k', linestyle='-')
# plt.xscale("log")
# plt.ylabel(r"$(w_p^s - w_p^r)/w_p^r$")
# plt.xlabel(r"$r_p[h^{-1}Mpc]$")
# plt.legend()
# plt.savefig("../output/plots/check_wp_pimax_red.jpg")

real_cat_path = ("{}/halo_lin_vel_{}_{}_{}_smoothed_v2"
                 ".dat".format(path_root, N_grid, smooth_type, 4.0))
real_cat = np.loadtxt(real_cat_path)
N = real_cat.shape[0]

wp_reals = []
for j, pimax in enumerate(pimaxs):
    wp_real = calc_wp(real_cat, pimax=pimax)
    wp_reals.append(wp_real)


# Generate randoms on the box
rand_N = 3*N
rand_X = np.random.uniform(0, boxsize, rand_N)
rand_Y = np.random.uniform(0, boxsize, rand_N)
rand_Z = np.random.uniform(0, boxsize, rand_N)
nthreads = 1
pimax = 80.0

# Auto pair counts in DD
autocorr = 1
DD_counts = DDrppi(autocorr, nthreads, pimax, bin_edges, real_cat[:, 4],
                   real_cat[:, 5], real_cat[:, 6],
                   periodic=True, verbose=True)

np.savetxt("../output/data_products/DD.dat", DD_counts)

# Cross pair counts in DR
autocorr = 0
DR_counts = DDrppi(autocorr, nthreads, pimax, bin_edges, real_cat[:, 4],
                   real_cat[:, 5], real_cat[:, 6],
                   X2=rand_X, Y2=rand_Y, Z2=rand_Z,
                   periodic=True, verbose=True)

np.savetxt("../output/data_products/DR.dat", DR_counts)

# Auto pairs counts in RR
autocorr = 1
RR_counts = DDrppi(autocorr, nthreads, pimax, bin_edges, rand_X, rand_Y,
                   rand_Z, periodic=True, verbose=True)

np.savetxt("../output/data_products/RR.dat", RR_counts)

# All the pair counts are done, get the correlation function
wp_pc = convert_rp_pi_counts_to_wp(N, N, rand_N, rand_N,
                                   DD_counts, DR_counts,
                                   DR_counts, RR_counts, 9, pimax)

np.savetxt("../output/data_products/wp_pc.dat", wp_pc)

plt.figure(1, figsize=(4., 3.), dpi=100)
for i in range(1, len(wp_reals)):
    plt.plot(seps, (wp_reals[i] - wp_reals[0])/wp_reals[0], color="C0",
             linestyle=linestyles[i], label=r"$\pi_{max}=$" + str(pimaxs[i]))
plt.plot(seps, (wp_pc - wp_reals[0])/wp_reals[0], color="C1",
         linestyle=":", label=r"PC $\pi_{max}=80$")
plt.axhline(y=0., color='k', linestyle='-')
plt.xscale("log")
plt.ylabel(r"$(w_p^r - w_p^{r,80})/w_p^{r,80}$")
plt.xlabel(r"$r_p[h^{-1}Mpc]$")
plt.legend()
plt.savefig("../output/plots/check_wp_pimax_real_pc.jpg")
