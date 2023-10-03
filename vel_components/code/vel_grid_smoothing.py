# Smooth linear velocity grid
# v0.2.0, 2022-08-31 - Added check to make sure position zeroing was correct

# Imports
import astropy.convolution as conv
import datetime as dt
import numpy as np
from scipy.stats import multivariate_normal
import sys


start_time = dt.datetime.now()
print("Starting at", start_time)


# Functions
def tophat_3d(R, l_grid):
    l_ker = int(np.ceil(R / l_grid))
    N = 2 * l_ker + 1
    win = np.zeros((N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if (np.sqrt((i - l_ker)**2 + (j - l_ker)**2 +
                            (k - l_ker)**2) <= float(R)):
                    win[i, j, k] = 1.
    # if norm:
    #     win /= 4. / 3. * np.pi * float(R)**3.

    return conv.Kernel(win)


def gaussian_3d(std, l_grid):
    l_ker = int(np.ceil(std / l_grid))
    N = 8 * l_ker + 1
    win = np.zeros((N, N, N))
    xs = np.linspace(-4 * l_ker, 4 * l_ker, N)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                win[i, j, k] = multivariate_normal.pdf([xs[i], xs[j], xs[k]],
                                                       cov=[[std**2., 0, 0],
                                                            [0, std**2., 0],
                                                            [0, 0, std**2.]])
    return conv.Kernel(win)


# Get user input parameters
boxid = sys.argv[1]
sim_name = sys.argv[2]
N_grid = int(sys.argv[3])
smooth_type = sys.argv[4]
R_smooth = float(sys.argv[5])
fft_smooth = sys.argv[6]

fft_smooth = True if fft_smooth == "True" else False

# Important global constants, that could be changed to parameters in the future
boxsize = 1100.

path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
             "{}_{}_products".format(sim_name, sim_name, boxid))

print("{}/vel_grid_{}_{}_{}_smoothed_v3".format(path_root, N_grid,
                                                smooth_type, R_smooth))

vel_grid = np.load("{}/vel_grid_{}_v3.npy".format(path_root, N_grid))

# Output to see if positions have been properly zeroed
if boxid == "00-0" and (sim_name == "AbacusCosmos_1100box_planck" and
                        N_grid == 1100 and smooth_type == "tophat" and
                        R_smooth == 5.):
    np.save("../output/data_products/test_unsmoothed_grid",
            vel_grid[806:817, 568:579, 802:813, :])

l_grid = boxsize / N_grid

if smooth_type == "tophat":
    kernel = tophat_3d(R_smooth, l_grid)

elif smooth_type == "gaussian":
    kernel = gaussian_3d(R_smooth, l_grid)

print("Starting smoothing, elapsed time", dt.datetime.now() - start_time)

if fft_smooth:
    vel_grid[:, :, :, 0] = conv.convolve_fft(vel_grid[:, :, :, 0], kernel,
                                             normalize_kernel=True,
                                             boundary="wrap", allow_huge=True)
    print("Finished vx, elapsed time", dt.datetime.now() - start_time)
    vel_grid[:, :, :, 1] = conv.convolve_fft(vel_grid[:, :, :, 1], kernel,
                                             normalize_kernel=True,
                                             boundary="wrap", allow_huge=True)
    print("Finished vy, elapsed time", dt.datetime.now() - start_time)
    vel_grid[:, :, :, 2] = conv.convolve_fft(vel_grid[:, :, :, 2], kernel,
                                             normalize_kernel=True,
                                             boundary="wrap", allow_huge=True)
    print("Finished vz, elapsed time", dt.datetime.now() - start_time)

    np.save("{}/vel_grid_{}_{}_{}_smoothed_fft_v3".format(path_root, N_grid,
                                                          smooth_type, R_smooth),
            vel_grid)

else:
    vel_grid[:, :, :, 0] = conv.convolve(vel_grid[:, :, :, 0], kernel,
                                         normalize_kernel=True,
                                         boundary="wrap")
    print("Finished vx, elapsed time", dt.datetime.now() - start_time)
    vel_grid[:, :, :, 1] = conv.convolve(vel_grid[:, :, :, 1], kernel,
                                         normalize_kernel=True,
                                         boundary="wrap")
    print("Finished vy, elapsed time", dt.datetime.now() - start_time)
    vel_grid[:, :, :, 2] = conv.convolve(vel_grid[:, :, :, 2], kernel,
                                         normalize_kernel=True,
                                         boundary="wrap")
    print("Finished vz, elapsed time", dt.datetime.now() - start_time)

    np.save("{}/vel_grid_{}_{}_{}_smoothed_v3".format(path_root, N_grid,
                                                      smooth_type, R_smooth),
            vel_grid)

# Output to see if positions have been properly zeroed
if boxid == "00-0" and (sim_name == "AbacusCosmos_1100box_planck" and
                        N_grid == 1100 and smooth_type == "tophat" and
                        R_smooth == 5.):
    # print("Test cell for proper position range (49, 1091, 1020):")
    # print("Smoothed grid values:", vel_grid[49, 1091, 1020, :])
    # print("Expected grid values: "
    #       "[160.86853854 205.63755677 273.31266465  10.        ]")
    # expected = np.array([160.86853854, 205.63755677, 273.31266465, 10.])
    # print("Difference:", vel_grid[49, 1091, 1020, :] - expected)

    print("Test cell for proper position range (811, 573, 807):")
    print("Smoothed grid values:", vel_grid[811, 573, 807, :])

print("Code complete, elapsed time", dt.datetime.now() - start_time)

# Change Log
# v0.1.0, 2022-04-21 - Code started
