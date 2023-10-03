import astropy.convolution as conv
import numpy as np
from scipy.stats import multivariate_normal
import sys


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


l_grid = 1.

for R in range(7):
    kernel = tophat_3d(R, l_grid)
    print("R={}".format(R))
    print(kernel.array[R, :, :])
