# Test FFT script and make sure everything is working right
# v0.1.0, 2022-04-07 - Code started

# Imports
import astropy.convolution as conv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np
from scipy.fft import fft, ifft, fft2, ifft2
from scipy.signal import general_gaussian
from scipy.stats import multivariate_normal, norm


# Functions


def plot_two_mats(mat0, mat1, title0, title1, output, vmin=-1, vmax=1):
    # Check that FT-IFT gives back same matrix
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.75))

    ais = axes[0].imshow(mat0, origin="lower", vmin=vmin, vmax=vmax)
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(ais, cax=cax)
    axes[0].set_title(title0)

    ais = axes[1].imshow(mat1, origin="lower", vmin=vmin, vmax=vmax)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(ais, cax=cax)
    axes[1].set_title(title1)

    plt.savefig(output)


def tophat_2d(R):
    N = int(2 * R + 1)
    win = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if np.sqrt((i - R)**2 + (j - R)**2) <= float(R):
                win[i, j] = 1.
    # win /= np.pi * float(R)**2.
    return conv.Kernel(win)


def plot_3d_slices(field, xslice, yslice, zslice, output, vmin=-1., vmax=1.):
    fig, axes = plt.subplots(1, 3, figsize=(19., 6.))

    ais = axes[0].imshow(field[xslice, :, :], origin="lower", vmin=vmin,
                         vmax=vmax)
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(ais, cax=cax)
    axes[0].set_title("Slice at x={}".format(xslice))

    ais = axes[1].imshow(field[:, yslice, :], origin="lower", vmin=vmin,
                         vmax=vmax)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(ais, cax=cax)
    axes[1].set_title("Slice at y={}".format(yslice))

    ais = axes[2].imshow(field[:, :, zslice], origin="lower", vmin=vmin,
                         vmax=vmax)
    divider = make_axes_locatable(axes[2])
    cax = divider.append_axes("right", size="5%", pad=0.1)
    plt.colorbar(ais, cax=cax)
    axes[2].set_title("Slice at z={}".format(zslice))

    plt.savefig(output)


def tophat_3d(R, norm=True):
    N = int(2 * R + 1)
    win = np.zeros((N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if np.sqrt((i - R)**2 + (j - R)**2 + (k - R)**2) <= float(R):
                    win[i, j, k] = 1.
    if norm:
        win /= 4. / 3. * np.pi * float(R)**3.
    return conv.Kernel(win)


def guassian_1d(std):
    xs = np.linspace(-4 * std, 4 * std, 8 * std + 1)
    win = norm.pdf(xs, scale=std)
    return conv.Kernel(win)


def guassian_2d(std):
    N = int(8 * std + 1)
    win = np.zeros((N, N))
    xs = np.linspace(-4 * std, 4 * std, N)
    for i in range(N):
        for j in range(N):
            win[i, j] = multivariate_normal.pdf([xs[i], xs[j]],
                                                cov=[[std**2., 0],
                                                     [0, std**2.]])
    return conv.Kernel(win)


def guassian_3d(std):
    N = 8 * std + 1
    win = np.zeros((N, N, N))
    xs = np.linspace(-4 * std, 4 * std, N)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                win[i, j, k] = multivariate_normal.pdf([xs[i], xs[j], xs[k]],
                                                       cov=[[std**2., 0, 0],
                                                            [0, std**2., 0],
                                                            [0, 0, std**2.]])
    return conv.Kernel(win)


def tophat_cont(R, l_grid):
    l_ker = int(np.ceil(R / l_grid))
    N = 2 * l_ker + 1
    win = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if (np.sqrt((i - l_ker)**2 + (j - l_ker)**2) <= float(R)):
                win[i, j] = 1.

    return conv.Kernel(win)


def guassian_cont(std, l_grid):
    l_ker = int(np.ceil(std / l_grid))
    N = 8 * l_ker + 1
    win = np.zeros((N, N))
    xs = np.linspace(-4 * l_ker, 4 * l_ker, N)
    for i in range(N):
        for j in range(N):
            win[i, j] = multivariate_normal.pdf([xs[i], xs[j]],
                                                cov=[[std**2., 0],
                                                     [0, std**2.]])
    return conv.Kernel(win)


"""
# 1D Test
N = 128
spike_array = np.zeros(N)
spike_locs = [N//6, N//5, N//5, N//4, N//3, N//2, N//2]
for loc in spike_locs:
    spike_array[loc] += 1.
# spike_array[N//4] = 1.
spike_array_ft = fft(spike_array)
spike_array_inv = ifft(spike_array_ft)

# L = N // 2
L = 2*N
win = np.roll(general_gaussian(L, p=1, sig=4), L//2)
# L = 16
# win = np.zeros(N)
# win[:L//2] = 1. / L
# win[-L//2:] = 1. / L

# plt.figure(figsize=(8., 6.))
# plt.plot(spike_array, label="Spike", color="C0")
# plt.plot(np.real(win), label="Window", color="C1")
# plt.legend()
# plt.show()

win_ft = fft(win, N)
transform_ft = win_ft * spike_array_ft
transform = ifft(transform_ft)

plt.figure(figsize=(8., 6.))
plt.plot(spike_array, label="Spike", color="C0")
plt.plot(np.real(transform), label="Transform", color="C2")
# plt.plot(np.real(spike_array_ft), label="Spike FT", color="C1")
# plt.plot(np.real(spike_array_inv), label="Spike Inv.", color="C1")
plt.legend()
plt.savefig("../output/plots/fft_test/2N_1d_gaus.jpg")
"""

"""
# 2D Test

N = 128
# spike_matrix = np.zeros((N, N))
# spike_matrix[N // 2, N // 2] = 1.
spike_matrix = (np.random.rand(N, N) - 0.5) * 2.

spike_matrix_ft = fft2(spike_matrix)
spike_matrix_inv = ifft2(spike_matrix_ft)

# plot_two_mats(spike_matrix, np.real(spike_matrix_inv), "Spike Field",
#               "Spike FT-IFT",
#               "../output/plots/fft_test/spike_ft-ift_check.jpg")

# plot_two_mats(spike_matrix, np.real(spike_matrix_ft), "Single spike",
#               "Spike FT",
#               "../output/plots/fft_test/spike_ft.jpg")

R = 4
win = np.zeros((N, N))
for i in range(-R, R + 1):
    for j in range(-R, R + 1):
        if np.sqrt(i**2 + j**2) <= float(R):
            win[i, j] = 1.
win /= np.pi * float(R)**2.
win_ft = fft2(win)
smooth_field_ft = win_ft * spike_matrix_ft
smooth_field = np.real(ifft2(smooth_field_ft))

plot_two_mats(spike_matrix, smooth_field, "Spike Field",
              "Smooth Field",
              "../output/plots/fft_test/smoothed_spike_field.jpg")

kernel = conv.Tophat2DKernel(R)
my_kernel = tophat_2d(R)
conv_field = conv.convolve(spike_matrix, kernel, boundary="wrap")
conv_fft_field = conv.convolve_fft(spike_matrix, kernel, boundary="wrap")
my_conv_field = conv.convolve(spike_matrix, my_kernel, boundary="wrap")
my_conv_fft_field = conv.convolve_fft(spike_matrix, my_kernel, boundary="wrap")

plot_two_mats(spike_matrix, conv_field, "Spike Field",
              "AP Convolved Field",
              "../output/plots/fft_test/convolved_spike_field.jpg")
plot_two_mats(spike_matrix, conv_fft_field, "Spike Field",
              "AP Convolved FFT Field",
              "../output/plots/fft_test/convolved_fft_spike_field.jpg")


diff_smooth_conv = smooth_field - conv_field
diff_smooth_fft = smooth_field - conv_fft_field

plot_two_mats(diff_smooth_conv, diff_smooth_fft, "Diff Smooth Conv",
              "Diff Smooth FFT",
              "../output/plots/fft_test/spike_field_diffs.jpg", vmin=-0.01,
              vmax=0.01)

my_diff_smooth_conv = my_conv_field - conv_field
my_diff_smooth_fft = my_conv_fft_field - conv_fft_field

plot_two_mats(my_diff_smooth_conv, my_diff_smooth_fft, "Diff Smooth Conv",
              "Diff Smooth FFT",
              "../output/plots/fft_test/my_spike_field_diffs.jpg", vmin=-0.001,
              vmax=0.001)

"""

"""
# 3D Test

# N = 16
# R = 4
# kernel = tophat_3d(R, norm=False)
# # print(kernel.array)

# spike_field = np.zeros((N, N, N))
# spike_field[N//2, N//2, N//2] = 1.
# field_name = "single_spike"

# conv_field = conv.convolve(spike_field, kernel, normalize_kernel=False,
#                            boundary="wrap")
# conv_fft_field = conv.convolve_fft(spike_field, kernel,
#                                    normalize_kernel=False,
#                                    boundary="wrap")

N = 128
R = 4
kernel = tophat_3d(R)

spike_field = (np.random.rand(N, N, N) - 0.5) * 2.
field_name = "spike_field"

conv_field = conv.convolve(spike_field, kernel, boundary="wrap")
conv_fft_field = conv.convolve_fft(spike_field, kernel, boundary="wrap")

plot_3d_slices(spike_field, N//2, N//2, N//2,
               "../output/plots/fft_test/{}_3d_slices"
               ".jpg".format(field_name))

plot_3d_slices(conv_field, N//2, N//2, N//2,
               "../output/plots/fft_test/{}_conv_3d_slices"
               ".jpg".format(field_name))
plot_3d_slices(conv_fft_field, N//2, N//2, N//2,
               "../output/plots/fft_test/{}_conv_fft_3d_slices"
               ".jpg".format(field_name))
"""

"""
# NaN Test


N = 3
R = 1
kernel = tophat_3d(R, norm=False)
# print(kernel.array)

spike_field = np.empty((N, N, N))
for i in range(N):
    for j in range(N):
        for k in range(N):
            spike_field[i, j, k] = i * N**2. + j * N + k
spike_field[0, 1, 2] = np.nan
spike_field[2, 0, 1] = np.nan
spike_field[1, 2, 0] = np.nan
spike_field[0, 2, 1] = np.nan
spike_field[1, 0, 2] = np.nan
spike_field[2, 1, 0] = np.nan

conv_field = conv.convolve(spike_field, kernel, normalize_kernel=True,
                           boundary="wrap")
conv_fft_field = conv.convolve_fft(spike_field, kernel, normalize_kernel=True,
                                   boundary="wrap")

print("Spike Field")
print(spike_field)

print("Conv Field")
print(conv_field)

print("Conv FFT Field")
print(conv_fft_field)
"""

"""
# Gaussian test 1D
std = 4

my_kernel = guassian_1d(std)
my_arr = my_kernel.array
print(my_arr)

ap_kernel = conv.Gaussian1DKernel(std)
ap_arr = ap_kernel.array
print(ap_arr)

diff = (my_arr - ap_arr) / ap_arr
print(diff)

plt.plot(diff)
plt.show()
"""

"""
# Gaussian test 2D
std = 4

my_kernel = guassian_2d(std)
my_arr = my_kernel.array
print(my_arr)

ap_kernel = conv.Gaussian2DKernel(std)
ap_arr = ap_kernel.array
print(ap_arr)

diff = (my_arr - ap_arr) / ap_arr
print(diff)

# fig = plt.figure(figsize=(6., 6.))
# axes = fig.add_axes
ais = plt.imshow(np.log10(diff))
axes = plt.gca()
divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="5%", pad=0.1)
plt.colorbar(ais, cax=cax)

plt.show()
"""

"""
# Gaussian smoothing test 2D
std = 4
N = 128
spike_field = (np.random.rand(N, N) - 0.5) * 2.

my_kernel = guassian_2d(std)
ap_kernel = conv.Gaussian2DKernel(std)

my_conv_field = conv.convolve(spike_field, my_kernel, normalize_kernel=True,
                              boundary="wrap")
ap_conv_field = conv.convolve(spike_field, ap_kernel, normalize_kernel=True,
                              boundary="wrap")

diff = (my_conv_field - ap_conv_field) / ap_conv_field
print(diff)

# fig = plt.figure(figsize=(6., 6.))
# axes = fig.add_axes
ais = plt.imshow(np.log10(diff))
axes = plt.gca()
divider = make_axes_locatable(axes)
cax = divider.append_axes("right", size="5%", pad=0.1)
plt.colorbar(ais, cax=cax)

plt.show()
"""

"""
# Gaussian smoothing test 3D
std = 2
N = 128
spike_field = (np.random.rand(N, N, N) - 0.5) * 2.
field_name = "spike_field_gaus"

kernel = guassian_3d(std)
conv_field = conv.convolve(spike_field, kernel, normalize_kernel=True,
                           boundary="wrap")

plot_3d_slices(spike_field, N//2, N//2, N//2,
               "../../output/plots/fft_test/{}_3d_slices"
               ".jpg".format(field_name))
plot_3d_slices(conv_field, N//2, N//2, N//2,
               "../../output/plots/fft_test/{}_conv_3d_slices"
               ".jpg".format(field_name))
"""


# Continuous smoothing scale test
R = 5.4
boxsize = 1100.
N = 1100

discrete_th = tophat_2d(R)
print("Discrete TH R={}".format(R))
print(discrete_th.array)

cont_th = tophat_cont(R, boxsize / N)
print("Continuous TH R={}".format(R))
print(cont_th.array)
