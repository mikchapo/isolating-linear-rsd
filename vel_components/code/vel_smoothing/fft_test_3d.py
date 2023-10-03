# Test FFt script and make sure everything is working right
# v0.1.0, 2022-04-07 - Code started

# Imports
import astropy.convolution as conv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np
from scipy.fft import fft2, ifft2


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
    N = 2 * R + 1
    win = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if np.sqrt((i - R)**2 + (j - R)**2) <= float(R):
                win[i, j] = 1.
    win /= np.pi * float(R)**2.
    return conv.Kernel(win)


def tophat_3d(R):
    N = 2 * R + 1
    win = np.zeros((N, N, N))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                if np.sqrt((i - R)**2 + (j - R)**2 + (k - R)**2) <= float(R):
                    win[i, j, k] = 1.
    win /= 4. / 3. * np.pi * float(R)**3.
    return conv.Kernel(win)


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
