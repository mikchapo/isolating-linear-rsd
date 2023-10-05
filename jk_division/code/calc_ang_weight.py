import numpy as np
import sys

input_root = sys.argv[1]
par_input = sys.argv[2]
fib_input = sys.argv[3]
norm_input = sys.argv[4]
output_root = sys.argv[5]
weight_output = sys.argv[6]
N_reg = int(sys.argv[7])

for i in range(N_reg+2):
    par = np.loadtxt(input_root + str(i) + par_input)
    fib = np.loadtxt(input_root + str(i) + fib_input)
    norms = np.loadtxt(input_root + str(i) + norm_input, skiprows=1)

    par[:, -1] = par[:, -1] / norms[0]
    fib[:, -1] = fib[:, -1] / norms[1]

    weight = np.copy(par)
    weight[:, -1] = weight[:, -1] / fib[:, -1]

    np.savetxt(output_root + str(i) + weight_output, weight)
