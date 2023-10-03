# Plot the progress of chains for each parameter to determine burn-in
# v0.1.0, 2022-01-14 - Code copied from RSD/fit/code/plot_burn_in.py and updated

# Imports
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import sys

start_time = dt.datetime.now()
print("Started at", start_time)

N_chains = int(sys.argv[1])
chain_path = sys.argv[2]
output = sys.argv[3]

print("Loading column headers", dt.datetime.now() - start_time)

first_chain = open("%s.1.txt" % (chain_path))
column_headers = first_chain.readline()
first_chain.close()
column_headers = column_headers[1:].split()

print(column_headers)
print(len(column_headers))

N_params = len(column_headers) - 6
fig, axes = plt.subplots(nrows=N_params, sharex=True, figsize=(6, 2*N_params))

for i in range(N_chains):
    print("Loading chain %d" % i, dt.datetime.now() - start_time)
    chain = np.loadtxt("%s.%i.txt" % (chain_path, i+1))
    for j in range(N_params):
        axes[j].plot(chain[:, j+2], linewidth=0.5)

print("Setting labels", dt.datetime.now() - start_time)
for i in range(N_params):
    axes[i].set_ylabel(column_headers[i+2])

print("Saving Figure", dt.datetime.now() - start_time)
plt.savefig(output)
