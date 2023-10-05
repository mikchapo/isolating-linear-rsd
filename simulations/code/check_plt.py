# Compares the PLT corrections to determine which are pure velocity scalings
# v0.1.1, 2021-10-25 - Fixed bug with tuple by using np.array

# Imports
import numpy as np
from struct import unpack_from


with open("ic_default/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)
    default_v = particle_data[-3:]

with open("ic_no_growth_corr/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)
    ngc_v = particle_data[-3:]

with open("ic_ngc_nplt/ic_0", "rb") as file:
    bdata = file.read()
    particle_data = unpack_from("3h6f", bdata, offset=0)
    ngc_nplt_v = particle_data[-3:]

print("Default to ngc ratio:")
print(np.array(default_v) / np.array(ngc_v))

print("Ngc to nplt, ngc ratio:")
print(np.array(ngc_v) / np.array(ngc_nplt_v))

print("Default to nplt, ngc ratio:")
print(np.array(default_v) / np.array(ngc_nplt_v))


# Change Log
# v0.1.0, 2021-10-25 - Code started with snippet from create_part_cat.py
