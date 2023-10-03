# Script to visualize velocity field
# v0.1.0, 2022-08-23

# Imports
import matplotlib.pyplot as plt
import numpy as np
import sys

from vis_lin_vel_funcs import *


sim_name = "AbacusCosmos_1100box_planck"
boxid = "00-0"
redshift = 0.7
field_height = 10.
N_cells = 100
z0 = 0.

field_len = float(sys.argv[1])
x0 = float(sys.argv[2])
y0 = float(sys.argv[3])

# cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
#            "{}_{}_products/{}_{}_FoF_halos/"
#            "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
#                             redshift))
# input_path = ("{}/vel_field/density_L-{}_H-{}_N-{}_o-{}-{}-{}"
#               ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
#                             y0, z0))
# output_path = ("{}/vel_field/density_L-{}_H-{}_N-{}_o-{}-{}-{}"
#                ".jpg".format(cat_dir, field_len, field_height, N_cells, x0,
#                              y0, z0))

# calc_density_field()
# plot_density_field(input_path, output_path)

# calc_part_vel_field()
# plot_density_vel_fields()

# calc_part_lin_vel_field()
# plot_density_lin_vel_fields()

# calc_part_smoothed_vel_field()
# plot_density_smoothed_vel_fields()

# calc_all_part_vel(field_len=field_len, x0=x0, y0=y0)
sample_name = "L-{}_H-{}_o-{}-{}-{}".format(100.0, 5.0, 100.0, 75.0, 0.0) # Load the largest calculated field, even when only plotting a subset
plot_all_part_vel_defence(field_len=field_len, x0=x0, y0=y0, sample_name=sample_name)
# plot_all_part_vel()
