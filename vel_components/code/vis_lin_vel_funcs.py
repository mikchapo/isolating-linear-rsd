# Functions to visualize different velocity components
# v0.1.2, 2023-09-04 - Added defence option

# Imports
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import os
from struct import unpack_from

from AbacusCosmos import Halos
from lin_vel_funcs import calc_vel_scaling


def load_uniform_subsample(cat_dir, field_len=100., field_height=10.,
                           x0=0., y0=0., z0=0., load_pids=False):
    # Load the halo catalogue
    cat = Halos.make_catalog_from_dir(dirname=cat_dir,
                                      load_uniform_subsample=True,
                                      load_halos=False, load_pids=load_pids)
    subsample = cat.uniform_subsample

    # These filters could be applied in a single step by multiplying them, but
    # since each removes about half of the remaining objects it will require
    # many fewer calculations to apply them one after another
    subsample = subsample[subsample['pos'][:, 0] >= x0]
    subsample = subsample[subsample['pos'][:, 0] < (x0 + field_len)]
    subsample = subsample[subsample['pos'][:, 1] >= y0]
    subsample = subsample[subsample['pos'][:, 1] < (y0 + field_len)]
    subsample = subsample[subsample['pos'][:, 2] >= (z0 - field_height / 2)]
    subsample = subsample[subsample['pos'][:, 2] < (z0 + field_height / 2)]

    return subsample


def calc_density_field(field_len=100., field_height=10., N_cells=100,
                       x0=0., y0=0., z0=0.,
                       sim_name="AbacusCosmos_1100box_planck", boxid="00-0",
                       redshift=0.7):
    cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
               "{}_{}_products/{}_{}_FoF_halos/"
               "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
                                redshift))
    load_kwargs = {"field_len": field_len, "field_height": field_height,
                   "x0": x0, "y0": y0, "z0": z0}
    subsample = load_uniform_subsample(cat_dir, **load_kwargs)

    cell_length = field_len / N_cells
    density_field = np.zeros((N_cells, N_cells))
    for i in range(subsample['pos'].shape[0]):
        density_field[int((subsample['pos'][i, 0] - x0) / cell_length),
                      int((subsample['pos'][i, 1] - y0) / cell_length)] += 1

    cell_vol = field_height * cell_length**2.
    density_field = density_field * 4.e10 / cell_vol

    output = ("{}/vel_field/density_L-{}_H-{}_N-{}_o-{}-{}-{}"
              ".dat".format(cat_dir, field_len, field_height, N_cells, x0, y0,
                            z0))
    np.savetxt(output, density_field)


def calc_part_vel_field(field_len=100., field_height=10., N_cells=20,
                        x0=0., y0=0., z0=0.,
                        sim_name="AbacusCosmos_1100box_planck", boxid="00-0",
                        redshift=0.7):
    cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
               "{}_{}_products/{}_{}_FoF_halos/"
               "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
                                redshift))
    load_kwargs = {"field_len": field_len, "field_height": field_height,
                   "x0": x0, "y0": y0, "z0": z0}
    subsample = load_uniform_subsample(cat_dir, **load_kwargs)

    cell_length = field_len / N_cells
    vel_field = np.zeros((N_cells, N_cells, 3))
    for i in range(subsample['pos'].shape[0]):
        vel_field[int((subsample['pos'][i, 0] - x0) / cell_length),
                  int((subsample['pos'][i, 1] - y0) / cell_length),
                  0] += subsample['vel'][i, 0]
        vel_field[int((subsample['pos'][i, 0] - x0) / cell_length),
                  int((subsample['pos'][i, 1] - y0) / cell_length),
                  1] += subsample['vel'][i, 1]
        vel_field[int((subsample['pos'][i, 0] - x0) / cell_length),
                  int((subsample['pos'][i, 1] - y0) / cell_length),
                  2] += 1

    vel_field[:, :, 0] = vel_field[:, :, 0] / vel_field[:, :, 2]
    vel_field[:, :, 1] = vel_field[:, :, 1] / vel_field[:, :, 2]

    np.savetxt(("{}/vel_field/vel-x_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
                              y0, z0)), vel_field[:, :, 0])
    np.savetxt(("{}/vel_field/vel-y_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
                              y0, z0)), vel_field[:, :, 1])


def calc_part_lin_vel_field(field_len=100., field_height=10., N_cells=20,
                            x0=0., y0=0., z0=0.,
                            sim_name="AbacusCosmos_1100box_planck",
                            boxid="00-0", redshift=0.7):
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                redshift))
    load_kwargs = {"field_len": field_len, "field_height": field_height,
                   "x0": x0, "y0": y0, "z0": z0}
    subsample = load_uniform_subsample(cat_dir, **load_kwargs, load_pids=True)

    pids = subsample["pid"]
    N_sample = pids.size
    subsample_array = np.empty((N_sample, 7))
    subsample_array[:, 0] = pids
    subsample_array[:, 1:4] = subsample["pos"]

    # Sort the particle IDs so only one loop of the particles is needed
    sort_map = np.argsort(pids)
    sorted_pids = pids[sort_map]
    last_sort_pid_index = 0
    tot_particles = 0
    # Boolean to make sure the final value isn't overwritten
    first_final = True

    ic_dir = ("{}/ic_ngc_nplt_z49.0".format(path_root))

    N_ic_files = len(os.listdir(ic_dir))
    for i in range(N_ic_files):
        print("Starting IC file {}".format(i))
        with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
            bdata = file.read()

            particle_data_first = unpack_from("3h6f", bdata, offset=0)
            particle_data_last = unpack_from("3h6f", bdata, offset=-32)
            N_particles = ((particle_data_last[0] -
                            particle_data_first[0]) * 1440**2 +
                           (particle_data_last[1] -
                            particle_data_first[1]) * 1440 +
                           (particle_data_last[2] -
                            particle_data_first[2]) + 1)

            for j in range(last_sort_pid_index, N_sample):
                if int(sorted_pids[j] - tot_particles) < N_particles:
                    if j == (N_sample-1):
                        print("Final particle, j={}".format(j))
                        last_sort_pid_index = j

                        if first_final:
                            particle_data = unpack_from("3h6f", bdata,
                                                        offset=int((sorted_pids[j] -
                                                                    tot_particles)*32))

                            subsample_array[sort_map[j], 4:7] = particle_data[-3:]
                            first_final = False

                    else:
                        particle_data = unpack_from("3h6f", bdata,
                                                    offset=int((sorted_pids[j] -
                                                                tot_particles)*32))

                        subsample_array[sort_map[j], 4:7] = particle_data[-3:]

                else:
                    print("Switched IC at {}".format(j))
                    last_sort_pid_index = j
                    break

        tot_particles += N_particles

    cell_length = field_len / N_cells
    vel_field = np.zeros((N_cells, N_cells, 3))
    for i in range(subsample_array.shape[0]):
        vel_field[int((subsample_array[i, 1] - x0) / cell_length),
                  int((subsample_array[i, 2] - y0) / cell_length),
                  0] += subsample_array[i, 4]
        vel_field[int((subsample_array[i, 1] - x0) / cell_length),
                  int((subsample_array[i, 2] - y0) / cell_length),
                  1] += subsample_array[i, 5]
        vel_field[int((subsample_array[i, 1] - x0) / cell_length),
                  int((subsample_array[i, 2] - y0) / cell_length),
                  2] += 1

    vel_field[:, :, 0] = vel_field[:, :, 0] / vel_field[:, :, 2]
    vel_field[:, :, 1] = vel_field[:, :, 1] / vel_field[:, :, 2]

    np.savetxt(("{}/vel_field/lin-vel-x_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
                              y0, z0)), vel_field[:, :, 0])
    np.savetxt(("{}/vel_field/lin-vel-y_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
                              y0, z0)), vel_field[:, :, 1])


def calc_part_smoothed_vel_field(field_len=100., field_height=10., N_cells=20,
                                 x0=0., y0=0., z0=0.,
                                 sim_name="AbacusCosmos_1100box_planck",
                                 boxid="00-0", redshift=0.7, N_grid=1100,
                                 smooth_type="tophat", R_smooth=5.):
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                redshift))

    vel_grid = np.load("{}/vel_grid_{}_{}_{}_smoothed_v2"
                       ".npy".format(path_root, N_grid, smooth_type, R_smooth))

    # Note: No position shift is necessary when using the v2 smoothed velocity
    # grid because the shift was excluded when generating the grid

    pshift = 0.

    cell_length = field_len / N_cells
    grid_length = 1100. / N_grid
    N_grid_field = int(field_len / grid_length)
    N_grid_height = int(field_height / grid_length)
    grid_per_cell = int(N_grid_field / N_cells)
    vel_field = np.zeros((N_cells, N_cells, 3))
    for i in range(N_grid_field):
        for j in range(N_grid_field):
            for k in range(-int(N_grid_height//2), int(N_grid_height//2)):
                if not np.isnan(vel_grid[int(pshift+x0+i), int(pshift+y0+j),
                                         int(pshift+z0+k), 0]):
                    vel_field[i//grid_per_cell, j//grid_per_cell,
                              0] += vel_grid[int(pshift+x0+i),
                                             int(pshift+y0+j),
                                             int(pshift+z0+k), 0]
                    vel_field[i//grid_per_cell, j//grid_per_cell,
                              1] += vel_grid[int(pshift+x0+i),
                                             int(pshift+y0+j),
                                             int(pshift+z0+k), 1]
                    vel_field[i//grid_per_cell, j//grid_per_cell,
                              2] += 1

    for i in range(N_cells):
        for j in range(N_cells):
            if vel_field[i, j, 2] != 0:
                vel_field[i, j, 0] = vel_field[i, j, 0] / vel_field[i, j, 2]
                vel_field[i, j, 1] = vel_field[i, j, 1] / vel_field[i, j, 2]

    np.savetxt(("{}/vel_field/smooth-vel-x_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
                              y0, z0)), vel_field[:, :, 0])
    np.savetxt(("{}/vel_field/smooth-vel-y_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, N_cells, x0,
                              y0, z0)), vel_field[:, :, 1])


def plot_density_field(input_path, output_path):
    density_field = np.loadtxt(input_path)

    plt.figure(figsize=(8., 8.))
    plt.imshow(density_field, origin="lower")
    plt.savefig(output_path)


def plot_density_vel_fields(field_len=100., field_height=10., dens_cells=100,
                            vel_cells=20, x0=0., y0=0., z0=0.,
                            sim_name="AbacusCosmos_1100box_planck",
                            boxid="00-0", redshift=0.7):
    cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
               "{}_{}_products/{}_{}_FoF_halos/"
               "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
                                redshift))
    dens_input = ("{}/vel_field/density_L-{}_H-{}_N-{}_o-{}-{}-{}"
                  ".dat".format(cat_dir, field_len, field_height, dens_cells,
                                x0, y0, z0))
    density_field = np.loadtxt(dens_input)
    vx_input = ("{}/vel_field/vel-x_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, vel_cells,
                              x0, y0, z0))
    vx_field = np.loadtxt(vx_input)
    vy_input = ("{}/vel_field/vel-y_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, vel_cells,
                              x0, y0, z0))
    vy_field = np.loadtxt(vy_input)

    vel_cell_length = dens_cells / vel_cells
    q_coords = np.linspace(vel_cell_length / 2.,
                           dens_cells - vel_cell_length / 2., vel_cells)

    plt.figure(figsize=(8., 8.))
    plt.imshow(density_field, origin="lower")
    plt.quiver(q_coords, q_coords, vx_field, vy_field)

    output_path = ("{}/vel_field/dens_vel_L-{}_H-{}_N-{}_o-{}-{}-{}"
                   ".jpg".format(cat_dir, field_len, field_height, vel_cells,
                                 x0, y0, z0))

    plt.savefig(output_path)


def plot_density_lin_vel_fields(field_len=100., field_height=10.,
                                dens_cells=100, vel_cells=20, x0=0., y0=0.,
                                z0=0., sim_name="AbacusCosmos_1100box_planck",
                                boxid="00-0", redshift=0.7):
    cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
               "{}_{}_products/{}_{}_FoF_halos/"
               "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
                                redshift))
    dens_input = ("{}/vel_field/density_L-{}_H-{}_N-{}_o-{}-{}-{}"
                  ".dat".format(cat_dir, field_len, field_height, dens_cells,
                                x0, y0, z0))
    density_field = np.loadtxt(dens_input)
    vx_input = ("{}/vel_field/lin-vel-x_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, vel_cells,
                              x0, y0, z0))
    vx_field = np.loadtxt(vx_input)
    vy_input = ("{}/vel_field/lin-vel-y_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, vel_cells,
                              x0, y0, z0))
    vy_field = np.loadtxt(vy_input)

    vel_cell_length = dens_cells / vel_cells
    q_coords = np.linspace(vel_cell_length / 2.,
                           dens_cells - vel_cell_length / 2., vel_cells)

    plt.figure(figsize=(8., 8.))
    plt.imshow(density_field, origin="lower")
    plt.quiver(q_coords, q_coords, vx_field, vy_field)

    output_path = ("{}/vel_field/dens_lin-vel_L-{}_H-{}_N-{}_o-{}-{}-{}"
                   ".jpg".format(cat_dir, field_len, field_height, vel_cells,
                                 x0, y0, z0))

    plt.savefig(output_path)


def plot_density_smoothed_vel_fields(field_len=100., field_height=10.,
                                     dens_cells=100, vel_cells=20, x0=0.,
                                     y0=0., z0=0.,
                                     sim_name="AbacusCosmos_1100box_planck",
                                     boxid="00-0", redshift=0.7):
    cat_dir = ("/home/mj3chapm/scratch/abacus/{}_products/"
               "{}_{}_products/{}_{}_FoF_halos/"
               "z{:.3f}".format(sim_name, sim_name, boxid, sim_name, boxid,
                                redshift))
    dens_input = ("{}/vel_field/density_L-{}_H-{}_N-{}_o-{}-{}-{}"
                  ".dat".format(cat_dir, field_len, field_height, dens_cells,
                                x0, y0, z0))
    density_field = np.loadtxt(dens_input)
    vx_input = ("{}/vel_field/smooth-vel-x_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, vel_cells,
                              x0, y0, z0))
    vx_field = np.loadtxt(vx_input)
    vy_input = ("{}/vel_field/smooth-vel-y_L-{}_H-{}_N-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, vel_cells,
                              x0, y0, z0))
    vy_field = np.loadtxt(vy_input)

    vel_cell_length = dens_cells / vel_cells
    q_coords = np.linspace(vel_cell_length / 2.,
                           dens_cells - vel_cell_length / 2., vel_cells)

    plt.figure(figsize=(8., 8.))
    plt.imshow(density_field, origin="lower")
    plt.quiver(q_coords, q_coords, vx_field, vy_field)

    output_path = ("{}/vel_field/dens_smooth-vel_L-{}_H-{}_N-{}_o-{}-{}-{}"
                   ".jpg".format(cat_dir, field_len, field_height, vel_cells,
                                 x0, y0, z0))

    plt.savefig(output_path)


def calc_all_part_vel(field_len=100., field_height=5., x0=50., y0=50., z0=0.,
                      sim_name="AbacusCosmos_1100box_planck", boxid="00-0",
                      z_cat=0.7, N_grid=1100, smooth_type="tophat",
                      R_smooth=5.):
    '''
    This function loads the uniform particle subsample of a given simulation
    slice, takes a spatial slice of given dimensions at a given location, and
    finds the velocity, linear velocity, and smoothed linear velocity of each
    particle in that slice. It then saves those values in a catalogue.

    Parameters:
    field_len - Length in the x and y directions of the spatial slice
    field_height - Total height of the field
    x0 - x coordinate of the origin of the field. Field includes [x, x+L_f]
    y0 - y coordinate of the origin of the field. Field includes [y, y+L_f]
    z0 - z coordinate of the origin of the field. Field includes [z, z + H_f]
    sim_name - The Abacus simulation being used
    boxid - The ID of the box being used
    z_cat - The redshift of the simulation slice being used
    N_grid - The number of grid cells in 1D for the velocity smoothing
    smooth_type - The kernel used for velocity smoothing
    R_smooth - The radius used for velocity smoothing

    Output:
    Saves the new catalogue at

    '''
    # Load the uniform subsample
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                z_cat))
    load_kwargs = {"field_len": field_len, "field_height": field_height,
                   "x0": x0, "y0": y0, "z0": z0}
    subsample = load_uniform_subsample(cat_dir, **load_kwargs, load_pids=True)

    # Convert the subsample object to a numpy array
    boxsize = 1100.
    pids = subsample["pid"]
    N_sample = pids.size
    subsample_array = np.empty((N_sample, 13))
    subsample_array[:, 0] = pids
    # Note: The AbacusCosmos particle loading function arranges the particles
    # between [-boxsize, boxsize] by default. Here I've hard coded a shift to
    # change the range to [0, boxsize]
    # subsample_array[:, 1:4] = subsample["pos"] + boxsize / 2.
    subsample_array[:, 1:4] = subsample["pos"]
    subsample_array[:, 4:7] = subsample["vel"]
    del subsample

    # Initialize some variables for loading the initial conditions
    # Sort the particle IDs so only one loop of the particles is needed
    sort_map = np.argsort(pids)
    sorted_pids = pids[sort_map]
    last_sort_pid_index = 0
    tot_particles = 0
    # Boolean to make sure the final value isn't overwritten
    first_final = True
    ic_dir = ("{}/ic_ngc_nplt_z49.0".format(path_root))
    N_ic_files = len(os.listdir(ic_dir))

    # Loop through the initial conditions files and find the linear velocity
    # for each particle
    for i in range(N_ic_files):
        print("Starting IC file {}".format(i))
        with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
            bdata = file.read()

            particle_data_first = unpack_from("3h6f", bdata, offset=0)
            particle_data_last = unpack_from("3h6f", bdata, offset=-32)
            N_particles = ((particle_data_last[0] -
                            particle_data_first[0]) * 1440**2 +
                           (particle_data_last[1] -
                            particle_data_first[1]) * 1440 +
                           (particle_data_last[2] -
                            particle_data_first[2]) + 1)

            for j in range(last_sort_pid_index, N_sample):
                if int(sorted_pids[j] - tot_particles) < N_particles:
                    if j == (N_sample-1):
                        print("Final particle, j={}".format(j))
                        last_sort_pid_index = j

                        if first_final:
                            particle_data = unpack_from("3h6f", bdata,
                                                        offset=int((sorted_pids[j] -
                                                                    tot_particles)*32))

                            subsample_array[sort_map[j], 7:10] = particle_data[-3:]
                            first_final = False

                    else:
                        particle_data = unpack_from("3h6f", bdata,
                                                    offset=int((sorted_pids[j] -
                                                                tot_particles)*32))

                        subsample_array[sort_map[j], 7:10] = particle_data[-3:]

                else:
                    print("Switched IC at {}".format(j))
                    last_sort_pid_index = j
                    break

        tot_particles += N_particles

    # Scale the initial condition linear velocities to the catalogue redshift
    # Load the cosmological parameters of the simulation box
    info_path = ("{}/{}_{}_rockstar_halos/"
                 "info".format(path_root, sim_name, boxid))
    # Calculate the scaling that needs to be applied to the initial condition
    # particle velocities
    z_ic = 49.0
    vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)
    subsample_array[:, 7:10] *= vel_scaling

    # Load the smoothed velocity grid
    vel_grid = np.load("{}/vel_grid_{}_{}_{}_smoothed_v2"
                       ".npy".format(path_root, N_grid, smooth_type, R_smooth))
    l_grid = boxsize / N_grid

    for i in range(subsample_array.shape[0]):
        try:
            subsample_array[i, 10] = vel_grid[int(subsample_array[i, 1] //
                                                  l_grid),
                                              int(subsample_array[i, 2] //
                                                  l_grid),
                                              int(subsample_array[i, 3] //
                                                  l_grid), 0]
            subsample_array[i, 11] = vel_grid[int(subsample_array[i, 1] //
                                                  l_grid),
                                              int(subsample_array[i, 2] //
                                                  l_grid),
                                              int(subsample_array[i, 3] //
                                                  l_grid), 1]
            subsample_array[i, 12] = vel_grid[int(subsample_array[i, 1] //
                                                  l_grid),
                                              int(subsample_array[i, 2] //
                                                  l_grid),
                                              int(subsample_array[i, 3] //
                                                  l_grid), 2]
        except IndexError:
            print("IndexError for subsample {} with position ({}, {}, {}"
                  ")".format(i, subsample_array[i, 1], subsample_array[i, 2],
                             subsample_array[i, 3]))
            x_index = int(subsample_array[i, 1] // l_grid)
            y_index = int(subsample_array[i, 2] // l_grid)
            z_index = int(subsample_array[i, 3] // l_grid)
            if x_index >= N_grid:
                x_index = -1
            if y_index >= N_grid:
                y_index = -1
            if z_index >= N_grid:
                z_index = -1
            subsample_array[i, 10] = vel_grid[x_index, y_index, z_index, 0]
            subsample_array[i, 11] = vel_grid[x_index, y_index, z_index, 1]
            subsample_array[i, 12] = vel_grid[x_index, y_index, z_index, 2]

    np.savetxt(("{}/vel_field/all_part_vel_L-{}_H-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, x0, y0, z0)),
               subsample_array)


'''
def plot_all_part_vel(field_len=100., field_height=5., x0=50., y0=50., z0=0.,
                      sim_name="AbacusCosmos_1100box_planck", boxid="00-0",
                      z_cat=0.7, sample_name=None):
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                z_cat))
    if sample_name is None:
        sample_input = ("{}/vel_field/all_part_vel_L-{}_H-{}_o-{}-{}-{}"
                        ".dat".format(cat_dir, field_len, field_height,
                                      x0, y0, z0))
    else:
        sample_input = ("{}/vel_field/all_part_vel_{}"
                        ".dat".format(cat_dir, sample_name))
    sample = np.loadtxt(sample_input)

    info_path = ("{}/{}_{}_rockstar_halos/"
                 "info".format(path_root, sim_name, boxid))
    cosmo_params = np.loadtxt("{}/cosmo_params.dat".format(info_path))
    Omega_M = (cosmo_params[1] + cosmo_params[2]) / (cosmo_params[0] / 100.)**2
    Omega_L = 1. - Omega_M
    # This is the inversion of the factor to convert velocities to redshift
    # space displacements, since the scale parameter is used to divide the
    # quiver length
    disp_conv = ((Omega_M * (1. + z_cat)**3. + Omega_L)**0.5 * 100. /
                 (1. + z_cat))
    print("Displacement conversion:", disp_conv)
    print("Mean absolute total x velocity:", np.mean(np.abs(sample[:, 4])))
    print("Mean absolute total x displacement:",
          np.mean(np.abs(sample[:, 4]))/disp_conv)
    print("Mean absolute linear x velocity:", np.mean(np.abs(sample[:, 7])))
    print("Mean absolute linear x displacement:",
          np.mean(np.abs(sample[:, 7]))/disp_conv)
    print("Mean absolute smoothed x velocity:", np.mean(np.abs(sample[:, 10])))
    print("Mean absolute smoothed x displacement:",
          np.mean(np.abs(sample[:, 10]))/disp_conv)

    fig, axes = plt.subplots(1, 3, figsize=(9, 3.4), dpi=300, sharey=True)

    if field_len == 25.:
        markersize = 4.
        head_scaling = 2.
    elif field_len == 50.:
        markersize = 2.
        head_scaling = 1.5
    else:
        markersize = 1.
        head_scaling = 1.

    axes[0].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0].quiver(sample[:, 1], sample[:, 2], sample[:, 4],
                   sample[:, 5], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[1].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1].quiver(sample[:, 1], sample[:, 2], sample[:, 7],
                   sample[:, 8], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[2].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[2].quiver(sample[:, 1], sample[:, 2], sample[:, 10],
                   sample[:, 11], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[0].set_xlim(x0, x0+field_len)
    axes[0].set_ylim(y0, y0+field_len)
    axes[0].set_title("Total Velocities")
    axes[1].set_xlim(x0, x0+field_len)
    axes[1].set_ylim(y0, y0+field_len)
    axes[1].set_title("Linear Velocities")
    axes[2].set_xlim(x0, x0+field_len)
    axes[2].set_ylim(y0, y0+field_len)
    axes[2].set_title("Smoothed Linear Velocities")

    axes[0].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[1].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[2].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")

    axes[0].set_ylabel(r"$y\,[h^{-1} {\rm Mpc}]$")

    plt.tight_layout()

    output_path = ("{}/vel_field/all_part-vel_L-{}_H-{}_o-{}-{}-{}"
                   ".jpg".format(cat_dir, field_len, field_height,
                                 x0, y0, z0))

    plt.savefig(output_path)
'''


def plot_all_part_vel_defence(field_len=100., field_height=5., x0=50., y0=50.,
                              z0=0., sim_name="AbacusCosmos_1100box_planck",
                              boxid="00-0", z_cat=0.7, sample_name=None):
    """Visualize particle velocities."""
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                z_cat))
    if sample_name is None:
        sample_input = ("{}/vel_field/all_part_vel_L-{}_H-{}_o-{}-{}-{}"
                        ".dat".format(cat_dir, field_len, field_height,
                                      x0, y0, z0))
    else:
        sample_input = ("{}/vel_field/all_part_vel_{}"
                        ".dat".format(cat_dir, sample_name))
    sample = np.loadtxt(sample_input)

    info_path = ("{}/{}_{}_rockstar_halos/"
                 "info".format(path_root, sim_name, boxid))
    cosmo_params = np.loadtxt("{}/cosmo_params.dat".format(info_path))
    Omega_M = (cosmo_params[1] + cosmo_params[2]) / (cosmo_params[0] / 100.)**2
    Omega_L = 1. - Omega_M
    # This is the inversion of the factor to convert velocities to redshift
    # space displacements, since the scale parameter is used to divide the
    # quiver length
    disp_conv = ((Omega_M * (1. + z_cat)**3. + Omega_L)**0.5 * 100. /
                 (1. + z_cat))
    print("Displacement conversion:", disp_conv)
    print("Mean absolute total x velocity:", np.mean(np.abs(sample[:, 4])))
    print("Mean absolute total x displacement:",
          np.mean(np.abs(sample[:, 4]))/disp_conv)
    print("Mean absolute linear x velocity:", np.mean(np.abs(sample[:, 7])))
    print("Mean absolute linear x displacement:",
          np.mean(np.abs(sample[:, 7]))/disp_conv)
    print("Mean absolute smoothed x velocity:", np.mean(np.abs(sample[:, 10])))
    print("Mean absolute smoothed x displacement:",
          np.mean(np.abs(sample[:, 10]))/disp_conv)

    fig, axes = plt.subplots(1, 3, figsize=(8, 3), dpi=300, sharey=True)

    if field_len <= 25.:
        markersize = 4.
        head_scaling = 2.
    elif field_len <= 50.:
        markersize = 2.
        head_scaling = 1.5
    else:
        markersize = 1.
        head_scaling = 1.

    axes[0].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0].quiver(sample[:, 1], sample[:, 2], sample[:, 4],
                   sample[:, 5], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling,
                   headaxislength=4.5*head_scaling)

    axes[1].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1].quiver(sample[:, 1], sample[:, 2], sample[:, 10],
                   sample[:, 11], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling,
                   headaxislength=4.5*head_scaling)

    axes[2].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[2].quiver(sample[:, 1], sample[:, 2], sample[:, 4] - sample[:, 10],
                   sample[:, 5] - sample[:, 11], angles="xy",
                   scale_units="xy", scale=disp_conv, alpha=0.5,
                   headwidth=3*head_scaling, headlength=5*head_scaling,
                   headaxislength=4.5*head_scaling)

    axes[0].set_xlim(x0, x0+field_len)
    axes[0].set_ylim(y0, y0+field_len)
    axes[0].set_title("Total Velocities")
    axes[1].set_xlim(x0, x0+field_len)
    axes[1].set_ylim(y0, y0+field_len)
    axes[1].set_title("Smoothed Linear Velocities")
    axes[2].set_xlim(x0, x0+field_len)
    axes[2].set_ylim(y0, y0+field_len)
    axes[2].set_title("Non-linear Velocities")

    axes[0].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[1].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[2].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")

    axes[0].set_ylabel(r"$y\,[h^{-1} {\rm Mpc}]$")

    plt.tight_layout()

    output_path = ("{}/vel_field/all_part-vel_NL_L-{}_H-{}_o-{}-{}-{}_defence"
                   ".png".format(cat_dir, field_len, field_height,
                                 x0, y0, z0))

    plt.savefig(output_path)


def plot_all_part_vel(field_len=100., field_height=5., x0=50., y0=50.,
                      z0=0., sim_name="AbacusCosmos_1100box_planck",
                      boxid="00-0", z_cat=0.7, sample_name=None,
                      thesis=False):
    """Visualize particle velocities."""
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                z_cat))
    if sample_name is None:
        sample_input = ("{}/vel_field/all_part_vel_L-{}_H-{}_o-{}-{}-{}"
                        ".dat".format(cat_dir, field_len, field_height,
                                      x0, y0, z0))
    else:
        sample_input = ("{}/vel_field/all_part_vel_{}"
                        ".dat".format(cat_dir, sample_name))
    sample = np.loadtxt(sample_input)

    info_path = ("{}/{}_{}_rockstar_halos/"
                 "info".format(path_root, sim_name, boxid))
    cosmo_params = np.loadtxt("{}/cosmo_params.dat".format(info_path))
    Omega_M = (cosmo_params[1] + cosmo_params[2]) / (cosmo_params[0] / 100.)**2
    Omega_L = 1. - Omega_M
    # This is the inversion of the factor to convert velocities to redshift
    # space displacements, since the scale parameter is used to divide the
    # quiver length
    disp_conv = ((Omega_M * (1. + z_cat)**3. + Omega_L)**0.5 * 100. /
                 (1. + z_cat))
    print("Displacement conversion:", disp_conv)
    print("Mean absolute total x velocity:", np.mean(np.abs(sample[:, 4])))
    print("Mean absolute total x displacement:",
          np.mean(np.abs(sample[:, 4]))/disp_conv)
    print("Mean absolute linear x velocity:", np.mean(np.abs(sample[:, 7])))
    print("Mean absolute linear x displacement:",
          np.mean(np.abs(sample[:, 7]))/disp_conv)
    print("Mean absolute smoothed x velocity:", np.mean(np.abs(sample[:, 10])))
    print("Mean absolute smoothed x displacement:",
          np.mean(np.abs(sample[:, 10]))/disp_conv)

    if thesis:
        fig_width = 6.375
    else:
        fig_width = 9

    fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_width), dpi=300, sharey=True,
                             sharex=True)

    if field_len <= 25.:
        markersize = 4.
        head_scaling = 2.
    elif field_len <= 50.:
        markersize = 2.
        head_scaling = 1.5
    else:
        markersize = 1.
        head_scaling = 1.

    axes[0, 0].plot(sample[:, 1], sample[:, 2], linestyle="none",
                    marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0, 0].quiver(sample[:, 1], sample[:, 2], sample[:, 4],
                      sample[:, 5], angles="xy", scale_units="xy",
                      scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                      headlength=5*head_scaling,
                      headaxislength=4.5*head_scaling)

    axes[0, 1].plot(sample[:, 1], sample[:, 2], linestyle="none",
                    marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0, 1].quiver(sample[:, 1], sample[:, 2], sample[:, 7],
                      sample[:, 8], angles="xy", scale_units="xy",
                      scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                      headlength=5*head_scaling,
                      headaxislength=4.5*head_scaling)

    axes[1, 0].plot(sample[:, 1], sample[:, 2], linestyle="none",
                    marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1, 0].quiver(sample[:, 1], sample[:, 2], sample[:, 10],
                      sample[:, 11], angles="xy", scale_units="xy",
                      scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                      headlength=5*head_scaling,
                      headaxislength=4.5*head_scaling)

    axes[1, 1].plot(sample[:, 1], sample[:, 2], linestyle="none",
                    marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1, 1].quiver(sample[:, 1], sample[:, 2], sample[:, 4] - sample[:, 10],
                      sample[:, 5] - sample[:, 11], angles="xy",
                      scale_units="xy", scale=disp_conv, alpha=0.5,
                      headwidth=3*head_scaling, headlength=5*head_scaling,
                      headaxislength=4.5*head_scaling)

    axes[0, 0].set_xlim(x0, x0+field_len)
    axes[0, 0].set_ylim(y0, y0+field_len)
    axes[0, 0].set_title("Total Velocities")
    axes[0, 1].set_xlim(x0, x0+field_len)
    axes[0, 1].set_ylim(y0, y0+field_len)
    axes[0, 1].set_title("Linear Velocities")
    axes[1, 0].set_xlim(x0, x0+field_len)
    axes[1, 0].set_ylim(y0, y0+field_len)
    axes[1, 0].set_title("Smoothed Linear Velocities")
    axes[1, 1].set_xlim(x0, x0+field_len)
    axes[1, 1].set_ylim(y0, y0+field_len)
    axes[1, 1].set_title("Non-linear Velocities")

    axes[1, 0].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[1, 1].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")

    axes[0, 0].set_ylabel(r"$y\,[h^{-1} {\rm Mpc}]$")
    axes[1, 0].set_ylabel(r"$y\,[h^{-1} {\rm Mpc}]$")

    plt.tight_layout()

    if thesis:
        output_path = ("{}/vel_field/all_part-vel_NL_L-{}_H-{}_o-{}-{}-{}"
                       "_thesis.png".format(cat_dir, field_len, field_height,
                                            x0, y0, z0))

    else:
        output_path = ("{}/vel_field/all_part-vel_NL_L-{}_H-{}_o-{}-{}-{}"
                       ".png".format(cat_dir, field_len, field_height,
                                     x0, y0, z0))

    plt.savefig(output_path)


"""
def plot_all_part_vel(field_len=100., field_height=5., x0=50., y0=50., z0=0.,
                      sim_name="AbacusCosmos_1100box_planck", boxid="00-0",
                      z_cat=0.7, sample_name=None):
    path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
                 "{}_{}_products".format(sim_name, sim_name, boxid))
    cat_dir = ("{}/{}_{}_FoF_halos/"
               "z{:.3f}".format(path_root, sim_name, boxid,
                                z_cat))
    if sample_name is None:
        sample_input = ("{}/vel_field/all_part_vel_L-{}_H-{}_o-{}-{}-{}"
                        ".dat".format(cat_dir, field_len, field_height,
                                      x0, y0, z0))
    else:
        sample_input = ("{}/vel_field/all_part_vel_{}"
                        ".dat".format(cat_dir, sample_name))
    sample = np.loadtxt(sample_input)

    info_path = ("{}/{}_{}_rockstar_halos/"
                 "info".format(path_root, sim_name, boxid))
    cosmo_params = np.loadtxt("{}/cosmo_params.dat".format(info_path))
    Omega_M = (cosmo_params[1] + cosmo_params[2]) / (cosmo_params[0] / 100.)**2
    Omega_L = 1. - Omega_M
    # This is the inversion of the factor to convert velocities to redshift
    # space displacements, since the scale parameter is used to divide the
    # quiver length
    disp_conv = ((Omega_M * (1. + z_cat)**3. + Omega_L)**0.5 * 100. /
                 (1. + z_cat))
    print("Displacement conversion:", disp_conv)
    print("Mean absolute total x velocity:", np.mean(np.abs(sample[:, 4])))
    print("Mean absolute total x displacement:",
          np.mean(np.abs(sample[:, 4]))/disp_conv)
    print("Mean absolute linear x velocity:", np.mean(np.abs(sample[:, 7])))
    print("Mean absolute linear x displacement:",
          np.mean(np.abs(sample[:, 7]))/disp_conv)
    print("Mean absolute smoothed x velocity:", np.mean(np.abs(sample[:, 10])))
    print("Mean absolute smoothed x displacement:",
          np.mean(np.abs(sample[:, 10]))/disp_conv)

    fig, axes = plt.subplots(2, 3, figsize=(9, 6), dpi=300, sharey=True,
                             sharex=True)

    if field_len == 25.:
        markersize = 4.
        head_scaling = 2.
    elif field_len == 50.:
        markersize = 2.
        head_scaling = 1.5
    else:
        markersize = 1.
        head_scaling = 1.

    axes[0, 0].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0, 0].quiver(sample[:, 1], sample[:, 2], sample[:, 4],
                   sample[:, 5], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[0, 1].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0, 1].quiver(sample[:, 1], sample[:, 2], sample[:, 7],
                   sample[:, 8], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[0, 2].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[0, 2].quiver(sample[:, 1], sample[:, 2], sample[:, 4] - sample[:, 7],
                   sample[:, 5] - sample[:, 8], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[1, 0].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1, 0].quiver(sample[:, 1], sample[:, 2], sample[:, 4],
                   sample[:, 5], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[1, 1].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1, 1].quiver(sample[:, 1], sample[:, 2], sample[:, 10],
                   sample[:, 11], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[1, 2].plot(sample[:, 1], sample[:, 2], linestyle="none",
                 marker=".", markersize=markersize, color="b", alpha=0.3)
    axes[1, 2].quiver(sample[:, 1], sample[:, 2], sample[:, 4] - sample[:, 10],
                   sample[:, 5] - sample[:, 11], angles="xy", scale_units="xy",
                   scale=disp_conv, alpha=0.5, headwidth=3*head_scaling,
                   headlength=5*head_scaling, headaxislength=4.5*head_scaling)

    axes[0, 0].set_xlim(x0, x0+field_len)
    axes[0, 0].set_ylim(y0, y0+field_len)
    axes[0, 0].set_title("Total Velocities")
    axes[0, 1].set_xlim(x0, x0+field_len)
    axes[0, 1].set_ylim(y0, y0+field_len)
    axes[0, 1].set_title("Linear Velocities")
    axes[0, 2].set_xlim(x0, x0+field_len)
    axes[0, 2].set_ylim(y0, y0+field_len)
    axes[0, 2].set_title("Non-linear Velocities")
    axes[0, 2].text(186, 142.5, "Raw Linear Velocities", size=12,
                    verticalalignment='center', rotation=270)

    axes[1, 0].set_xlim(x0, x0+field_len)
    axes[1, 0].set_ylim(y0, y0+field_len)
    axes[1, 1].set_xlim(x0, x0+field_len)
    axes[1, 1].set_ylim(y0, y0+field_len)
    axes[1, 2].set_xlim(x0, x0+field_len)
    axes[1, 2].set_ylim(y0, y0+field_len)
    axes[1, 2].text(186, 142.5, "Smoothed Linear Velocities", size=12,
                    verticalalignment='center', rotation=270)

    axes[1, 0].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[1, 1].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")
    axes[1, 2].set_xlabel(r"$x\,[h^{-1} {\rm Mpc}]$")

    axes[0, 0].set_ylabel(r"$y\,[h^{-1} {\rm Mpc}]$")
    axes[1, 0].set_ylabel(r"$y\,[h^{-1} {\rm Mpc}]$")

    plt.tight_layout()

    output_path = ("{}/vel_field/all_part-vel_NL_L-{}_H-{}_o-{}-{}-{}"
                   ".jpg".format(cat_dir, field_len, field_height,
                                 x0, y0, z0))

    plt.savefig(output_path)
"""


def test_part_calcs(field_len=100., field_height=5., x0=50., y0=50., z0=0.,
                    sim_name="AbacusCosmos_1100box_planck", boxid="00-0",
                    z_cat=0.7, N_grid=1100, smooth_type="tophat",
                    R_smooth=5.):
    # path_root = ("/home/mj3chapm/scratch/abacus/{}_products/"
    #              "{}_{}_products".format(sim_name, sim_name, boxid))
    # cat_dir = ("{}/{}_{}_FoF_halos/"
    #            "z{:.3f}".format(path_root, sim_name, boxid,
    #                             z_cat))
    # load_kwargs = {"field_len": field_len, "field_height": field_height,
    #                "x0": x0, "y0": y0, "z0": z0}
    # subsample = load_uniform_subsample(cat_dir, **load_kwargs, load_pids=True)
    subsample = np.array([([], [], 0)], dtype=[('name', 'U10'), ('age', 'i4'), ('weight', 'f4')])

    # Convert the subsample object to a numpy array
    boxsize = 1100.
    pids = subsample["pid"]
    N_sample = pids.size
    subsample_array = np.empty((N_sample, 13))
    subsample_array[:, 0] = pids
    # Note: The AbacusCosmos particle loading function arranges the particles
    # between [-boxsize, boxsize] by default. Here I've hard coded a shift to
    # change the range to [0, boxsize]
    subsample_array[:, 1:4] = subsample["pos"] + boxsize / 2.
    subsample_array[:, 4:7] = subsample["vel"]
    del subsample

    # Initialize some variables for loading the initial conditions
    # Sort the particle IDs so only one loop of the particles is needed
    sort_map = np.argsort(pids)
    sorted_pids = pids[sort_map]
    last_sort_pid_index = 0
    tot_particles = 0
    # Boolean to make sure the final value isn't overwritten
    first_final = True
    ic_dir = ("{}/ic_ngc_nplt_z49.0".format(path_root))
    N_ic_files = len(os.listdir(ic_dir))

    # Loop through the initial conditions files and find the linear velocity
    # for each particle
    for i in range(N_ic_files):
        print("Starting IC file {}".format(i))
        with open("{}/ic_{}".format(ic_dir, i), "rb") as file:
            bdata = file.read()

            particle_data_first = unpack_from("3h6f", bdata, offset=0)
            particle_data_last = unpack_from("3h6f", bdata, offset=-32)
            N_particles = ((particle_data_last[0] -
                            particle_data_first[0]) * 1440**2 +
                           (particle_data_last[1] -
                            particle_data_first[1]) * 1440 +
                           (particle_data_last[2] -
                            particle_data_first[2]) + 1)

            for j in range(last_sort_pid_index, N_sample):
                if int(sorted_pids[j] - tot_particles) < N_particles:
                    if j == (N_sample-1):
                        print("Final particle, j={}".format(j))
                        last_sort_pid_index = j

                        if first_final:
                            particle_data = unpack_from("3h6f", bdata,
                                                        offset=int((sorted_pids[j] -
                                                                    tot_particles)*32))

                            subsample_array[sort_map[j], 7:10] = particle_data[-3:]
                            first_final = False

                    else:
                        particle_data = unpack_from("3h6f", bdata,
                                                    offset=int((sorted_pids[j] -
                                                                tot_particles)*32))

                        subsample_array[sort_map[j], 7:10] = particle_data[-3:]

                else:
                    print("Switched IC at {}".format(j))
                    last_sort_pid_index = j
                    break

        tot_particles += N_particles

    # Scale the initial condition linear velocities to the catalogue redshift
    # Load the cosmological parameters of the simulation box
    info_path = ("{}/{}_{}_rockstar_halos/"
                 "info".format(path_root, sim_name, boxid))
    # Calculate the scaling that needs to be applied to the initial condition
    # particle velocities
    z_ic = 49.0
    vel_scaling = calc_vel_scaling(z_ic, z_cat, info_path)
    subsample_array[:, 7:10] *= vel_scaling

    # Load the smoothed velocity grid
    vel_grid = np.load("{}/vel_grid_{}_{}_{}_smoothed_v2"
                       ".npy".format(path_root, N_grid, smooth_type, R_smooth))
    l_grid = boxsize / N_grid

    for i in range(subsample_array.shape[0]):
        try:
            subsample_array[i, 10] = vel_grid[int(subsample_array[i, 1] //
                                                  l_grid),
                                              int(subsample_array[i, 2] //
                                                  l_grid),
                                              int(subsample_array[i, 3] //
                                                  l_grid), 0]
            subsample_array[i, 11] = vel_grid[int(subsample_array[i, 1] //
                                                  l_grid),
                                              int(subsample_array[i, 2] //
                                                  l_grid),
                                              int(subsample_array[i, 3] //
                                                  l_grid), 1]
            subsample_array[i, 12] = vel_grid[int(subsample_array[i, 1] //
                                                  l_grid),
                                              int(subsample_array[i, 2] //
                                                  l_grid),
                                              int(subsample_array[i, 3] //
                                                  l_grid), 2]
        except IndexError:
            print("IndexError for subsample {} with position ({}, {}, {}"
                  ")".format(i, subsample_array[i, 1], subsample_array[i, 2],
                             subsample_array[i, 3]))
            x_index = int(subsample_array[i, 1] // l_grid)
            y_index = int(subsample_array[i, 2] // l_grid)
            z_index = int(subsample_array[i, 3] // l_grid)
            if x_index >= N_grid:
                x_index = -1
            if y_index >= N_grid:
                y_index = -1
            if z_index >= N_grid:
                z_index = -1
            subsample_array[i, 10] = vel_grid[x_index, y_index, z_index, 0]
            subsample_array[i, 11] = vel_grid[x_index, y_index, z_index, 1]
            subsample_array[i, 12] = vel_grid[x_index, y_index, z_index, 2]

    np.savetxt(("{}/vel_field/all_part_vel_L-{}_H-{}_o-{}-{}-{}"
                ".dat".format(cat_dir, field_len, field_height, x0, y0, z0)),
               subsample_array)

# Change Log
# v0.1.1, 2023-06-07 - Added thesis formatting
# v0.1.0, 2022-08-23 - Code started
