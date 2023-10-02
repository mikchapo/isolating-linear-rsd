import numpy as np

from vp_plot_funcs import pairwise_vel_plot


zs = [0., 0.7, 1.5, 5., 20.]
lr_names = ["z0", "z0.7", "z1.5", "z5", "z20"]
cosmo_params = np.loadtxt("/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/info/cosmo_params.dat")

for i, z in enumerate(zs):
    vp_path1 = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_no_growth_corr_vp_v7.dat".format(z)
    vp_path2 = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_no_growth_corr_{}_vp_v7.dat".format(z, lr_names[i])
    output_path = "../output/plots/abacus_part_z{}_vp_v7.jpg".format(z)
    if z < 1.:
        pairwise_vel_plot([vp_path1, vp_path2], z, cosmo_params, output_path)
    else:
        pairwise_vel_plot([vp_path1, vp_path2], z, cosmo_params, output_path, static=False)

z = 49.
vp_path1 = "/home/mj3chapm/scratch/abacus/AbacusCosmos_1100box_products/AbacusCosmos_1100box_planck_products/AbacusCosmos_1100box_planck_FoF_halos/z0.700/particle_z{}_100000_ic_no_growth_corr_vp_v7.dat".format(z)
output_path = "../output/plots/abacus_part_z{}_vp_v7.jpg".format(z)
pairwise_vel_plot([vp_path1], z, cosmo_params, output_path, static=False)