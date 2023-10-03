import os
from struct import unpack_from


for i in range(41):
    if i == 40:
        ic_dir = ("/home/mj3chapm/scratch/abacus/"
                  "AbacusCosmos_1100box_products/"
                  "AbacusCosmos_1100box_planck_products/"
                  "ic_ngc_nplt_z49.0")

    elif i < 10:
        ic_dir = ("/home/mj3chapm/scratch/abacus/"
                  "AbacusCosmos_1100box_products/"
                  "AbacusCosmos_1100box_0{}_products/"
                  "ic_ngc_nplt_z49.0".format(i))

    else:
        ic_dir = ("/home/mj3chapm/scratch/abacus/"
                  "AbacusCosmos_1100box_products/"
                  "AbacusCosmos_1100box_{}_products/"
                  "ic_ngc_nplt_z49.0".format(i))

    tot_particles = 0

    tot_ic_size = os.path.getsize("{}/ic_0".format(ic_dir))
    N_particles = int(tot_ic_size / 32)

    with open("{}/ic_0".format(ic_dir), "rb") as file:
        bdata = file.read()
        particle_data_b = unpack_from("3h6f", bdata, offset=0)
        particle_data_e = unpack_from("3h6f", bdata, offset=-32)

    N_index = ((particle_data_e[0] - particle_data_b[0]) * 1440**2 +
               (particle_data_e[1] - particle_data_b[1]) * 1440 +
               (particle_data_e[2] - particle_data_b[2]) + 1)

    print("Box {}, N_len: {},".format(i, N_particles),
          "N_ind: {}, Ratio: {}".format(N_index, N_particles / N_index))
