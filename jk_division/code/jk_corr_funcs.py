import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import sys
from seaborn import color_palette
from astropy.table import Table, vstack


def landy_szalay(DD, DR, RR):
    corr_func = np.copy(DD)
    corr_func[:,2] = (np.divide(DD[:,2], RR[:,2], out=np.zeros_like(DD[:,2]), where=RR[:,2]!=0) - 2. *
                      np.divide(DR[:,2], RR[:,2], out=np.zeros_like(DR[:,2]), where=RR[:,2]!=0) + 1.)
    return corr_func


def L_0(x):
    return 1.


def L_2(x):
    return 0.5 * (3. * x * x - 1.)


def L_4(x):
    return (35.* x * x * x * x - 30. * x * x + 3.) / 8.


L_polynomials = [L_0, L_2, L_4]


def corr_func_multipole_calc(l, corr_func, mu_bins=100):
    L = L_polynomials[np.floor(l / 2).astype(int)]
    
    corr_func_multipole = np.zeros((np.unique(corr_func[:,0]).size, 2))
    corr_func_multipole[:,0] = np.unique(corr_func[:,0])
    for multi in corr_func_multipole:
        for corr in corr_func:
            if corr[0] == multi[0]:
    
                multi[1] += (1./ float(mu_bins)) * corr[2] * L(corr[1] + 0.5 / float(mu_bins))
    
    corr_func_multipole[:,1] = corr_func_multipole[:,1] * (2. * l + 1.)
    
    return corr_func_multipole


def corr_func_projected_calc(corr_func, dr_pi=1., r_pi_max=80.):
    corr_func_projected = np.zeros((np.unique(corr_func[:,0]).size, 2))
    corr_func_projected[:,0] = np.unique(corr_func[:,0])
#     corr_func_cutoff = np.zeros((np.unique(corr_func[:,0]).size, 2))
#     corr_func_cutoff[:,0] = np.unique(corr_func[:,0])
    
    for proj in corr_func_projected:
        for corr in corr_func:
            if corr[0] == proj[0] and corr[1] < r_pi_max:
                proj[1] += dr_pi * corr[2]
    
    corr_func_projected[:,1] = corr_func_projected[:,1] * 2.
    
    return corr_func_projected


def load_pair_counts(filename):
    return np.loadtxt(filename, dtype=np.float, delimiter='\t')


def calc_step_size(binning):
    step_sizes = []
    for dim in binning:
        if dim[3] == 'log':
            step_size = (np.log10(dim[1]) - np.log10(dim[0])) / dim[2]
            
        else:
            step_size = (dim[1] - dim[0]) / dim[2]
        
        step_sizes.append(step_size)
        
    return step_sizes


def generate_empty_bins(binning):
    binning = np.array(binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
    column_bases = []
    for dim in binning:
        if dim['scaling'] == 'log':
            column_bases.append(np.logspace(np.log10(dim['min']), np.log10(dim['max']), dim['N_bin'], endpoint=False))
            
        elif dim['scaling'] == 'lin':
            column_bases.append(np.linspace(dim['min'], dim['max'], dim['N_bin'], endpoint=False))
                       
        else:
            print("Invalid scaling type")
            return False

    N_bins_tot = np.prod(binning['N_bin'])
    step_sizes = calc_step_size(binning)
    
    columns = []
    for i in range(binning.shape[0]):
#         print("Starting dim %i" % i)            
        column = np.tile(column_bases[i], (np.prod(binning['N_bin'][i+1:]),1))
#         print("Shape after initial tile:", column.shape)
        column = column.flatten(order='F')
#         print("Shape after flatten:", column.shape)
        column = np.tile(column, np.prod(binning['N_bin'][:i]))
#         print("Shape after final tile:", column.shape)
            
        columns.append(column)

    columns.append(np.zeros(N_bins_tot))
    
    bin_array = np.column_stack(columns)

    return bin_array

                    
def rebin_pair_counts(original, compression_factor, binning):
    binning = np.array(binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
    new_binning = np.copy(binning)
    new_binning["N_bin"][0] = binning["N_bin"][0] / compression_factor
    rebinned = generate_empty_bins(new_binning)
    for i in range(binning["N_bin"][0]):
        i_new = int(i / compression_factor)
        for j in range(binning["N_bin"][1]):
            rebinned[i_new*new_binning["N_bin"][1]+j, 2] += (
                        original[i*binning["N_bin"][1]+j, 2])
    return rebinned


def jk_corr_funcs_test(input_root, dv_input, run_name, references=[], reference_labels=[], N_reg=200, include_full=False, output_path="../output/jk_test_"):
    textsize = 18
    figsize = [12., 9.]
    linewidth = 2
    s_pow = 0.
    
    all_regions = np.loadtxt(input_root + str(N_reg) + dv_input)
    if include_full:
        full = np.loadtxt(input_root + str(N_reg+1) + dv_input)
    if references:
        for i in range(len(references)):
            references[i] = np.loadtxt(references[i])
    
    log_dsep = np.log10(all_regions[1,0]) - np.log10(all_regions[0,0])
    seps = np.power(10., (np.log10(all_regions[:9,0]) + log_dsep / 2.))

    fig = plt.figure(figsize=figsize) 
    for i in range(N_reg):
        region = np.loadtxt(input_root + str(i) + dv_input)
        if i==0:
            plt.semilogx(seps, np.power(seps, 2)*region[:9,1], label="Mono JK Realizations",
                     linestyle="-", color=color_palette('colorblind')[2], linewidth=linewidth, alpha=0.3)
            plt.semilogx(seps, np.power(seps, 2)*region[9:18,1], label="Quad JK Realizations",
                     linestyle="-", color=color_palette('colorblind')[3], linewidth=linewidth, alpha=0.3)
        else:
            plt.semilogx(seps, np.power(seps, 2)*region[:9,1],
                     linestyle="-", color=color_palette('colorblind')[2], linewidth=linewidth, alpha=0.3)
            plt.semilogx(seps, np.power(seps, 2)*region[9:18,1],
                     linestyle="-", color=color_palette('colorblind')[3], linewidth=linewidth, alpha=0.3)
    
    plt.semilogx(seps, np.power(seps, 2)*all_regions[:9,1],
                 linestyle="-", label="Mono All Regions", color=color_palette('colorblind')[0], linewidth=2)
    plt.semilogx(seps, np.power(seps, 2)*all_regions[9:18,1],
                 linestyle="-", label="Quad All Regions",color=color_palette('colorblind')[1], linewidth=2)
    
    if references:
        for i in range(len(references)):
            plt.semilogx(seps, np.power(seps, 2)*references[i][:9,1],
                         linestyle="--", label="Mono " + reference_labels[i], color=color_palette('colorblind')[6 + i*2], linewidth=2)
            plt.semilogx(seps, np.power(seps, 2)*references[i][9:18,1],
                         linestyle="--", label="Quad " + reference_labels[i], color=color_palette('colorblind')[7 + i*2], linewidth=2)

    if include_full:
        plt.semilogx(seps, np.power(seps, 2)*full[:9,1],
                     linestyle=":", label="Mono Full", color=color_palette('colorblind')[4], linewidth=2)
        plt.semilogx(seps, np.power(seps, 2)*full[9:18,1],
                     linestyle=":", label="Quad Full",color=color_palette('colorblind')[5], linewidth=2)
    
    plt.title('%s JK correlation function multipoles' % run_name, fontsize=textsize)
    plt.xlabel(r'$s[h^{-1}Mpc]$', fontsize=textsize)
    plt.ylabel(r'$s^2*\xi_l(s)$', fontsize=textsize)
    plt.xlim(0.1,60)
    plt.xscale('log')
    plt.legend(fontsize=textsize)
    plt.savefig(output_path + "mps.png")
    
    fig = plt.figure(figsize=figsize)
    for i in range(N_reg):
        region = np.loadtxt(input_root + str(i) + dv_input)
        if i==0:
            plt.semilogx(seps, region[18:,1], linestyle="-", color=color_palette('colorblind')[1], linewidth=0.5,
                        alpha=0.3, label="JK Realizations")
        else:
            plt.semilogx(seps, region[18:,1], linestyle="-", color=color_palette('colorblind')[1], linewidth=0.5,
                        alpha=0.3)
    plt.semilogx(seps, all_regions[18:,1], linestyle="--",
                    label="All Regions", color=color_palette('colorblind')[0], linewidth=2)
    
    if references:
        for i in range(len(references)):
            plt.semilogx(seps, references[i][18:,1], linestyle="--",
                            label=reference_labels[i], color=color_palette('colorblind')[3 + i*2], linewidth=2)

    if include_full:
        plt.semilogx(seps, full[18:,1], linestyle=":",
                        label="Full", color=color_palette('colorblind')[2], linewidth=2)

    plt.title('%s JK projected correlation functions' % run_name, fontsize=textsize)
    plt.xlabel(r'$r_\perp[h^{-1}Mpc]$', fontsize=textsize)
    plt.ylabel(r'$w_p(r_\perp)$', fontsize=textsize)
    plt.xlim(0.1,60)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(fontsize=textsize)
    plt.savefig(output_path + "wp.png")

    fig = plt.figure(figsize=figsize)
    mean_mono = np.zeros(9)
    for i in range(N_reg):
        region = np.loadtxt(input_root + str(i) + dv_input)
        mean_mono += region[:9,1]
        if i==0:
            plt.semilogx(seps, (region[:9,1]-all_regions[:9,1])/all_regions[:9,1], label="JK Realizations",
                     linestyle="-", color=color_palette('colorblind')[1], linewidth=0.5, alpha=0.5)
        else:
            plt.semilogx(seps, (region[:9,1]-all_regions[:9,1])/all_regions[:9,1],
                     linestyle="-", color=color_palette('colorblind')[1], linewidth=0.5, alpha=0.5)

    mean_mono = mean_mono / N_reg
    plt.semilogx(seps, (mean_mono-all_regions[:9,1])/all_regions[:9,1], label="JK Mean",
                     linestyle="-", color=color_palette('colorblind')[0], linewidth=2)
    if references:
        for i in range(len(references)):
            plt.semilogx(seps, (references[i][:9,1]-all_regions[:9,1])/all_regions[:9,1], label=reference_labels[i],
                         linestyle="-", color=color_palette('colorblind')[3 + i], linewidth=2)

    if include_full:
        plt.semilogx(seps, (full[:9,1]-all_regions[:9,1])/all_regions[:9,1], label="Full",
                     linestyle="--", color=color_palette('colorblind')[2], linewidth=2)

    plt.semilogx(seps, [0,0,0,0,0,0,0,0,0], "k-", linewidth=2)
    plt.title('%s JK correlation function monopole relative difference from all regions' % run_name,
              fontsize=textsize)
    plt.xlabel(r'$s[h^{-1}Mpc]$', fontsize=textsize)
    plt.ylabel(r'$\frac{\xi_{0,JK}-\xi_{0,A}}{\xi_{0,A}}$', fontsize=textsize)
    plt.xlim(0.1,60)
    plt.ylim(-0.2, 0.2)
    plt.xscale('log')
    plt.legend(fontsize=textsize)
    plt.savefig(output_path + "relative_mono.png")

    fig = plt.figure(figsize=figsize)
    mean_quad = np.zeros(9)
    for i in range(N_reg):
        region = np.loadtxt(input_root + str(i) + dv_input)
        mean_quad += region[9:18,1]
        if i==0:
            plt.semilogx(seps, (region[9:18,1]-all_regions[9:18,1])/all_regions[9:18,1], label="JK Realizations",
                     linestyle="-", color=color_palette('colorblind')[1], linewidth=0.5, alpha=0.5)
        else:
            plt.semilogx(seps, (region[9:18,1]-all_regions[9:18,1])/all_regions[9:18,1],
                     linestyle="-", color=color_palette('colorblind')[1], linewidth=0.5, alpha=0.5)

    mean_quad = mean_quad / N_reg
    plt.semilogx(seps, (mean_quad-all_regions[9:18,1])/all_regions[9:18,1], label="JK Mean",
                     linestyle="-", color=color_palette('colorblind')[0], linewidth=2)
    if references:
        for i in range(len(references)):
            plt.semilogx(seps, (references[i][9:18,1]-all_regions[9:18,1])/all_regions[9:18,1], label=reference_labels[i],
                         linestyle="-", color=color_palette('colorblind')[3 + i], linewidth=2)

    if include_full:
        plt.semilogx(seps, (full[9:18,1]-all_regions[9:18,1])/all_regions[9:18,1], label="Full",
                     linestyle="--", color=color_palette('colorblind')[2], linewidth=2)


    plt.semilogx(seps, [0,0,0,0,0,0,0,0,0], "k-", linewidth=2)
    plt.title('%s JK correlation function quadrupole relative difference from all regions' % run_name,
              fontsize=textsize)
    plt.xlabel(r'$s[h^{-1}Mpc]$', fontsize=textsize)
    plt.ylabel(r'$\frac{\xi_{2,JK}-\xi_{2,A}}{\xi_{2,A}}$', fontsize=textsize)
    plt.xlim(0.1,60)
    plt.ylim(-0.2, 0.2)
    plt.xscale('log')
    plt.legend(fontsize=textsize)
    plt.savefig(output_path + "relative_quad.png")
    
    fig = plt.figure(figsize=figsize)
    mean_wp = np.zeros(9)
    for i in range(N_reg):
        region = np.loadtxt(input_root + str(i) + dv_input)
        mean_wp += region[18:,1]
        if i==0:
            plt.semilogx(seps, (region[18:,1]-all_regions[18:,1])/all_regions[18:,1], linestyle="-",
                         color=color_palette('colorblind')[1], linewidth=0.5, alpha=0.5, label="JK Realizations")
        else:
            plt.semilogx(seps, (region[18:,1]-all_regions[18:,1])/all_regions[18:,1], linestyle="-",
                         color=color_palette('colorblind')[1], linewidth=0.5, alpha=0.5)

    mean_wp = mean_wp / N_reg
    plt.semilogx(seps, (mean_wp-all_regions[18:,1])/all_regions[18:,1], linestyle="-",
                 color=color_palette('colorblind')[0], linewidth=2, label="JK Mean")
    if references:
        for i in range(len(references)):
            plt.semilogx(seps, (references[i][18:,1]-all_regions[18:,1])/all_regions[18:,1], label=reference_labels[i],
                         linestyle="-", color=color_palette('colorblind')[3 + i], linewidth=2)

    if include_full:
        plt.semilogx(seps, (full[18:,1]-all_regions[18:,1])/all_regions[18:,1], label="Full",
                     linestyle="--", color=color_palette('colorblind')[2], linewidth=2)

    plt.semilogx(seps, [0,0,0,0,0,0,0,0,0],"k-")
    plt.title('%s JK projected correlation functions relative difference from all regions' % run_name,
              fontsize=textsize)
    plt.xlabel(r'$r_\perp[h^{-1}Mpc]$', fontsize=textsize)
    plt.ylabel(r'$\frac{w_{p,JK}-w_{p,A}}{w_{p,A}}$', fontsize=textsize)
    plt.xlim(0.1,60)
    plt.ylim(-0.2, 0.2)
    plt.xscale('log')
    plt.legend(fontsize=textsize)
    plt.savefig(output_path + "relative_wp.png")

    # if include_full:
    #     full_comparison = np.empty((27,5))
    #     full_comparison[:, 0] = full[:, 0]
    #     full_comparison[:, 1] = full[:, 1]
    #     full_comparison[:, 2] = references[i][:, 1]
    #     full_comparison[:, 3] = full[:, 1] - references[i][:, 1]
    #     full_comparison[:, 4] = (full[:, 1] - references[i][:, 1]) / references[i][:, 1]

    #     np.savetxt(output_path + "full_comparison.txt", full_comparison)


mps_binning = [(0.1, 60., 180, 'log'),
          (0., 1., 100, 'lin')]

wp_binning = [(0.1, 60., 180, 'log'),
          (0., 150., 150, 'lin')]

mps_binning = np.array(mps_binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
wp_binning = np.array(wp_binning, dtype=[('min', 'f4'), ('max', 'f4'), ('N_bin', 'i4'), ('scaling', 'U5')])
wp_step_sizes = calc_step_size(wp_binning)


input_root = sys.argv[1]
output_root = sys.argv[2]
dv_output = sys.argv[3]
N_reg = int(sys.argv[4])

cfs = ["mps", "wp"]
caps = ["NGC", "SGC"]
pcs = ["DD", "DR", "RR"]

alpha = N_reg / (2. + np.sqrt(2.) * (N_reg - 1.))

D_eff_all_reg_NGC = np.loadtxt("{}{}/DR_norm_eBOSS_LRG_NGC_v7_2_pip.dat".format(input_root, N_reg), skiprows=1)[0]
D_eff_all_reg_SGC = np.loadtxt("{}{}/DR_norm_eBOSS_LRG_SGC_v7_2_pip.dat".format(input_root, N_reg), skiprows=1)[0]
D_eff_all_reg = D_eff_all_reg_NGC + D_eff_all_reg_SGC

R_eff_all_reg_NGC = np.loadtxt("{}{}/DR_norm_eBOSS_LRG_NGC_v7_2_pip.dat".format(input_root, N_reg), skiprows=1)[1]
R_eff_all_reg_SGC = np.loadtxt("{}{}/DR_norm_eBOSS_LRG_SGC_v7_2_pip.dat".format(input_root, N_reg), skiprows=1)[1]
R_eff_all_reg = R_eff_all_reg_NGC + R_eff_all_reg_SGC

R2_eff_all_reg_NGC = np.loadtxt("{}{}/RR_norm_eBOSS_LRG_NGC_v7_2_pip.dat".format(input_root, N_reg), skiprows=1)[1]
R2_eff_all_reg_SGC = np.loadtxt("{}{}/RR_norm_eBOSS_LRG_SGC_v7_2_pip.dat".format(input_root, N_reg), skiprows=1)[1]
R2_eff_all_reg = R2_eff_all_reg_NGC + R2_eff_all_reg_SGC

for k in range(N_reg+2):
    mps_pair_counts = []
    wp_pair_counts = []
    norms = []
    cap_norms = []
    for cap in caps:
        for pc in pcs:
            mps_pair_counts.append(load_pair_counts(input_root + str(k) + "/%s_eBOSS_LRG_%s_v7_2_mps_pip%s.dat" % (pc, cap, "" if pc=="RR" else "_ang")))
            wp_pair_counts.append(load_pair_counts(input_root + str(k) + "/%s_eBOSS_LRG_%s_v7_2_wp_pip%s.dat" % (pc, cap, "" if pc=="RR" else "_ang")))
        cap_norms.append(np.loadtxt(input_root + str(k) + "/DD_norm_eBOSS_LRG_%s_v7_2_pip.dat" % cap, skiprows=1))
        cap_norms.append(np.loadtxt(input_root + str(k) + "/DR_norm_eBOSS_LRG_%s_v7_2_pip.dat" % cap, skiprows=1)[0])
        cap_norms.append(np.loadtxt(input_root + str(k) + "/DR_norm_eBOSS_LRG_%s_v7_2_pip.dat" % cap, skiprows=1)[1])
        cap_norms.append(np.loadtxt(input_root + str(k) + "/RR_norm_eBOSS_LRG_%s_v7_2_pip.dat" % cap, skiprows=1)[0])
        cap_norms.append(np.loadtxt(input_root + str(k) + "/RR_norm_eBOSS_LRG_%s_v7_2_pip.dat" % cap, skiprows=1)[1])

    DD_norm_cross = np.loadtxt(input_root + str(k) + "/DD_cross_norm_eBOSS_LRG_v7_2_pip.dat", skiprows=1)
    norms.append(cap_norms[0] + cap_norms[5] + DD_norm_cross)

    if k == N_reg+1:
        DR_norm = (cap_norms[1] + cap_norms[6]) * (cap_norms[2] + cap_norms[7])
        RR_norm = ((cap_norms[3] + cap_norms[8])**2. - (cap_norms[4] + cap_norms[9])) / 2.

    else:
        D_eff_rem = D_eff_all_reg - (cap_norms[1] + cap_norms[6])
        R_eff_rem = R_eff_all_reg - (cap_norms[2] + cap_norms[7])
        DR_norm = (D_eff_all_reg * R_eff_all_reg - alpha *
                   (D_eff_rem * R_eff_all_reg + D_eff_all_reg * R_eff_rem - 2 *
                    D_eff_rem * R_eff_rem) - D_eff_rem * R_eff_rem)

        R2_eff_rem = R2_eff_all_reg - (cap_norms[4] + cap_norms[9])
        RR_norm = ((R_eff_all_reg**2. - R2_eff_all_reg) / 2. - alpha *
                   (R_eff_rem*R_eff_all_reg) - (R_eff_rem**2. - R2_eff_rem) /
                   2.)

    norms.append(DR_norm)
    norms.append(RR_norm)

    if k == 0:
        print(cap_norms)
        print(norms)

    for i in range(6):
        mps_pair_counts[i] = rebin_pair_counts(mps_pair_counts[i], 20, mps_binning)
        wp_pair_counts[i] = rebin_pair_counts(wp_pair_counts[i], 20, wp_binning)

        if k == 0:
            print("Last sep bin, 10 rp bin PC {}:".format(i))
            print(wp_pair_counts[i][-150:-140, -1])

    for i in range(3):
        mps_pair_counts[i][:,-1] += mps_pair_counts[i+3][:,-1]
        wp_pair_counts[i][:,-1] += wp_pair_counts[i+3][:,-1]

        if k == 0:
            print("Last sep bin, 10 rp bin PC {}, cap combined:".format(i))
            print(wp_pair_counts[i][-150:-140, -1])

        mps_pair_counts[i][:,-1] = mps_pair_counts[i][:,-1] / norms[i]
        wp_pair_counts[i][:,-1] = wp_pair_counts[i][:,-1] / norms[i]

        if k == 0:
            print("Last sep bin, 10 rp bin PC {}, normed:".format(i))
            print(wp_pair_counts[i][-150:-140, -1])

    # del mps_pair_counts[3:]
    # del wp_pair_counts[3:]

    mps_corr_func = landy_szalay(mps_pair_counts[0], mps_pair_counts[1], mps_pair_counts[2])
    wp_corr_func = landy_szalay(wp_pair_counts[0], wp_pair_counts[1], wp_pair_counts[2])

    if k == 0:
        print("Last sep bin, 10 rp bin LS:")
        print(wp_corr_func[-150:-140, -1])

    mono_corr_func = corr_func_multipole_calc(0, mps_corr_func, mu_bins=mps_binning['N_bin'][1])
    quad_corr_func = corr_func_multipole_calc(2, mps_corr_func, mu_bins=mps_binning['N_bin'][1])
    proj_corr_func = corr_func_projected_calc(wp_corr_func, dr_pi=wp_step_sizes[1], r_pi_max=80.)

    if k == 0:
        print("Last sep bin wp:")
        print(proj_corr_func[-1, 1])

    data_vector = np.concatenate((mono_corr_func, quad_corr_func, proj_corr_func))

    np.savetxt(output_root + str(k) + dv_output, data_vector)
    print("Finished %i" % k)

jk_corr_funcs_test(output_root, dv_output, "eBOSS LRG v7_2 PIP+ANG", references=["../data/dv_eBOSS_LRG_comb_v7_2_pip_ang_ab.dat"],
                   reference_labels=["Reference"], N_reg=N_reg, include_full=True, output_path="../output/jk_test_")
