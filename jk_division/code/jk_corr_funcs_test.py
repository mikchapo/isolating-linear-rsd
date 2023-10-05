import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
import sys
from seaborn import color_palette
from astropy.table import Table, vstack


def jk_corr_funcs_test(input_root, dv_input, run_name, references=[], reference_labels=[], N_reg=200, include_full=False, output_path="../output/jk_test_"):
    textsize = 18
    figsize = [14., 6.]
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
    axes = [0, 0, 0, 0, 0, 0]
    for i in range(3):
        axes[2*i] = fig.add_axes([float(i)/3., 0.25, 1./3., 0.75]) 
        axes[2*i+1] = fig.add_axes([float(i)/3., 0., 1./3., 0.25]) 
    for i in range(N_reg):
        region = np.loadtxt(input_root + str(i) + dv_input)
        if i==0:
            axes(seps, np.power(seps, 2)*region[:9,1], label="Mono JK Realizations",
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


input_root = sys.argv[1]
dv_output = sys.argv[2]
N_reg = int(sys.argv[3])

jk_corr_funcs_test(input_root, dv_output, "eBOSS LRG v7_2 PIP+ANG", references=["../data/dv_eBOSS_LRG_comb_v7_2_pip_ang_ab.dat"],
                   reference_labels=["Reference"], N_reg=N_reg, include_full=True, output_path="../output/jk_test_")
