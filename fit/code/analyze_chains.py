# Analyze chains from Abacus emulator fits
# v1.0.5, 2023-09-03 - Added option for gamma_l only plot for defence

from collections import OrderedDict
from getdist.chains import ParamError
from getdist.mcsamples import loadMCSamples
import getdist.plots as gdplt
import matplotlib as mpl
import numpy as np
import os
import seaborn as sns
import sys
import yaml


print("Starting...")
mpl.rcParams['figure.dpi'] = 300

def omegam(z, omegam0=0.315):
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
    l = 0.83 + 0.19 / (omegam0**0.8)
    beta = 0.98 + 0.01 / (omegam0**1.2)
    gamma = 0.64 + 0.09 / (omegam0**0.4)
    print("lambda:", l, "beta:", beta, "gamma:", gamma)
    return l * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


# Load Params
input_file = sys.argv[1]
with open(input_file) as file:
    input_dict = yaml.load(file, Loader=yaml.FullLoader)

if input_dict["z_data"] == 0.7:
    fsig8_ref = 0.4638823655118819 # z=0.7
else:
    fsig8_ref = fsigma8_approximate(input_dict["z_data"])

print("fsig8_ref:", fsig8_ref)

default_vals = OrderedDict([("omegam", 0.3159),
                ("ombh2", 0.0224),
                ("sigma8_aem", 0.8111),
                ("H0", 67.36),
                ("ns", 0.965),
                ("logM_sat", 14.93),
                ("alpha", 0.43),
                ("logM_cut", 11.62),
                ("sigma_logM", 0.81),
                ("v_bc", 0.),
                ("v_bs", 1.),
                ("c_vir", 1.),
                ("gamma_n", 1.),
                ("fmax", 1.),
                ("gamma_l", 1.)])

samples = []
samples_reduced = []
for a, path in enumerate(input_dict["paths"]):
    folder, name = os.path.split(os.path.abspath(path))
    labels = input_dict["params"].copy()
    if name == "planck_test":
        sample = loadMCSamples(path, no_cache=True)
    else:
        if isinstance(input_dict["ignore_rows"], list):
            sample = loadMCSamples(path, settings={"ignore_rows": input_dict["ignore_rows"][a]}, no_cache=True)
        else:
            sample = loadMCSamples(path, settings={"ignore_rows": input_dict["ignore_rows"]}, no_cache=True)
        p = sample.getParams()
        min_aem_chi2 = np.min(p.chi2__emulator_likelihood)
        min_index = np.argmin(p.chi2__emulator_likelihood)
        min_tot_chi2 = p.chi2[min_index]
        attributes = dir(p)
        min_params = []
        for param in default_vals.keys():
            if param in attributes:
                min_params.append(getattr(p, param)[min_index])
            else:
                min_params.append(default_vals[param])

        print("Minimum emulator chi^2:", min_aem_chi2)
        print("Minimum Total chi^2:", min_tot_chi2)
        try:
            if not input_dict["constraint_test"]:
                print("Best fit fsig8_comp:", p.fsigma8_comp[min_index])
        except KeyError:
            print("Best fit fsig8_comp:", p.fsigma8_comp[min_index])
        print("Best fit params:", min_params)

        while True:
            try:
                means = sample.mean(labels)
                cov = sample.cov(labels)

            except ParamError as err:
                bad_label = err.args[0].split()[-1]
                print(err.args)
                print("KeyError: Removing %s from labels." % (bad_label))
                labels.remove(bad_label)
                continue

            break

        print(name)
        for i in range(len(labels)):
            print(labels[i], means[i], cov[i,i]**0.5)
        print("\n")

        samples_reduced.append(sample)

    samples.append(sample)


if input_dict["colour_type"] == "sequential":
    contour_colors = sns.color_palette("viridis")

else:
    contour_colors = ["tab:orange", "tab:green", "tab:red", "tab:purple", "tab:brown", "tab:pink", "tab:gray", "tab:olive"]

if input_dict["baseline_index"] != None:
    contour_colors.insert(input_dict["baseline_index"], "tab:blue")

reduced_fig_width = input_dict["reduced_fig_width"]
fig_width = input_dict["fig_width"]

try:
    markers = {**input_dict["markers"]}

except KeyError:
    markers = {"fsigma8_comp": fsig8_ref, "gamma_n": 1., "gamma_l": 1.}

try:
    if input_dict["test_hod"] is not None:
        markers = {"omegam": 0.314153, "ombh2": 0.02222, "sigma8": 0.83,
                   "H0": 67.26, "ns": 0.9652}
        hod_models = np.loadtxt("/home/mj3chapm/P2/fit/data/test_clustering/"
                                "HOD_parameters_test.dat")
        hod_params = hod_models[input_dict["test_hod"], :]
        markers["logM_sat"] = np.log10(hod_params[0])
        markers["alpha"] = hod_params[1]
        markers["logM_cut"] = np.log10(hod_params[2])
        markers["sigma_logM"] = hod_params[3]
        markers["fmax"] = hod_params[8]
        markers["v_bc"] = hod_params[4]
        markers["v_bs"] = hod_params[5]
        markers["c_vir"] = hod_params[6]
        markers["gamma_l"] = hod_params[9]
        markers["gamma_n"] = hod_params[7]
        markers["fsigma8_comp"] = hod_params[9] * fsig8_ref


except KeyError:
    print("Not a test HOD run.")

try:
    param_limits = input_dict["param_limits"]

except KeyError:
    param_limits = {}


if input_dict["reduced_tp"]:
    gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=reduced_fig_width)
    gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
    gdplot.settings.scaling_factor = 1.25
    gdplot.triangle_plot(samples_reduced, input_dict["params_reduced"], filled=True, markers=markers, legend_labels=input_dict["legend_labels"],
                         contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits)
    gdplot.export(input_dict["output_root"] + "_reduced_tp.png")


if input_dict["fsig8_plot"]:
    # gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=reduced_fig_width)
    # gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
    # gdplot.triangle_plot(samples_reduced, ["v_bc", "fsigma8_comp"], filled=True, markers={"fsigma8_comp": fsig8_ref}, legend_labels=input_dict["legend_labels"],
    #                      contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits)
    # gdplot.export(input_dict["output_root"] + "_fsig8_tp.png")

    gdplot = gdplt.get_single_plotter(width_inch=reduced_fig_width)
    gdplot.plot_1d(samples_reduced, "gamma_l",
                   colors=contour_colors[:len(input_dict["legend_labels"])])
    gdplot.add_legend(input_dict["legend_labels"], legend_loc='upper right')
    gdplot.add_x_marker(1.)
    gdplot.export(input_dict["output_root"] + "_gamma_l_single.png")


if input_dict["gamma_l_plot"]:
    gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=reduced_fig_width)
    gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
    gdplot.triangle_plot(samples_reduced, ["v_bc", "gamma_l"], filled=True, markers={"gamma_l": 1.}, legend_labels=input_dict["legend_labels"],
                         contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits)
    gdplot.export(input_dict["output_root"] + "_gamma_l_v_bc_tp.png")

    gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=reduced_fig_width)
    gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
    gdplot.triangle_plot(samples_reduced, ["gamma_l", "fsigma8_comp"], filled=True, markers={"gamma_l": 1., "fsigma8_comp": fsig8_ref}, legend_labels=input_dict["legend_labels"],
                         contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits)
    gdplot.export(input_dict["output_root"] + "_gamma_l_fsig8_tp.png")

    gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=reduced_fig_width)
    gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
    gdplot.triangle_plot(samples_reduced, ["gamma_n", "fsigma8_comp"], filled=True, markers={"gamma_n": 1., "fsigma8_comp": fsig8_ref}, legend_labels=input_dict["legend_labels"],
                         contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits)
    gdplot.export(input_dict["output_root"] + "_gamma_n_fsig8_tp.png")


if input_dict["gammas_plot"]:
    gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=reduced_fig_width)
    gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
    gdplot.triangle_plot(samples_reduced, ["v_bc", "gamma_l", "gamma_n"], filled=True, markers={"gamma_n": 1., "gamma_l": 1.}, legend_labels=input_dict["legend_labels"],
                         contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits, plot_args={"dpi": 600})
    gdplot.export(input_dict["output_root"] + "_gammas_tp.png")


try:
    if "planck_test" in input_dict["paths"][1]:
        contour_colors.insert(1, "tab:cyan")
        input_dict["legend_labels"].insert(1, "Planck18")

except IndexError:
    print("Only one chain, so assuming no Planck")


if input_dict["full_tp"]:
    gdplot_settings = gdplt.GetDistPlotSettings(fig_width_inch=fig_width)
    gdplot = gdplt.get_subplot_plotter(settings=gdplot_settings)
# gdplot = gdplt.get_subplot_plotter()
    gdplot.triangle_plot(samples, input_dict["params"], filled=True, markers=markers, legend_labels=input_dict["legend_labels"],
                         contour_colors=contour_colors[:len(input_dict["legend_labels"])], param_limits=param_limits)
    gdplot.export(input_dict["output_root"] + "_full_tp.png")



# Change Log
# v1.0.4, 2023-06-26 - Added option for constraint test
# v1.0.3, 2023-06-07 - Increased resolution, switched to png
# v1.0.2, 2022-12-01 - Added markers for test_hod runs
# v1.0.1, 2022-06-29 - Updated 'aemulus' to emulator
# v1.0.0, 2022-01-14 - Copied from RSD/fit/code/analyze_afmax.py
