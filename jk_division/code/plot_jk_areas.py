# Imports
import matplotlib.pyplot as plt
import numpy as np
import sys
from seaborn import color_palette
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.io import ascii as ap_ascii
from pathlib import Path


data_file_ngc = sys.argv[1]
rand_file_ngc = sys.argv[2]
data_file_sgc = sys.argv[3]
rand_file_sgc = sys.argv[4]
output_dir = sys.argv[5]
output_root = sys.argv[6]
N_regions = int(sys.argv[7])
load_regions = sys.argv[8]


def generate_rows(min_dec, max_dec, square_side, buffer_array=np.array([0.5, 0.5])):
    dec_range = max_dec - min_dec
    N_rows = np.ceil(dec_range / square_side)
    row_buffer = square_side * N_rows - dec_range
    buffer_array = row_buffer*buffer_array
    print("N_rows+1", N_rows+1, "int(N_rows+1)", int(N_rows+1))
    row_edges = np.linspace(min_dec - buffer_array[0], max_dec + buffer_array[1], int(N_rows + 1))
        
    return N_rows, row_edges


def divide_catalog(rand_cat, N_regions=200, area_scaling=1., row_start="Top", start_position=-2.):
    square_area = 6000. / N_regions * area_scaling
    square_side = np.sqrt(square_area)
    
    if row_start == "Top":
        buffer_array = np.array([1., 0.])
        N_rows, row_edges = generate_rows(np.min(rand_cat['DEC']), np.max(rand_cat['DEC']), square_side,
                                          buffer_array=buffer_array)
     
    if row_start == "Mid":
        buffer_array_low = np.array([1., 0.])
        N_rows_low, row_edges_low = generate_rows(np.min(rand_cat['DEC']), start_position, square_side,
                                          buffer_array=buffer_array_low)
        buffer_array_high = np.array([0., 1.])
        N_rows_high, row_edges_high = generate_rows(start_position, np.max(rand_cat['DEC']), square_side,
                                          buffer_array=buffer_array_high)
        N_rows = N_rows_low + N_rows_high
        row_edges = np.concatenate((row_edges_low[:-1], row_edges_high))

    square_widths = []
    N_squares = []
    square_edges_list = []
    for i in range(int(N_rows)):
        row_centre = (row_edges[i+1] + row_edges[i]) / 2.
        square_width = square_side / np.cos(row_centre * np.pi / 180.)
        square_widths.append(square_width)
        row_mask_1 = rand_cat['DEC'] >= row_edges[i]
        rand_row = rand_cat[row_mask_1]
        row_mask_2 = rand_row['DEC'] < row_edges[i+1]
        rand_row = rand_row[row_mask_2]
        ra_max = np.max([ra-360. if ra>300. else ra for ra in rand_row['RA']])
        ra_min = np.min([ra-360. if ra>300. else ra for ra in rand_row['RA']])
        ra_range = ra_max - ra_min
        N_row_squares = np.ceil(ra_range / square_width)
        N_squares.append(N_row_squares)
        ra_buffer = (square_width * N_row_squares - ra_range) / 2.
        if i==0:
            print("N_row_squares+1", N_row_squares+1, "int(N_row_squares+1)", int(N_row_squares+1))
        square_edges_list.append(np.linspace(ra_min - ra_buffer, ra_max + ra_buffer, int(N_row_squares + 1)))
        
    return N_rows, row_edges, N_squares, square_edges_list


def assign_region(catalog, N_rows, row_edges, N_squares, square_edges_list, cadence=500000,
                  initial_index=0, removed_regions=[]):
    
    catalog['JK_REGION'] = -1
    for i in range(len(catalog)):
        for j in range(int(N_rows)):
            if catalog['DEC'][i] > row_edges[j] and catalog['DEC'][i] <= row_edges[j+1]:
                for k in range(int(N_squares[j])):
                    if catalog['RA'][i] > 300.:
                        if ((catalog['RA'][i] - 360.) > square_edges_list[j][k] and
                            (catalog['RA'][i] - 360.) <= square_edges_list[j][k+1]):
                            
                            if (int(np.sum(N_squares[:j])) + k + initial_index) not in removed_regions:
                                if catalog['JK_REGION'][i] != -1:
                                    print('Object reassigned from', catalog['JK_REGION'][i], "to",
                                          int(np.sum(N_squares[:j])) + k + initial_index)
                                catalog['JK_REGION'][i] = int(np.sum(N_squares[:j])) + k + initial_index
                                break
                                
                            else:
                                break
                            
                    else:
                        if (catalog['RA'][i] > square_edges_list[j][k] and 
                            catalog['RA'][i] <= square_edges_list[j][k+1]):

                            if (int(np.sum(N_squares[:j])) + k + initial_index) not in removed_regions:
                                if catalog['JK_REGION'][i] != -1:
                                    print('Object reassigned from', catalog['JK_REGION'][i], "to",
                                          int(np.sum(N_squares[:j])) + k + initial_index)
                                catalog['JK_REGION'][i] = int(np.sum(N_squares[:j])) + k + initial_index
                                break
                                
                            else:
                                break
                break
        
        if (i % cadence)==0:
            print("Assigned %i" % i)
        
    return catalog


def filled_areas_division(data_cat, rand_cat, N_regions=200, row_start="Top", area_scaling=1., occupation_scaling=0.8,
                          mark_regions=False, initial_index=0):
    N_rows, row_edges, N_squares, square_edges_list = divide_catalog(rand_cat, N_regions=N_regions,
                                                                 area_scaling=area_scaling, row_start=row_start)
    rand_cat = assign_region(rand_cat, N_rows, row_edges, N_squares, square_edges_list,
                             initial_index=initial_index)
    unassigned_mask = (rand_cat["JK_REGION"]==1)
    print("Unassigned Objects:", len(rand_cat[unassigned_mask]))

    N = int(np.sum(N_squares))
    binned_regions, edges = np.histogram(rand_cat['JK_REGION'], N, weights=(rand_cat["WEIGHT_SYSTOT"]*
                                         rand_cat["WEIGHT_CP"]*rand_cat["WEIGHT_NOZ"]*rand_cat["WEIGHT_FKP"]))
    
    max_occupation = np.max(binned_regions)
    
    print("Max Occupation:", max_occupation)
    
    removed_regions = []
    accepted_regions = 0
    for i in range(len(binned_regions)):
        if binned_regions[i] >= max_occupation*occupation_scaling:
            accepted_regions += 1

        else:
            binned_regions[i] = 0
            removed_regions.append(i + initial_index)

    for obj in rand_cat:
        if obj["JK_REGION"] in removed_regions:
            obj["JK_REGION"] = -1
            
    print("Number of Accepted Regions:", accepted_regions)
    print("Removed Regions:", removed_regions)
    print("Desired Number of Regions:", N_regions)
        
    return N_rows, row_edges, N_squares, square_edges_list, rand_cat, accepted_regions, removed_regions


def plot_catalog(cat_ngc, cat_sgc, ngc_region_params, sgc_region_params, removed_regions=[],
                 mark_regions=False, output_filename="../output/divided_catalog.png",
                 initial_index=0):

    N_rows_ngc, row_edges_ngc, N_squares_ngc, square_edges_list_ngc = ngc_region_params
    N_rows_sgc, row_edges_sgc, N_squares_sgc, square_edges_list_sgc = sgc_region_params

    N_rows_list = [N_rows_ngc, N_rows_sgc]
    row_edges_list = [row_edges_ngc, row_edges_sgc]
    N_squares_list = [N_squares_ngc, N_squares_sgc]
    square_edges_list_list = [square_edges_list_ngc, square_edges_list_sgc]
    cats = [cat_ngc, cat_sgc]

    # fig, axes = plt.subplots(figsize=(4.75, 6.75))

    # for a in (0, 1):
    #     axes[a].scatter([ra-360. if ra>300. else ra for ra in cats[a]['RA']], cats[a]['DEC'], s=0.2)
    #     axes[a].set_xlabel("RA [deg]")
    #     axes[a].set_ylabel("DEC [deg]")

    #     for i in range(int(np.sum(N_squares_list[a]))):
    #         if i + initial_index not in removed_regions:
    #             for j in range(int(N_rows_list[a])):
    #                 if i < np.sum(N_squares_list[a][:j+1]):                         
#                         square_index = int(i - np.sum(N_squares_list[a][:j]))
#                         axes[a].plot([square_edges_list_list[a][j][square_index], square_edges_list_list[a][j][square_index]],
#                                  [row_edges_list[a][j], row_edges_list[a][j+1]], 'k-')
#                         axes[a].plot([square_edges_list_list[a][j][square_index], square_edges_list_list[a][j][square_index+1]],
#                                  [row_edges_list[a][j+1], row_edges_list[a][j+1]], 'k-')
#                         axes[a].plot([square_edges_list_list[a][j][square_index+1], square_edges_list_list[a][j][square_index+1]],
#                                  [row_edges_list[a][j], row_edges_list[a][j+1]], 'k-')
#                         axes[a].plot([square_edges_list_list[a][j][square_index], square_edges_list_list[a][j][square_index+1]],
#                                  [row_edges_list[a][j], row_edges_list[a][j]], 'k-')
#                         if mark_regions:
#                             axes[a].plot((square_edges_list_list[a][j][square_index]+square_edges_list_list[a][j][square_index+1])/2.,
#                                      (row_edges_list[a][j]+row_edges_list[a][j+1])/2., color=color_palette('colorblind')[2],
#                                      marker='o')
#                         break

    fig_width = 10.
    plt.figure(figsize=(fig_width, fig_width / 2.), dpi=300)

    for a in (0, 1):
        plt.scatter([ra-360. if ra>300. else ra for ra in cats[a]['RA']], cats[a]['DEC'], s=0.01, color=color_palette('muted')[a])
        plt.xlabel("RA [deg]")
        plt.ylabel("DEC [deg]")

        for i in range(int(np.sum(N_squares_list[a]))):
            if i + initial_index not in removed_regions:
                for j in range(int(N_rows_list[a])):
                    if i < np.sum(N_squares_list[a][:j+1]):
                        square_index = int(i - np.sum(N_squares_list[a][:j]))
                        plt.plot([square_edges_list_list[a][j][square_index], square_edges_list_list[a][j][square_index]],
                                 [row_edges_list[a][j], row_edges_list[a][j+1]], 'k-')
                        plt.plot([square_edges_list_list[a][j][square_index], square_edges_list_list[a][j][square_index+1]],
                                 [row_edges_list[a][j+1], row_edges_list[a][j+1]], 'k-')
                        plt.plot([square_edges_list_list[a][j][square_index+1], square_edges_list_list[a][j][square_index+1]],
                                 [row_edges_list[a][j], row_edges_list[a][j+1]], 'k-')
                        plt.plot([square_edges_list_list[a][j][square_index], square_edges_list_list[a][j][square_index+1]],
                                 [row_edges_list[a][j], row_edges_list[a][j]], 'k-')
                        if mark_regions:
                            plt.plot((square_edges_list_list[a][j][square_index]+square_edges_list_list[a][j][square_index+1])/2.,
                                     (row_edges_list[a][j]+row_edges_list[a][j+1])/2., color=color_palette('colorblind')[2],
                                     marker='o')
                        break

        initial_index = int(np.sum(N_squares_list[a]))

    plt.savefig(output_filename)
    

data_cat_ngc = Table.read(data_file_ngc, format="ascii")
data_cat_sgc = Table.read(data_file_sgc, format="ascii")

if load_regions == "True":
    N_rows_array = np.loadtxt(output_dir + output_root + "_N_rows.dat", skiprows=1, delimiter="\t")
    N_rows_ngc, N_rows_sgc = list(N_rows_array)
    
    row_edges_ngc = np.loadtxt(output_dir + output_root + "_row_edges_ngc.dat", delimiter="\t")
    row_edges_sgc = np.loadtxt(output_dir + output_root + "_row_edges_sgc.dat", delimiter="\t")

    N_squares_ngc = list(np.loadtxt(output_dir + output_root + "_N_squares_ngc.dat", delimiter="\t"))
    N_squares_sgc = list(np.loadtxt(output_dir + output_root + "_N_squares_sgc.dat", delimiter="\t"))

    removed_regions = np.loadtxt(output_dir + output_root + "_removed_regions.dat", delimiter="\t")

    square_edges_list_ngc = []
    square_edges_list_sgc = []

    for i in range(int(N_rows_ngc)):
        square_edges_list_ngc.append(np.loadtxt("%ssquare_edges/%s_ngc.%i.dat" % (output_dir, output_root, i), delimiter="\t"))

    for i in range(int(N_rows_sgc)):
        square_edges_list_sgc.append(np.loadtxt("%ssquare_edges/%s_sgc.%i.dat" % (output_dir, output_root, i), delimiter="\t"))

else:
    rand_cat_ngc = Table.read(rand_file_ngc, format="ascii")
    rand_cat_sgc = Table.read(rand_file_sgc, format="ascii")

    kwargs = {"N_regions": N_regions, "occupation_scaling": 0.8,
               "area_scaling": 0.59049}

    print("Dividing NGC")
    (N_rows_ngc, row_edges_ngc, N_squares_ngc, square_edges_list_ngc, rand_cat_ngc,
     accepted_regions_ngc, removed_regions_ngc) = filled_areas_division(data_cat_ngc, rand_cat_ngc, **kwargs,
                                                                       row_start="Top")
    total_regions_ngc = accepted_regions_ngc + len(removed_regions_ngc)

    print("Dividing SGC")
    (N_rows_sgc, row_edges_sgc, N_squares_sgc, square_edges_list_sgc, rand_cat_sgc,
     accepted_regions_sgc, removed_regions_sgc) = filled_areas_division(data_cat_sgc, rand_cat_sgc, **kwargs,
                                                                       initial_index=total_regions_ngc,
                                                                       row_start="Mid")

    accepted_regions = accepted_regions_ngc + accepted_regions_sgc
    removed_regions = removed_regions_ngc + removed_regions_sgc

    print("Number of Accepted Regions Total:", accepted_regions)
    print("Removed Regions Total:", removed_regions)
    print("Desired Number of Regions:", N_regions)

    N = accepted_regions + len(removed_regions)
    print("N+2", N+2, "int(N+2)", int(N+2))
    bins = np.linspace(-1.5, N-0.5, int(N+2))
    while accepted_regions > N_regions:
        weight_table = np.concatenate((rand_cat_ngc["WEIGHT_SYSTOT"]*rand_cat_ngc["WEIGHT_CP"]*
                               rand_cat_ngc["WEIGHT_NOZ"]*rand_cat_ngc["WEIGHT_FKP"],
                               rand_cat_sgc["WEIGHT_SYSTOT"]*rand_cat_sgc["WEIGHT_CP"]*
                               rand_cat_sgc["WEIGHT_NOZ"]*rand_cat_sgc["WEIGHT_FKP"]))
        binned_regions, edges = np.histogram(np.concatenate((rand_cat_ngc['JK_REGION'],
                                                             rand_cat_sgc['JK_REGION'])), bins=bins,
                                             weights = weight_table)
        
        smallest = np.min(binned_regions[binned_regions>0.][1:])
        smallest_index = np.where(binned_regions==smallest)[0][0]
        # Need to subtract one from this index when comparing to regions since there is the unassigned bin
        smallest_region = smallest_index - 1
        smallest_mask_ngc = rand_cat_ngc['JK_REGION'] == smallest_region
        smallest_objects_ngc = rand_cat_ngc[smallest_mask_ngc]
        smallest_mask_sgc = rand_cat_sgc['JK_REGION'] == smallest_region
        smallest_objects_sgc = rand_cat_sgc[smallest_mask_sgc]
        region_occupation = (np.sum(smallest_objects_ngc["WEIGHT_SYSTOT"]*smallest_objects_ngc["WEIGHT_CP"]*
                                    smallest_objects_ngc["WEIGHT_NOZ"]*smallest_objects_ngc["WEIGHT_FKP"]) +
                             np.sum(smallest_objects_sgc["WEIGHT_SYSTOT"]*smallest_objects_sgc["WEIGHT_CP"]*
                                    smallest_objects_sgc["WEIGHT_NOZ"]*smallest_objects_sgc["WEIGHT_FKP"]))
        print("Lowest occupation region:", smallest_region)
        print("Weighted occupation of region:", smallest)
        print("Occupation of removed region:", binned_regions[smallest_index])
        print("Occupation by index:", region_occupation)
        binned_regions[smallest_index] = 0
        removed_regions.append(smallest_region)
        
        if smallest_region < total_regions_ngc:
            for obj in rand_cat_ngc:
                if obj["JK_REGION"] == smallest_region:
                    obj["JK_REGION"] = -1
        
        else:
            for obj in rand_cat_sgc:
                if obj["JK_REGION"] == smallest_region:
                    obj["JK_REGION"] = -1
                    
                
        accepted_regions -= 1

    N_rows_array = [N_rows_ngc, N_rows_sgc]
    np.savetxt(output_dir + output_root + "_N_rows.dat", N_rows_array, header="N_rows_ngc\tN_rows_sgc", delimiter="\t")

    np.savetxt(output_dir + output_root + "_row_edges_ngc.dat", row_edges_ngc, delimiter="\t")
    np.savetxt(output_dir + output_root + "_row_edges_sgc.dat", row_edges_sgc, delimiter="\t")

    np.savetxt(output_dir + output_root + "_N_squares_ngc.dat", N_squares_ngc, delimiter="\t")
    np.savetxt(output_dir + output_root + "_N_squares_sgc.dat", N_squares_sgc, delimiter="\t")

    np.savetxt(output_dir + output_root + "_removed_regions.dat", removed_regions, delimiter="\t")

    Path(output_dir + "square_edges").mkdir(exist_ok=True)
    for i in range(len(square_edges_list_ngc)):
        np.savetxt("%ssquare_edges/%s_ngc.%i.dat" % (output_dir, output_root, i), square_edges_list_ngc[i], delimiter="\t")

    for i in range(len(square_edges_list_sgc)):
        np.savetxt("%ssquare_edges/%s_sgc.%i.dat" % (output_dir, output_root, i), square_edges_list_sgc[i], delimiter="\t")

ngc_region_params = [N_rows_ngc, row_edges_ngc, N_squares_ngc, square_edges_list_ngc]
sgc_region_params = [N_rows_sgc, row_edges_sgc, N_squares_sgc, square_edges_list_sgc]

plot_catalog(data_cat_ngc, data_cat_sgc, ngc_region_params, sgc_region_params,
             removed_regions=removed_regions, mark_regions=False,
             output_filename=output_dir + output_root + "_jk_areas.png")
