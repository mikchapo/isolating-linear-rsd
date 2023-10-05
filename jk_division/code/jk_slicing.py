# Imports
import matplotlib.pyplot as plt
import numpy as np
import sys
from seaborn import color_palette
from astropy.table import Table, vstack
from astropy.io import fits
from astropy.io import ascii as ap_ascii
from pathlib import Path

# Required input. Call script as:
# python /path/to/jk_areas.py data_file_ngc rand_file_ngc data_file_sgc rand_file_sgc output_root N_reg area_scaling
data_file_ngc = sys.argv[1]         # Path to file containing the data NGC catalogue in ascii format. See below for column names
rand_file_ngc = sys.argv[2]
data_file_sgc = sys.argv[3]
rand_file_sgc = sys.argv[4]
output_root = sys.argv[5]           # The root of the path where output is stored. For example, the divided data NGC catalogue is written
                                    # to "%s_NGC_jk_%i.dat" % (output_root, N_reg)     
N_reg = sys.argv[6]                 # The number of regions to have in the final division. 200 in fiducial analysis
area_scaling = float(sys.argv[7])   # Initial size scaling to apply to the regions to reduce the number of steps to find the correct size.
                                    # Default should be 1 if unknown. For eBOSS LRG area_scaling=0.59049

# Necessary column names:
# RA DEC WEIGHT_SYSTOT WEIGHT_CP WEIGHT_NOZ WEIGHT_FKP


# Assign functions
def generate_rows(min_dec, max_dec, square_side, buffer_array=np.array([0.5, 0.5]), verbose=True):
    dec_range = max_dec - min_dec
    N_rows = np.ceil(dec_range / square_side)
    row_buffer = square_side * N_rows - dec_range
    buffer_array = row_buffer*buffer_array
    print("N_rows+1", N_rows+1, "int(N_rows+1)", int(N_rows+1))
    row_edges = np.linspace(min_dec - buffer_array[0], max_dec + buffer_array[1], int(N_rows + 1))
    
    if verbose:
        print("Rows:", N_rows)
        print("Row Buffer:", row_buffer)
        print("Row Edges:", row_edges)
        
    return N_rows, row_edges

def divide_catalog(rand_cat, N_regions=200, area_scaling=1., row_start="Top", start_position=-2., verbose=True):
    square_area = 6000. / N_regions * area_scaling
    square_side = np.sqrt(square_area)

    if verbose:
        print("Area:", square_area)
        print("Side Length:", square_side)
    
    if row_start == "Top":
        buffer_array = np.array([1., 0.])
        N_rows, row_edges = generate_rows(np.min(rand_cat['DEC']), np.max(rand_cat['DEC']), square_side,
                                          buffer_array=buffer_array, verbose=verbose)
     
    if row_start == "Mid":
        buffer_array_low = np.array([1., 0.])
        N_rows_low, row_edges_low = generate_rows(np.min(rand_cat['DEC']), start_position, square_side,
                                          buffer_array=buffer_array_low, verbose=verbose)
        buffer_array_high = np.array([0., 1.])
        N_rows_high, row_edges_high = generate_rows(start_position, np.max(rand_cat['DEC']), square_side,
                                          buffer_array=buffer_array_high, verbose=verbose)
        N_rows = N_rows_low + N_rows_high
        row_edges = np.concatenate((row_edges_low[:-1], row_edges_high))

    square_widths = []
    N_squares = []
    square_edges_list = []
    row_ends = []
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
        row_ends.append([ra_min - ra_buffer, ra_max + ra_buffer])
        if i==0:
            print("N_row_squares+1", N_row_squares+1, "int(N_row_squares+1)", int(N_row_squares+1))
        square_edges_list.append(np.linspace(ra_min - ra_buffer, ra_max + ra_buffer, int(N_row_squares + 1)))

    if verbose:
        print("Square Widths:", square_widths)
        print("N Squares:", N_squares)
        print("Row Ends:", row_ends)
        print("Square Edges:", square_edges_list)
        
    return N_rows, row_edges, N_squares, square_edges_list, row_ends


def plot_catalog(data_cat, rand_cat, N_rows, row_edges, N_squares, square_edges_list, row_ends, removed_regions=[],
                 save_plot=False, mark_regions=False, output_filename="../output/divided_catalog.png",
                 initial_index=0):
    plt.figure(figsize=(16,12))
    plt.scatter([ra-360. if ra>300. else ra for ra in rand_cat['RA']], rand_cat['DEC'],
                color=color_palette('colorblind')[1], s=0.02, label="Rand")
    plt.scatter([ra-360. if ra>300. else ra for ra in data_cat['RA']], data_cat['DEC'],
                color=color_palette('colorblind')[0], s=0.1, label="Data")
#     for i in range(int(N_rows + 1)):
#         if i==0:
#             row_end = row_ends[i]

#         elif i==N_rows_ngc:
#             row_end = row_ends[i-1]

#         else:
#             row_end = [np.min([row_ends[i-1][0], row_ends[i][0]]), np.max([row_ends[i-1][1], row_ends[i][1]])]

#         plt.plot(row_end, [row_edges[i], row_edges[i]], 'k-')

#     for i in range(int(N_rows)):
#         for j in range(int(N_squares[i]) + 1):
#             plt.plot([square_edges_list[i][j], square_edges_list[i][j]], [row_edges[i], row_edges[i+1]], 'k-')

    for i in range(int(np.sum(N_squares))):
        if i + initial_index not in removed_regions:
            for j in range(int(N_rows)):
                if i < np.sum(N_squares[:j+1]):
                    if i == 355:
                        square_index = int(i - np.sum(N_squares[:j]))
                        plt.plot([square_edges_list[j][square_index], square_edges_list[j][square_index]],
                                 [row_edges[j], row_edges[j+1]], linestyle='-', color=color_palette('colorblind')[3])
                        plt.plot([square_edges_list[j][square_index], square_edges_list[j][square_index+1]],
                                 [row_edges[j+1], row_edges[j+1]], linestyle='-', color=color_palette('colorblind')[3])
                        plt.plot([square_edges_list[j][square_index+1], square_edges_list[j][square_index+1]],
                                 [row_edges[j], row_edges[j+1]], linestyle='-', color=color_palette('colorblind')[3])
                        plt.plot([square_edges_list[j][square_index], square_edges_list[j][square_index+1]],
                                 [row_edges[j], row_edges[j]], linestyle='-', color=color_palette('colorblind')[3])
                        if mark_regions:
                            plt.plot((square_edges_list[j][square_index]+square_edges_list[j][square_index+1])/2.,
                                     (row_edges[j]+row_edges[j+1])/2., color=color_palette('colorblind')[3], marker='o')
                        break
                        
                    else:
                        square_index = int(i - np.sum(N_squares[:j]))
                        plt.plot([square_edges_list[j][square_index], square_edges_list[j][square_index]],
                                 [row_edges[j], row_edges[j+1]], 'k-')
                        plt.plot([square_edges_list[j][square_index], square_edges_list[j][square_index+1]],
                                 [row_edges[j+1], row_edges[j+1]], 'k-')
                        plt.plot([square_edges_list[j][square_index+1], square_edges_list[j][square_index+1]],
                                 [row_edges[j], row_edges[j+1]], 'k-')
                        plt.plot([square_edges_list[j][square_index], square_edges_list[j][square_index+1]],
                                 [row_edges[j], row_edges[j]], 'k-')
                        if mark_regions:
                            plt.plot((square_edges_list[j][square_index]+square_edges_list[j][square_index+1])/2.,
                                     (row_edges[j]+row_edges[j+1])/2., color=color_palette('colorblind')[2],
                                     marker='o')
                        break
    plt.title("Sky distribution of randoms and data with accepted jackknife regions")
    plt.xlabel("RA [deg]")
    plt.ylabel("DEC [deg]")
    plt.legend()

    if save_plot:
        plt.savefig(output_filename)
        
    else:
        plt.show()
        
        
def assign_region(catalog, N_rows, row_edges, N_squares, square_edges_list, verbose=True, cadence=500000,
                  initial_index=0, removed_regions=[]):
    if verbose:
        print("Starting Assignment")
    
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
    
    if verbose:
        print("Assignment Finished")
        
    return catalog
        
        
# def test_division(catalog, N_regions, N_squares, ratio, total_objects):
#     N = int(np.sum(N_squares))
#     binned_regions, edges = np.histogram(catalog['JK_REGION'], N)

#     print("Total Number of Regions:", N)
#     print("Desired Number of Regions:", round(ratio*N_regions, 2))
#     print("Mean Occupation:", round(np.mean(binned_regions), 1))
#     print("Desired Mean:", round(total_objects / N_regions, 1))
#     print("Standard Deviation:", round(np.std(binned_regions), 1))
#     print("Relative Standard Deviation:", round(np.std(binned_regions) / np.mean(binned_regions), 2))


#     plt.figure(figsize=(16,12))
#     plt.hist(binned_regions, 60)
#     plt.show()
    

def filled_areas_division(data_cat, rand_cat, N_regions=200, row_start="Top", area_scaling=1., occupation_scaling=0.8,
                          verbose=False, save_plot=False, mark_regions=False, initial_index=0):
    N_rows, row_edges, N_squares, square_edges_list, row_ends = divide_catalog(rand_cat, N_regions=N_regions,
                                                                 area_scaling=area_scaling, row_start=row_start,
                                                                 verbose=verbose)
    rand_cat = assign_region(rand_cat, N_rows, row_edges, N_squares, square_edges_list, verbose=verbose,
                             initial_index=initial_index)
    unassigned_mask = (rand_cat["JK_REGION"]==1)
    print("Unassigned Objects:", len(rand_cat[unassigned_mask]))

    N = int(np.sum(N_squares))
    binned_regions, edges = np.histogram(rand_cat['JK_REGION'], N, weights=(rand_cat["WEIGHT_SYSTOT"]*
                                         rand_cat["WEIGHT_CP"]*rand_cat["WEIGHT_NOZ"]*rand_cat["WEIGHT_FKP"]))
    if verbose:
        print("Initial division bins:")
        print(binned_regions)
        print("Initial division binning edges:")
        print(edges)
    
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

    if verbose:
        print("Bins after removing low occupation:")
        print(binned_regions)
        print("Edges after removing low occupation:")
        print(edges)
            
    print("Number of Accepted Regions:", accepted_regions)
    print("Removed Regions:", removed_regions)
    print("Desired Number of Regions:", N_regions)
        
    return N_rows, row_edges, N_squares, square_edges_list, row_ends, rand_cat, accepted_regions, removed_regions


def repeated_division(data_cat_ngc, data_cat_sgc, rand_cat_ngc, rand_cat_sgc, N_regions=200, area_scaling=1.,
                      occupation_scaling=0.8, verbose=False, save_plot=False, mark_regions=False, row_start_ngc="Top",
                      row_start_sgc="Mid", output_root="../output/eBOSS_LRG_v7", write=False, output_format="fits"):
    # Sample kwargs
    # kwargs_ngc = {"output_filename": "../output/divided_catalog.png", "row_start": "Top"}
    
    # Updated kwargs for ngc and sgc. Organized to prevent a disconnect betweent he two regions while limiting the
    # number of places params are changed
    kwargs = {"N_regions": N_regions, "occupation_scaling": occupation_scaling, "verbose": verbose,
              "save_plot": save_plot, "mark_regions": mark_regions}
    
    insufficient_regions = True
    while insufficient_regions:
        kwargs["area_scaling"] = area_scaling
        
        print("Dividing NGC")
        (N_rows_ngc, row_edges_ngc, N_squares_ngc, square_edges_list_ngc, row_ends_ngc, rand_cat_ngc,
         accepted_regions_ngc, removed_regions_ngc) = filled_areas_division(data_cat_ngc, rand_cat_ngc, **kwargs,
                                                                           row_start=row_start_ngc)
        total_regions_ngc = accepted_regions_ngc + len(removed_regions_ngc)
        
        print("Dividing SGC")
        (N_rows_sgc, row_edges_sgc, N_squares_sgc, square_edges_list_sgc, row_ends_sgc, rand_cat_sgc,
         accepted_regions_sgc, removed_regions_sgc) = filled_areas_division(data_cat_sgc, rand_cat_sgc, **kwargs,
                                                                           initial_index=total_regions_ngc,
                                                                           row_start=row_start_sgc)
        
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
            
            if verbose:
                plt.figure(figsize=(16,12))
                plt.hist(binned_regions[binned_regions>0.][1:], 20)
                plt.title("Histogram of regions by weighted randoms")
                plt.xlabel("Weighted randoms contained by region")
                plt.ylabel("Number of regions")
                plt.show()
                print("Initial combined bins:")
                print(binned_regions)
                print("Initial combined edges:")
                print(edges)
            
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
            
#             removed_counter = 1
            if smallest_region < total_regions_ngc:
                for obj in rand_cat_ngc:
                    if obj["JK_REGION"] == smallest_region:
                        obj["JK_REGION"] = -1
#                         if (removed_counter % 1000)==0:
#                             print("Removed ", removed_counter)
                            
#                         removed_counter += 1
            
            else:
                for obj in rand_cat_sgc:
                    if obj["JK_REGION"] == smallest_region:
                        obj["JK_REGION"] = -1
                        
#                         if (removed_counter % 5000)==0:
#                             print("Removed ", removed_counter)
                            
#                         removed_counter += 1
                        
#             removed_counter = 0
                    
            accepted_regions -= 1
            
        if accepted_regions == N_regions:
            insufficient_regions = False
            
        else:
            area_scaling = area_scaling*0.9

    # THIS CODE NEW
    # I've tested it in a similar script but not here, so there is the possibility of formatting errors. If there is a bug let me know and I can fix it
    # - Mike    
    N_rows_array = [N_rows_ngc, N_rows_sgc]
    np.savetxt(output_root + "_N_rows.dat", N_rows_array, header="N_rows_ngc\tN_rows_sgc", delimiter="\t")

    np.savetxt(output_root + "_row_edges_ngc.dat", row_edges_ngc, delimiter="\t")
    np.savetxt(output_root + "_row_edges_sgc.dat", row_edges_sgc, delimiter="\t")

    np.savetxt(output_root + "_N_squares_ngc.dat", N_squares_ngc, delimiter="\t")
    np.savetxt(output_root + "_N_squares_sgc.dat", N_squares_sgc, delimiter="\t")

    np.savetxt(output_root + "_removed_regions.dat", removed_regions, delimiter="\t")

    Path(output_root + "_square_edges").mkdir(exist_ok=True)
    for i in range(len(square_edges_list_ngc)):
        np.savetxt("%s_square_edges/square_edges_ngc.%i.dat" % (output_root, i), square_edges_list_ngc[i], delimiter="\t")

    for i in range(len(square_edges_list_sgc)):
        np.savetxt("%s_square_edges/square_edges_sgc.%i.dat" % (output_root, i), square_edges_list_sgc[i], delimiter="\t")

    print("Assigning data NGC")
    data_cat_ngc = assign_region(data_cat_ngc, N_rows_ngc, row_edges_ngc, N_squares_ngc, square_edges_list_ngc,
                                 verbose=verbose, removed_regions=removed_regions)
    
    print("Assigning data SGC")
    data_cat_sgc = assign_region(data_cat_sgc, N_rows_sgc, row_edges_sgc, N_squares_sgc, square_edges_list_sgc,
                                 verbose=verbose, removed_regions=removed_regions, initial_index=total_regions_ngc)
    
    # END OF NEW CODE

    print("N+1", N+1, "int(N+1)", int(N+1))
    plot_bins = np.linspace(-0.5, N-0.5, int(N+1))
    
    plt.figure(figsize=(16,12))
    plt.hist(rand_cat_ngc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per accepted region in NGC")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    plt.figure(figsize=(16,12))
    plt.hist(rand_cat_sgc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per accepted region in SGC")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    plt.figure(figsize=(16,12))
    plt.hist(data_cat_ngc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per aceepted region in NGC")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    plt.figure(figsize=(16,12))
    plt.hist(data_cat_sgc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per accepted region in SGC")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    total_rand_objects = np.sum(weight_table)
    print("\n", "Overall Statistics:")
    print("Final Number of Regions:", accepted_regions)
    print("Desired Number of Regions:", N_regions)
    print("Removed Regions:", removed_regions)
    print("Number of Removed Regions:", len(removed_regions))
    
    print("\n", "Randoms Statistics:")
    print("Mean Occupation:", round(np.sum(binned_regions[1:])/accepted_regions, 1))
    print("Desired Mean:", round(total_rand_objects / N_regions, 1))
    print("Standard Deviation:", round(np.std(binned_regions[binned_regions>0.][1:]), 1))
    print("Relative Standard Deviation:", round(np.std(binned_regions[binned_regions>0.][1:]) /
                                                np.mean(binned_regions[binned_regions>0.][1:]), 2))


    plt.figure(figsize=(16,12))
    plt.hist(binned_regions[binned_regions>0.][1:], 20)
    plt.title("Histogram of regions by weighted randoms")
    plt.xlabel("Weighted randoms contained by region")
    plt.ylabel("Number of regions")
    if save_plot:
        filename = output_root + "_rand_region_hist.png"
        plt.savefig(filename)
        
    else:
        plt.show()
    
    data_weight_table = np.concatenate((data_cat_ngc["WEIGHT_SYSTOT"]*data_cat_ngc["WEIGHT_CP"]*
                           data_cat_ngc["WEIGHT_NOZ"]*data_cat_ngc["WEIGHT_FKP"],
                           data_cat_sgc["WEIGHT_SYSTOT"]*data_cat_sgc["WEIGHT_CP"]*
                           data_cat_sgc["WEIGHT_NOZ"]*data_cat_sgc["WEIGHT_FKP"]))
    data_binned_regions, data_edges = np.histogram(np.concatenate((data_cat_ngc['JK_REGION'],
                                                         data_cat_sgc['JK_REGION'])), bins=bins,
                                         weights = data_weight_table)
    total_data_objects = np.sum(data_weight_table)
    print("\n", "Data Statistics:")
    print("Mean Occupation:", round(np.sum(data_binned_regions[1:])/accepted_regions, 1))
    print("Desired Mean:", round(total_data_objects / N_regions, 1))
    print("Standard Deviation:", round(np.std(data_binned_regions[data_binned_regions>0.][1:]), 1))
    print("Relative Standard Deviation:", round(np.std(data_binned_regions[data_binned_regions>0.][1:]) /
                                                np.mean(data_binned_regions[data_binned_regions>0.][1:]), 2))


    plt.figure(figsize=(16,12))
    plt.hist(data_binned_regions[data_binned_regions>0.][1:], 20)
    plt.title("Histogram of regions by weighted data objects")
    plt.xlabel("Weighted data objects contained by region")
    plt.ylabel("Number of regions")
    if save_plot:
        filename = output_root + "_data_region_hist.png"
        plt.savefig(filename)
        
    else:
        plt.show()
    
    filename_ngc = output_root + "_NGC_division.png"
    filename_sgc = output_root + "_SGC_division.png"
    plot_catalog(data_cat_ngc, rand_cat_ngc, N_rows_ngc, row_edges_ngc, N_squares_ngc, square_edges_list_ngc,
                 row_ends_ngc, removed_regions=removed_regions, save_plot=save_plot, mark_regions=mark_regions,
                 output_filename=filename_ngc)
    
    plot_catalog(data_cat_sgc, rand_cat_sgc, N_rows_sgc, row_edges_sgc, N_squares_sgc, square_edges_list_sgc,
                 row_ends_sgc, removed_regions=removed_regions, save_plot=save_plot, mark_regions=mark_regions,
                 output_filename=filename_sgc, initial_index=total_regions_ngc)
    
    if write:
        print("Shifting Indices")
        for i in range(N_regions):
            if i in removed_regions:
                shift_index = i+1
                while shift_index in removed_regions:
                    shift_index += 1

                if shift_index < total_regions_ngc:
                    for obj in rand_cat_ngc:
                        if obj["JK_REGION"] == shift_index:
                            obj["JK_REGION"] = i

                    for obj in data_cat_ngc:
                        if obj["JK_REGION"] == shift_index:
                            obj["JK_REGION"] = i

                else:
                    for obj in rand_cat_sgc:
                        if obj["JK_REGION"] == shift_index:
                            obj["JK_REGION"] = i

                    for obj in data_cat_sgc:
                        if obj["JK_REGION"] == shift_index:
                            obj["JK_REGION"] = i

                removed_regions.pop(removed_regions.index(i))
                removed_regions.append(shift_index)

        if output_format=="ascii":
            print("Saving Data NGC")
            ap_ascii.write(data_cat_ngc, output="%s_NGC_jk_%i.dat" % (output_root, N_regions),
                        delimiter="\t", format="commented_header")
            print("Saving Data SGC")
            ap_ascii.write(data_cat_sgc, output="%s_SGC_jk_%i.dat" % (output_root, N_regions),
                        delimiter="\t", format="commented_header")
            print("Saving Rand NGC")
            ap_ascii.write(rand_cat_ngc, output="%s_NGC_jk_%i.rand" % (output_root, N_regions),
                        delimiter="\t", format="commented_header")
            print("Saving Rand SGC")
            ap_ascii.write(rand_cat_sgc, output="%s_SGC_jk_%i.rand" % (output_root, N_regions),
                        delimiter="\t", format="commented_header")

        else:
            print("Saving Data NGC")
            data_cat_ngc.write("%s_NGC_jk_%i_data.fits" % (output_root, N_regions), format=output_format)
            print("Saving Data SGC")
            data_cat_sgc.write("%s_SGC_jk_%i_data.fits" % (output_root, N_regions), format=output_format)
            print("Saving Rand NGC")
            rand_cat_ngc.write("%s_NGC_jk_%i_rand.fits" % (output_root, N_regions), format=output_format)
            print("Saving Rand SGC")
            rand_cat_sgc.write("%s_SGC_jk_%i_rand.fits" % (output_root, N_regions), format=output_format)

    print("N+1", N+1, "int(N+1)", int(N+1))
    plot_bins = np.linspace(-0.5, N-0.5, int(N+1))
    
    plt.figure(figsize=(16,12))
    plt.hist(rand_cat_ngc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per accepted region in NGC after index shift")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    plt.figure(figsize=(16,12))
    plt.hist(rand_cat_sgc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per accepted region in SGC after index shift")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    plt.figure(figsize=(16,12))
    plt.hist(data_cat_ngc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per aceepted region in NGC after index shift")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    plt.figure(figsize=(16,12))
    plt.hist(data_cat_sgc["JK_REGION"], bins=plot_bins)
    plt.title("Weighted randoms per accepted region in SGC after index shift")
    plt.xlabel("Region ID")
    plt.ylabel("Counts")
    plt.show()
    
    return data_cat_ngc, data_cat_sgc, rand_cat_ngc, rand_cat_sgc, removed_regions

data_table_ngc = Table.read(data_file_ngc)
rand_table_ngc = Table.read(rand_file_ngc)
data_table_sgc = Table.read(data_file_sgc)
rand_table_sgc = Table.read(rand_file_sgc)

data_table_ngc = data_table_ngc.filled()
rand_table_ngc = rand_table_ngc.filled()
data_table_sgc = data_table_sgc.filled()
rand_table_sgc = rand_table_sgc.filled()

# RA DEC WEIGHT_SYSTOT WEIGHT_CP WEIGHT_NOZ WEIGHT_FKP
weight_cols = ["WEIGHT_SYSTOT", "WEIGHT_CP", "WEIGHT_NOZ", "WEIGHT_FKP"]
for weight_col in weight_cols:
    if weight_col not in data_table_ngc.colnames:
        data_table_ngc[weight_col] = 1.

    if weight_col not in rand_table_ngc.colnames:
        rand_table_ngc[weight_col] = 1.

    if weight_col not in data_table_sgc.colnames:
        data_table_sgc[weight_col] = 1.

    if weight_col not in rand_table_sgc.colnames:
        rand_table_sgc[weight_col] = 1.

# data_table_ngc = Table.read(data_file_ngc, format="ascii")
# rand_table_ngc = Table.read(rand_file_ngc, format="ascii")
# data_table_sgc = Table.read(data_file_sgc, format="ascii")
# rand_table_sgc = Table.read(rand_file_sgc, format="ascii")

(data_cat_ngc_divided, data_cat_sgc_divided, rand_cat_ngc_divided, rand_cat_sgc_divided,
 removed_regions) = repeated_division(data_table_ngc, data_table_sgc, rand_table_ngc, rand_table_sgc, save_plot=True,
                                      mark_regions=True, write=True, area_scaling=area_scaling, verbose=False,
                                      output_root=output_root, N_regions=int(N_reg))
