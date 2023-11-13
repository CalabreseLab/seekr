#########################################################################
### Description:
# Heatmap with both rows and columns dendrograms for input dataframe
# visualize results of either pearson correlation r values of seekr.pearson or p-values of find_pval

### Details:
# Customizeable heatmap for easier visualization of the results of seekr.pearson or find_pval
# takes in a dataframe with row and column names
# performs hierarchial clustering on both rows and columns, and then reorder the cols and rows based on the dendrograms
# can also use kmer_dendrgram to plot only the dendrograms (partial or full) to get a better idea of the clustering
# for interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot and kmer_textplot for further analysis

### Input:
# df: a pandas dataframe with row and column names
# df can be either the output of seekr.pearson (r values) or find_pval (p values)
# datamin: the minimum possible value of the data. For r values, this is -1. For p values, this is 0.
# datamax: the maximum possible value of the data. For r values, this is 1. For p values, this is 1.
# thresh_value: the threshold value. Default is 0.05. If using a 3 color palette, this corresponds to the middle color. 
# thresh_value can be used to separate the data into two colors based on a threshold
# for example, if using r values, thresh_value=0 will separate the data into positive and negative correlations
# if using p values, thresh_value=0.05 will separate the data into significant and non-significant p values
# color_range: the color palette to use for the heatmap, has to be a list of 2 or 3 hex color strings. 
# color_range Default is ['#1b7837', '#ffffff', '#c51b7d'] which is [green, white, magenta]
# if using three colors, the middle color will corresponds to the thresh_value
# if using two colors, the thresh_value will be ignored
# the first and last colors correspond to the datamin and datamax respectively
# distmetric: the distance metric to use for the dendrogram. Default is 'correlation'
# common distmetric options are: 'euclidean', 'cityblock', 'correlation', 'cosine', 'jaccard', 'hamming'
# Check the documentation for scipy.spatial.distance.pdist for the full list of valid distmetrics.
# linkmethod: the linkage method to use for the dendrogram. Default is 'complete'
# common linkmethod options are: 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'
# Check the documentation for scipy.cluster.hierarchy.linkage for the full list of valid linkmethods.
# hmapw_ratio: a ratio factor to control the heatmap width based on number of columns. Default is 0.3. 
# hmapw_ratio has to be a positive number (>0). if wider heatmap is desired, increase the ratio. If narrower heatmap is desired, decrease the ratio.
# hmaph_ratio: a ratio factor to control the heatmap height based on number of rows. Default is 0.3.
# hmaph_ratio has to be a positive number (>0). if taller heatmap is desired, increase the ratio. If shorter heatmap is desired, decrease the ratio.
# x_tick_size: the font size for the x tick labels or the column labels. Default is 16
# y_tick_size: the font size for the y tick labels or the row labels. Default is 16
# cbar_font_size: the font size for the colorbar ticks. Default is 16
# outputname: the path and name to save the output heatmap
# outputname Default, 'test_kmer_heatmap', which is saved under current working directory, another example is '/Users/username/Desktop/test_kmer_heatmap'
# hformat: the format of the output heatmap. Default is 'pdf'. Other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'
# hdpi: the dpi of the output heatmap. Default is 300.

### Output:
# user defined format of file (pdf as default) of a heatmap with both rows and columns dendrograms with the name outputname.hformat

### Example:
# kmer_heatmap(df, datamin=0, datamax=1, thresh_value=0.05,
#             color_range=['#1b7837', '#ffffff', '#c51b7d'],
#             distmetric='correlation', linkmethod='complete', 
#             hmapw_ratio=0.3, hmaph_ratio=0.3, x_tick_size=6, y_tick_size=6, cbar_font_size=16, 
#             outputname='test_kmer_heatmap',hformat='pdf', hdpi=300)

#########################################################################################################
import numpy as np
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt

import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
from scipy.spatial.distance import pdist
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import re
import os

def is_hex_color(s):
    return re.fullmatch(r'#[0-9a-fA-F]{6}', s) is not None

def check_hex_colors(lst):
    return all(is_hex_color(color) for color in lst)

def kmer_heatmap(df, datamin, datamax, thresh_value=0.05,
                 color_range=['#1b7837', '#ffffff', '#c51b7d'], 
                 distmetric='correlation', linkmethod='complete',  
                 hmapw_ratio=0.3, hmaph_ratio=0.3, 
                 x_tick_size=16, y_tick_size=16, cbar_font_size=16, 
                 outputname='test_kmer_heatmap',hformat='pdf', hdpi=300):

    # get matrix out of the dataframe
    data = df.values
    xheaders = df.columns
    yheaders = df.index

    # Create custom colormap
    turnval=(thresh_value-datamin)/(datamax-datamin)
    if not check_hex_colors(color_range):
        print("color_range must be a list of valid hex colors (for example '#ffffff').")
        print("Use default color_range instead: ['#1b7837', '#ffffff', '#c51b7d']")
        color_range=['#1b7837', '#ffffff', '#c51b7d']
        cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', [(0, color_range[0]), (turnval, color_range[1]), (1, color_range[2])])
    if (len(color_range) < 2 or len(color_range) >3):
        print("color_range must have 2 or 3 colors. Check color_range list length.")
        print("Use default color_range instead: ['#1b7837', '#ffffff', '#c51b7d']")
        color_range=['#1b7837', '#ffffff', '#c51b7d']
        cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', [(0, color_range[0]), (turnval, color_range[1]), (1, color_range[2])])
    elif (len(color_range) == 2):
        cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', [(0, color_range[0]), (1, color_range[1])])
    else:
        cmap = LinearSegmentedColormap.from_list(
        'custom_cmap', [(0, color_range[0]), (turnval, color_range[1]), (1, color_range[2])])


    # Perform hierarchical clustering on rows
    try:
        row_linkage = linkage(pdist(data, metric=distmetric), method=linkmethod)
    except ValueError as e:
        if "Unknown Distance Metric" in str(e):
            print(f"The specified distance metric '{distmetric}' is not supported.")
            print("Check the documentation for scipy.spatial.distance.pdist for a list of valid metrics.")
            raise(e)
        elif "Invalid method" in str(e):
            print(f"The specified linkage method '{linkmethod}' is not supported.")
            print("Check the documentation for scipy.cluster.hierarchy.linkage for a list of valid methods.")
            raise(e)
        else:
            raise e

    row_order = leaves_list(row_linkage)

    # Perform hierarchical clustering on columns
    try:
        col_linkage = linkage(pdist(data.T, metric=distmetric), method=linkmethod)
    except ValueError as e:
        if "Unknown Distance Metric" in str(e):
            print(f"The specified distance metric '{distmetric}' is not supported.")
            print("Check the documentation for scipy.spatial.distance.pdist for a list of valid metrics.")
            raise(e)
        elif "Invalid method" in str(e):
            print(f"The specified linkage method '{linkmethod}' is not supported.")
            print("Check the documentation for scipy.cluster.hierarchy.linkage for a list of valid methods.")
            raise(e)
        else:
            raise e
        
    col_order = leaves_list(col_linkage)

    # Reorder data
    data_clustered = data[row_order, :][:, col_order]

    if hmapw_ratio<=0:
        print("hmapw_ratio must be a positive number (>0). Use default hmapw_ratio instead: 0.3")
        hmapw_ratio=0.3

    if hmaph_ratio<=0:
        print("hmaph_ratio must be a positive number (>0). Use default hmaph_ratio instead: 0.3")
        hmaph_ratio=0.3

    fx = round(len(xheaders)*hmapw_ratio)
    fy = round(len(yheaders)*hmaph_ratio)

    # Create a figure and gridspec layout
    # plot the main heatmap and dendrogram in the left plot
    # plot the colorbar in the right plot
    plt.figure(figsize=(fx+3, fy+1))
    gs = GridSpec(1, 2, width_ratios=[fx+1, 2])

    # use Airal for fonts
    # Specify the path to your .ttf font file
    # Get the directory where the current file is located.
    current_dir = os.path.dirname(os.path.realpath(__file__))

    # Construct the path to the .ttf file.
    font_path = os.path.join(current_dir, 'data', 'arial.ttf')
    # Register the font with matplotlib
    font_manager.fontManager.addfont(font_path)
    # Set the font properties
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()

    # make sure the plot font in pdf can be editted in using illustrator
    mpl.rcParams['pdf.fonttype'] = 42


    # Create the main subplot for heatmap and dendrograms (left)
    ax_main = plt.subplot(gs[0])
    # Add the row dendrogram (left of heatmap)
    ax_row_dendrogram = ax_main.inset_axes([0.05, 0.1, 0.2, 0.65])
    # The x-coordinate of the left side of the axes (5% from the left edge of the figure).
    # The y-coordinate of the bottom side of the axes (10% from the bottom edge of the figure).
    # The width of the axes, which is 20% of the figure's width.
    # The height of the axes, which is 60% of the figure's height.
    dendrogram(row_linkage, orientation='left', ax=ax_row_dendrogram, color_threshold=0)
    ax_row_dendrogram.set_axis_off()

    # Add the column dendrogram (above heatmap)
    ax_col_dendrogram = ax_main.inset_axes([0.26, 0.76, 0.65, 0.2])
    dendrogram(col_linkage, ax=ax_col_dendrogram, color_threshold=0)
    ax_col_dendrogram.set_axis_off()


    # Plot the heatmap within the main subplot using fig.add_axes
    ax_heatmap = ax_main.inset_axes([0.26, 0.1, 0.65, 0.65])
    sns.heatmap(data_clustered, cmap=cmap, vmin=datamin, vmax=datamax,
                yticklabels=np.array(yheaders)[row_order],
                xticklabels=np.array(xheaders)[col_order], cbar=False,
                ax=ax_heatmap)
    ax_heatmap.yaxis.set_ticks_position('right')
    # Set rotation and font size for y-tick labels
    ax_heatmap.tick_params(axis='y', rotation=0, labelsize=y_tick_size)
    # Set rotation and font size for x-tick labels
    ax_heatmap.tick_params(axis='x', rotation=90, labelsize=x_tick_size)

    # hide ticks and boundary lines for the main figure
    ax_main.set_xticks([])
    ax_main.set_yticks([])
    for spine in ax_main.spines.values():
        spine.set_visible(False)

    # Create the second subplot for colorbar (right)
    ax_cbar_main = plt.subplot(gs[1])

    #Define the region within this subplot for the colorbar (left, bottom, width, height)
    # ax_cbar = ax_cbar_main.inset_axes([ax_heatmap.get_position().x0,  # Adjust left position
    #                                    ax_heatmap.get_position().y0 + 0.1,  # Adjust bottom position
    #                                    1,  # Width (e.g., 20% of subplot width)
    #                                    0.65]) # Height (e.g., 40% of subplot height)

    ax_cbar = ax_cbar_main.inset_axes([0.3,0.1,1,0.65]) 
    # Add colorbar to the right panel
    # set the height of the colorbar to be the same as the heatmap
    # set the width of the colorbar to be 100% of the total figure width

    cbar = plt.colorbar(ax_heatmap.collections[0], ax=ax_cbar, fraction=1, pad=0, anchor=(0,0), aspect=30)
    #ax4 = fig.add_axes([0.95, 0.1, 0.01, 0.65])
    #cbar.set_label('Correlation')
    # Set tick font size for colorbar
    cbar.ax.tick_params(labelsize=cbar_font_size)
    
    # Get current ticks from colorbar
    current_ticks = cbar.get_ticks()
    # Add user input to ticks if it doesn't already exist
    if thresh_value not in current_ticks:
        new_ticks = np.append(current_ticks, thresh_value)
        new_ticks = np.sort(new_ticks)  # Sort the ticks
        cbar.set_ticks(new_ticks)

    # Set the colorbar layer to be behind other elements 
    ax_cbar_main.set_zorder(-1) 


    # Hide ticks and boundary lines for the colorbar axis
    ax_cbar_main.set_xticks([])
    ax_cbar_main.set_yticks([])
    for spine in ax_cbar_main.spines.values():
        spine.set_visible(False)

    ax_cbar.set_xticks([])
    ax_cbar.set_yticks([])
    for spine in ax_cbar.spines.values():
        spine.set_visible(False)

    #plt.subplots_adjust(bottom=0.2)
    #plt.savefig(f'{outputname}.{hformat}', format=hformat, dpi=hdpi, bbox_inches='tight')
    formatlist = list(plt.gcf().canvas.get_supported_filetypes())
    if hformat in formatlist: 
        plt.savefig(f'{outputname}.{hformat}', format=hformat, dpi=hdpi, bbox_inches='tight')
    else: 
        print("plotformat not supported. use default 'pdf' now. other common formats are: 'png', 'jpg', 'svg', 'eps', 'tif', 'tiff', 'ps', 'webp'")
        plt.savefig(f'{outputname}.pdf', format='pdf', dpi=hdpi, bbox_inches='tight')



