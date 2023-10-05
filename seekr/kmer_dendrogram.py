#########################################################################
### Description:
# Dendrpgram for the hierarchical clustering of the rows or columns for the input dataframe
# visualize clustering results of either pearson correlation r values of seekr.pearson or p-values of find_pval

### Details:
# Customizeable dendrograms for easier visualization of the hierarchical clustering results of seekr.pearson or find_pval
# takes in a dataframe with row and column names
# performs hierarchial clustering on either rows or columns, is a better way to visualize the whole or partial clusters in kmer_heatmap
# for interested subgroups, use kmer_leiden, kmer_count_barplot, kmer_msd_barplot and kmer_textplot for further analysis

### Input:
# df: a pandas dataframe with row and column names
# df can be either the output of seekr.pearson (r values) or find_pval (p values)
# dendro_direct: the direction to perform hierarchical clustering. Can only be either 'row' (default) or 'column'
# distmetric: the distance metric to use for the dendrogram. Default is 'correlation'
# common distmetric options are: 'euclidean', 'cityblock', 'correlation', 'cosine', 'jaccard', 'hamming'
# Check the documentation for scipy.spatial.distance.pdist for the full list of valid distmetrics.
# linkmethod: the linkage method to use for the dendrogram. Default is 'complete'
# common linkmethod options are: 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward'
# Check the documentation for scipy.cluster.hierarchy.linkage for the full list of valid linkmethods.
# plot_ht: height of the dendrogram plot. must be a positive number (>0). Default is 8.  
# wd_ratio: a ratio factor to control the dendrgram width based on number of leaves. Default is 0.5.
# wd_ratio has to be a positive number (>0). if wider dendrogram is desired, increase the ratio. If narrower dendrogram is desired, decrease the ratio.
# leaf_font_size: the font size of the leaves, which corresponds to either the row or column names of the input dataframe. Default is 16
# outputname: the name and path to save the output dendrogram will automatically be combined with the trailing part _kmer_dendrogram_row/column depending on dendro_direct. 
# outputname Default 'test' will be saved under current working directory as test_kmer_dendrogram_row or test_kmer_dendrogram_column. another example is '/Users/username/Desktop/test'
# pformat: the format of the output dendrogram. Default is 'pdf'. Other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'
# pdpi: the dpi of the output dendrogram. Default is 300.

### Output:
# user defined format of file (pdf as default) of a dendrogram representing the hierarchical clustering of either the row or the column of input dataframe
# with the name outputname_kmer_dendrogram_row/column.pformat

### Example:
# kmer_dendrogram(df, dendro_direct='row', distmetric='correlation', linkmethod='complete', 
#                 plot_ht=8, wd_ratio=0.5, leaf_font_size = 16, 
#                 outputname='test', pformat='pdf', pdpi=300)

#########################################################################################################################
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt
import os


def kmer_dendrogram (df, dendro_direct='row', distmetric='correlation', linkmethod='complete', 
                     plot_ht=8, wd_ratio=0.5, leaf_font_size = 16, 
                     outputname='test', pformat='pdf', pdpi=300):
    
    if dendro_direct == 'row':

        # Perform hierarchical clustering on rows
        row_linkage = linkage(pdist(df, distmetric), linkmethod)
        # row_order = leaves_list(row_linkage)
        
        fx = round(df.shape[0]*wd_ratio)

        # Plot dendrogram
        plt.figure(figsize=(fx, plot_ht))  

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
        #plt.rcParams["axes.labelsize"] = 30
        # make sure the plot font in pdf can be editted in using illustrator
        mpl.rcParams['pdf.fonttype'] = 42

        dendrogram(row_linkage, labels=df.index, distance_sort=True, leaf_rotation=90, leaf_font_size=leaf_font_size)
        # add this to avoid truncation of x axis tick labels
        # plt.subplots_adjust(bottom=0.2)
        # save the plot
        formatlist = list(plt.gcf().canvas.get_supported_filetypes())
        if pformat in formatlist: 
            plt.savefig(f'{outputname}_kmer_dendrogram_{dendro_direct}.{pformat}', format=pformat, dpi=pdpi, bbox_inches='tight')
        else: 
            print("plotformat not supported. use default 'pdf' now. other common formats are: 'png', 'jpg', 'svg', 'eps', 'tif', 'tiff', 'ps', 'webp'")
            plt.savefig(f'{outputname}_kmer_dendrogram_{dendro_direct}.pdf', format='pdf', dpi=pdpi, bbox_inches='tight')

    elif dendro_direct == 'column':

        # Perform hierarchical clustering on columns
        col_linkage = linkage(pdist(df.T.values, distmetric), linkmethod)
        #col_order = leaves_list(col_linkage)

        if wd_ratio <= 0:
            print("wd_ratio must be a positive number (>0). Use default wd_ratio instead: 0.5")
            wd_ratio = 0.5

        fx = round(df.shape[1]*wd_ratio)

        if plot_ht <= 0:
            print("plot_ht must be a positive number (>0). Use default plot_ht instead: 8")
            plot_ht = 8

        # Plot dendrogram
        plt.figure(figsize=(fx, plot_ht))  

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
        #plt.rcParams["axes.labelsize"] = 30
        # make sure the plot font in pdf can be editted in using illustrator
        mpl.rcParams['pdf.fonttype'] = 42

        dendrogram(col_linkage, labels=df.columns, distance_sort=True, leaf_rotation=90, leaf_font_size=leaf_font_size)
        # add this to avoid truncation of x axis tick labels
        # plt.subplots_adjust(bottom=0.2)
        # save the plot
        formatlist = list(plt.gcf().canvas.get_supported_filetypes())
        if pformat in formatlist: 
            plt.savefig(f'{outputname}_kmer_dendrogram_{dendro_direct}.{pformat}', format=pformat, dpi=pdpi, bbox_inches='tight')
        else: 
            print("plotformat not supported. use default 'pdf' now. other common formats are: 'png', 'jpg', 'svg', 'eps', 'tif', 'tiff', 'ps', 'webp'")
            plt.savefig(f'{outputname}_kmer_dendrogram_{dendro_direct}.pdf', format='pdf', dpi=pdpi, bbox_inches='tight')

    else:
        print("dendro_direct must be either 'row' or 'column'. Please check and rerun.")
        return

    


