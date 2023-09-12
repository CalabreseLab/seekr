###################################################################################################
### Description:
# barplot of the mean of the transformed or raw z-score for kmer words across input sequences
# error bar is the standard deviation of the transformed or raw z-score for kmer words across input sequences

### Details:
# the z-score that is the output of BasicCounter, where user can define whether and how to do log2 transform
# get the mean and standard deviation (sd) of the z-score for each kmer word across input sequences
# can choose to order kmer words by mean or sd, can either be descending or ascending
# user define to plot the top x kmer words (could be messy if there are too many words)
# x axis label is the kmer words

### Input:
# inputfile: fasta file with unique headers for each sequences
# mean: the route to your normalization mean vector file
# std: the route to your normalization std vector file
# k: the kmer size used in seekr, which should be consistent with the mean and std file
# log2: whether to do log2 transform, possible options are: 'Log2.post','Log2.pre' and 'Log2.none', default is 'Log2.post'
# Log2.post -- Log2 transformation post-standardization: self.counts += np.abs(np.min(self.counts)) self.counts += 1 self.counts = np.log2(self.counts)
# Log2.pre -- Log2 transformation pre-standardization: self.counts += np.abs(np.min(self.counts)) self.counts += 1
# Log2.none -- no log2 transformation
# sortstat: whether to sort the kmer words based on mean or sd, possible options: 'mean' (default) and 'sd'
# sortmethod: how to sort the kmer words based on selected sortstat, possible options: 'ascending' and 'descending' (default)
# topkmernumber: the number of sorted kmer words to plot, default is 10, if input number is more than total number of kmer words, plot all words
# if you want to plot all kmer words, set topkmernumber to a extreme number such as 1000000
# xlabelsize: x axis label font size, default is 20
# ylabelsize: y axis label font size, default is 20
# xticksize: x tick label font size, default is 20
# yticksize: y tick label font size, default is 20
# lengendsize: legend font size, default is 12
# outputpath: the path to save the output file, default is current working directory, example: '/Users/username/Desktop/'
# outputname: the name of the output file, default is 'test'
# pformat: the format of the output file, default is 'pdf'. Other options are: 'eps', 'jpeg', 'jpg', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff', 'webp'
# pdpi: the dpi of the output file, default is 300

### Output:
# pdf file of the barplot named as outputname_kmer_msd_barplot.pdf

### Example:
# kmer_msd_barplot(inputfile='test.fa', mean='gencodevM25_unique_mean_4mers.npy',
#                  std = 'gencodevM25_unique_std_4mers.npy', log2 = 'Log2.post',
#                  k=4, sortstat='sd', sortmethod='ascending', topkmernumber=10,
#                  xlabelsize=20, ylabelsize=20, xticksize=20, yticksize=20, lengendsize=12,
#                  outputpath='', outputname='test', pformat='pdf', pdpi=300)


#######################################################################################################################
import matplotlib.font_manager as font_manager
import matplotlib as mpl
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.fasta_reader import Reader as seekrReader


def kmer_msd_barplot(inputfile, mean, std, k, log2='Log2.post', sortstat='mean',sortmethod='descending', 
                     topkmernumber=10, xlabelsize=20, ylabelsize=20, xticksize=20, yticksize=20, lengendsize=12, 
                     outputpath='', outputname='test', pformat='pdf', pdpi=300):

    t1 = seekrBasicCounter(inputfile, mean=mean, std=std, log2=log2, k=k,silent=True)

    t1.make_count_file()
    header1=seekrReader(inputfile).get_headers()
    # remove the '>' in the header
    header1 = [i[1:] for i in header1]

    # get counts
    tt1=t1.counts

    # get kmers and convert to dataframe
    kmer1=t1.kmers
    df = pd.DataFrame(tt1, index=header1, columns=kmer1)

    # Calculate the mean and sd for each column
    column_means = df.mean()
    column_sd = df.std()

    # Sort by the means or standard deviations in descending or ascending order
    if sortstat == 'mean': 
        if sortmethod == 'ascending':
            sorted_indices = column_means.sort_values(ascending=True).index
        elif sortmethod == 'descending':
            sorted_indices = column_means.sort_values(ascending=False).index
        else:
            print("Please choose a sorting method: 'ascending' or 'descending', use default 'descending' now")
            sorted_indices = column_means.sort_values(ascending=False).index
    elif sortstat == 'sd':
        if sortmethod == 'ascending':
            sorted_indices = column_sd.sort_values(ascending=True).index
        elif sortmethod == 'descending':
            sorted_indices = column_sd.sort_values(ascending=False).index
        else:
            print("Please choose a sorting method: 'ascending' or 'descending', use default 'descending' now")
            sorted_indices = column_sd.sort_values(ascending=False).index
    else:
        print("Please choose a sorting stat: 'mean' or 'sd', use default 'mean' and default sortmethod'descending' now")
        sorted_indices = column_means.sort_values(ascending=False).index


    # Reorder the dataframe based on the summed absolute difference from the mean
    df_reordered = df[sorted_indices]


    # Melt the dataframe to a long format
    df_melted = df_reordered.reset_index().melt(id_vars='index', value_vars=df_reordered.columns)
    df_melted.columns = ['Sample', 'Kword', 'Value']


    # calculate the number of input sequences
    ttlnum = topkmernumber*len(header1)

    if ttlnum > len(df_melted):
        print(f"Only {int(len(df_melted)/len(header1))} kmer words, less than {topkmernumber} words you want to plot, plot all words")
        df_melted_plot=df_melted
    else:
        df_melted_plot=df_melted[:ttlnum]

    # Plot
    plt.figure(figsize=((topkmernumber*2), 8))

    # use Airal for fonts
    # Specify the path to your .ttf font file
    font_path = 'arial.ttf'
    # Register the font with matplotlib
    font_manager.fontManager.addfont(font_path)
    # Set the font properties
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = prop.get_name()
    # make sure the plot font in pdf can be editted in using illustrator
    mpl.rcParams['pdf.fonttype'] = 42

    # plot the barplot
    sns.barplot(x='Kword', y='Value', data=df_melted_plot, palette='tab10', errorbar='sd',capsize=0.2)

    # Setting x and y axis labels
    plt.xlabel('Kmer Words', fontsize=xlabelsize)
    plt.ylabel('z-score (transformed or raw)', fontsize=ylabelsize)

    # Adjusting tick font sizes
    plt.xticks(rotation=90, fontsize=xticksize)
    plt.yticks(fontsize=yticksize)

    # Place the legend to the right of the plot
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=lengendsize)


    #plt.show()
    formatlist = list(plt.gcf().canvas.get_supported_filetypes())
    if pformat in formatlist:
        plt.savefig(f'{outputpath}{outputname}_kmer_msd_barplot.{pformat}', format=pformat, dpi=pdpi, bbox_inches='tight')
    else:
        print("plotformat not supported. use default 'pdf' now. other common formats are: 'png', 'jpg', 'svg', 'eps', 'tif', 'tiff', 'ps', 'webp'")
        plt.savefig(f'{outputpath}{outputname}_kmer_msd_barplot.pdf', format='pdf', dpi=pdpi, bbox_inches='tight')


