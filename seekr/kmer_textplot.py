###################################################################################################
### Description:
# highlight the input kmer words in the input sequences

### Details:
# plot 2 input sequences and several interested words (max 10)
# first input seq is plot in black and second input seq is plot in grey
# limit the number of interested words to 10 as it is hard to distinguish more than 10 colors
# align the sequence at the beginning and label sequence positions on the bottom
# highlight in colors and bold the words of interest
# for multiple words, if the words overlap, the overlapped characters will be highlighted in the color of the first word in the list
# therefore arrange the words in the list based on the priority of the words, with the most important word at the beginning of the list

### Input:
# seq1file, seq2file: input sequence 1 and 2, both in fasta format, if seq1 and seq2 includes more than one sequence, only the first sequence will be plotted
# seq1file will be plot in black while seq2file will be plot in grey
# words: word of interest, max 10 words in the format of a list, example: ['ATTA','AAAA','ACTC','CCTT','GGCC'], if more than 10 words, only the first 10 will be plotted
# color_vec: list of hex color for the words of interest (for example: ['#d62728','#e377c2','#ff7f0e']), default will use the 'tab10' pallette with rearrangements into a quasi-rainbow order
# the order of default color: red, pink, orange, olive, green, cyan, blue, purple, brown, grey
# if you want to use your own color vector, please make sure the length of the color vector is the same as the length of the words list
# wraplen: wrapping length, how many words to wrap in each line, default is 60
# char_spacing: space between characters in the plot, default is 0.1
# line_spacing: line space between seq1, seq2 and number, default is 0.2
# seqfontsize: sequence character font size, default is 72
# numfontsize: sequence position number font size, default is 40
# colorblockh: the height of the highlight color block, default is 0.15, change this when change seqfontsize
# outpoutname: the name and path of the output file, default is under current directory 'test_kmer_textplot', other example: '/Users/username/Desktop/test_kmer_textplot'
# plotformat: output format, default is 'pdf', other common options are 'png', 'jpg', 'svg', 'eps', 'tif', 'tiff', 'ps', 'webp'
# you can check all the options here: https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html
# or by running this code: matplotlib.pyplot.gcf().canvas.get_supported_filetypes()
# plotdpi: resolution of the output file, default is 300

### Output:
# user defined format of file (pdf as default) of the highlighted text plot named as outputname.plotformat

### Example:
# kmer_textplot(seq1file='test1.fa', seq2file='test2.fa', words=['TAAA','GGTG','TCCA','GGTG','ACCT','ATAC','AGAC','GTCC','TTTT','ACCT'], 
#               color_vec='default', wraplen=60, char_spacing=0.1, line_spacing=0.2,
#               seqfontsize=72, numfontsize=40, colorblockh=0.15, outputname='test_kmer_textplot', plotformat='pdf', plotdpi=300)


#######################################################################################################################
import matplotlib.font_manager as font_manager
import matplotlib as mpl
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
import os

from seekr.fasta_reader import Reader as seekrReader

# loop thru seq and find unique coordinates for the word that matches the input word
def find_word_coordinates(seq, inputword):
    word_length = len(inputword)
    coords = []
    
    # Loop through each character of seq
    for i in range(len(seq) - word_length + 1):
        seqword = seq[i:i + word_length]
        
        # If there's a match, record the coordinates
        if seqword == inputword:
            coords.extend(list(range(i, i + word_length)))

    # Convert to numpy array and get unique coordinates
    unique_coords = np.unique(np.array(coords))

    return unique_coords

# generate color for input coord 
def ass_color(coord, matched_seq, color_vec):

    for n in range(len(matched_seq)):
        if coord in matched_seq[n]:
            return color_vec[n]



def kmer_textplot(seq1file, seq2file, words, color_vec='default', wraplen=60, char_spacing=0.1, line_spacing=0.2, 
                  seqfontsize=72, numfontsize=40, colorblockh=0.15, outputname='test_kmer_textplot', plotformat='pdf', plotdpi=300):

    # read in seq1 and seq2
    seq1 = seekrReader(seq1file).get_seqs()[0]
    seq2 = seekrReader(seq2file).get_seqs()[0]
        

    header1=seekrReader(seq1file).get_headers()[0]
    header2=seekrReader(seq2file).get_headers()[0]
    # remove the '>' in the header
    header1=header1[1:]
    header2=header2[1:]

    # find matched coordinates for each word
    if len(words)<=10:
        matched_seq1 = [find_word_coordinates(seq1, word) for word in words]
        matched_seq2 = [find_word_coordinates(seq2, word) for word in words]
    else: 
        print("length of words list exceeds 10, plotting the first 10 only")
        words = words[0:10]
        matched_seq1 = [find_word_coordinates(seq1, word) for word in words]
        matched_seq2 = [find_word_coordinates(seq2, word) for word in words]

    # convert the list of lists to a 1D np array and only get the unique values
    matched_seq1_fl = np.unique(np.array([item for sublist in matched_seq1 for item in sublist]))
    matched_seq2_fl = np.unique(np.array([item for sublist in matched_seq2 for item in sublist]))

    if color_vec == 'default':
        color_vec = ['#d62728','#e377c2','#ff7f0e','#bcbd22','#2ca02c','#17becf','#1f77b4','#9467bd','#8c564b','#7f7f7f']
        print('default color order: red, pink, orange, olive, green, cyan, blue, purple, brown, grey')
    else: 
        if len(color_vec) != len(words):
            print('the length of color vector is not the same as the length of the words list, use default color now')
            print('default color order: red, pink, orange, olive, green, cyan, blue, purple, brown, grey')
            color_vec = ['#d62728','#e377c2','#ff7f0e','#bcbd22','#2ca02c','#17becf','#1f77b4','#9467bd','#8c564b','#7f7f7f']

    wrapped_seq1 = [seq1[i:i+wraplen] for i in range(0, len(seq1), wraplen)]
    wrapped_seq2 = [seq2[i:i+wraplen] for i in range(0, len(seq2), wraplen)]

    # Determine the total number of lines to plot
    total_lines = max(len(wrapped_seq1), len(wrapped_seq2))

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(wraplen/4, total_lines*1))  # 3 rows for each line


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

    # Hide the x and y axis
    ax.axis('off')

    # Plot sequences and numbers
    # for matched characters use red color and bold font
    for i in range(total_lines):
        for j in range(wraplen):
            
            rect_height = colorblockh 
            rect_width = char_spacing 
            
            if i < len(wrapped_seq1) and j < len(wrapped_seq1[i]):
                color = ass_color(i*wraplen+j, matched_seq1, color_vec) if i*wraplen+j in matched_seq1_fl else 'none'
                rect1 = Rectangle((j*char_spacing - rect_width/2, total_lines*1 - 1*i - 0.1 - rect_height/2), rect_width, rect_height, color=color, linewidth=0)
                rect1.set_clip_on(False) # enforce plotting patch even if outside of ax
                ax.add_patch(rect1)
                weight = 'bold' if i*wraplen+j in matched_seq1_fl else 'normal'
                ax.text(j*char_spacing, total_lines*1 - 1*i - 0.1, wrapped_seq1[i][j], fontsize=seqfontsize, color='#000000',ha='center', va='center', weight=weight)
                
                
            if i < len(wrapped_seq2) and j < len(wrapped_seq2[i]):
                color = ass_color(i*wraplen+j, matched_seq2, color_vec) if i*wraplen+j in matched_seq2_fl else 'none'
                rect2 = Rectangle((j*char_spacing - rect_width/2, total_lines*1 - 1*i - 0.1 - line_spacing - rect_height/2), rect_width, rect_height, color=color, linewidth=0)
                rect2.set_clip_on(False) # enforce plotting patch even if outside of ax
                ax.add_patch(rect2)
                weight = 'bold' if i*wraplen+j in matched_seq2_fl else 'normal'                
                ax.text(j*char_spacing, total_lines*1 - 1*i - 0.1 - line_spacing, wrapped_seq2[i][j], fontsize=seqfontsize, ha='center', va='center', color='#838383', weight=weight)
                
            if i*wraplen + j < max(len(seq1), len(seq2)):
                ax.text(j*char_spacing, total_lines*1 - 1*i - 0.1-line_spacing*2, str(i*wraplen+j+1), fontsize=numfontsize, ha='center', va='center')
    
    #plt.show()
    formatlist = list(plt.gcf().canvas.get_supported_filetypes())
    if plotformat in formatlist: 
        plt.savefig(f'{outputname}.{plotformat}', format=plotformat, dpi=plotdpi, bbox_inches='tight')
    else: 
        print("plotformat not supported. use default 'pdf' now. other common formats are: 'png', 'jpg', 'svg', 'eps', 'tif', 'tiff', 'ps', 'webp'")
        plt.savefig(f'{outputname}.pdf', format='pdf', dpi=plotdpi, bbox_inches='tight')

