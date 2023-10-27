###################################################################################################
### Description: 
# calculte p values of the seekr.pearson correlation values for input sequence 1 vs input sequence 2 
# p value is based on the output of find_dist, which is either a list of distributions or a npy array

### Details:
# this function connects the output of find_dist and the input sequnces of interests (2 input sequences)
# given the background sequencs, find_dist calculate all possible pairwise seekr.pearson values
# and then outputs either a list of fitted distributions or the npy array of the actual data 
# find_pval firstly calculates the seekr.pearson of the two input sequences, which produces a correlation matrix 
# find_pval then calculate the p value for each r value in the pearson correlation matrix based on the output of find_dist
# the output of find_pval is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences

### Input:
# seq1file: the path to the fasta file of input sequence 1, could be a single sequence or multiple sequences
# seq2file: the path to the fasta file of input sequence 2, could be a single sequence or multiple sequences
# se1file and seq2file could be the same file, in this case, the p values are calculated for the seekr correlation between each pair of sequences in the file
# mean_path: the route to your normalization mean vector file
# std_path: the route to your normalization std vector file
# k_mer: the kmer size used in seekr, which should be consistent with the mean and std file
# fitres: the output of find_dist, which is either a list of distributions or a npy array
# log2: whether to do log2 transform, possible options are: 'Log2.post','Log2.pre' and 'Log2.none', default is 'Log2.post'
# Log2.post -- Log2 transformation post-standardization: self.counts += np.abs(np.min(self.counts)) self.counts += 1 self.counts = np.log2(self.counts)
# Log2.pre -- Log2 transformation pre-standardization: self.counts += np.abs(np.min(self.counts)) self.counts += 1
# Log2.none -- no log2 transformation
# bestfit: the index of user determined best fit distribution in the list of distributions returned by find_dist, using normal index starting from 1 (not 0), default is 1
# outputname: full path of the output file, default is None, which means the output is not saved to a csv file
# outputname can be set to 'test' under current directory or another folder: '/Users/username/Desktop/test', the trailing part _pval.csv will be automatically added
# progress_bar: whether to show the progress bar during p value caculation, showing how many rows of the final dataframe have been processed, default is True

### Output:
# a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences
# if save_file is True, the dataframe is saved to a csv file named as outputname (default is 'testpval.csv')
# if the fitres is not in the correct format, the output is None
# correct format of fitres is either a npy array 
# or a list of distributions, in the format of tuples (string, number, tuple of numbers) corresponds to (distribution name, deviance, parameters)

### Example:
# fitres = find_dist(inputseq='default', k_mer=4, log2='Log2.post', models='common10', subsetting=True, subset_size = 10000, fit_model=True, statsmethod='ks',progress_bar=True, plotfit=None, outputname='test')
# pvals=find_pval(seq1file='test1.fa', seq2file='test2.fa', mean_path='bkg_mean_4mers.npy', std_path='bkg_std_4mers.npy', k_mer=4, fitres=fitres, log2='Log2.post', bestfit=1, outputname='test', progress_bar=True)


########################################################################################################
import numpy as np
from scipy import stats
from tqdm import tqdm
import pandas as pd

from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

# check if fitres is a list consisting of tuples of (string, number, tuple of numbers)
# th number here could be float64 or float32
def is_float_type(x):
    return isinstance(x, float) or np.isscalar(x)

def check_tuple_format(tup):
    if len(tup) != 3:
        return False
    return isinstance(tup[0], str) and \
           is_float_type(tup[1]) and \
           isinstance(tup[2], tuple) and \
           all(is_float_type(x) for x in tup[2])

def check_main_list(main_list):
    return all(check_tuple_format(tup) for tup in main_list)


def find_pval(seq1file, seq2file, mean_path, std_path, k_mer, fitres, log2='Log2.post', bestfit=1, outputname=None, progress_bar=True):

    t1 = seekrBasicCounter(seq1file, mean=mean_path, std=std_path, k=k_mer, log2=log2, silent=True)
    t2 = seekrBasicCounter(seq2file, mean=mean_path, std=std_path, k=k_mer, log2=log2, silent=True)

    t1.make_count_file()
    t2.make_count_file()

    sim = seekrPearson(t1.counts,t2.counts)

    header1=seekrReader(seq1file).get_headers()
    header2=seekrReader(seq2file).get_headers()
    # remove the '>' in the header
    header1 = [i[1:] for i in header1]
    header2 = [i[1:] for i in header2]

    # check if header1 only contains unique elements
    if len(header1) != len(set(header1)):
        print('The headers of seq1file is not unique.')
        print('Be carefule during further analysis as there are potential indexing problems.')

    # check if header2 only contains unique elements
    if len(header2) != len(set(header2)):
        print('The headers of seq2file is not unique.')
        print('Be carefule during further analysis as there are potential indexing problems.')
        
    # check if fitres is a distribution or a npy array
    if isinstance(fitres, list):
        # check if fitres is a list consisting of (string, number, tuple of numbers)
        if not check_main_list(fitres):
            print('The format of fitres is wrong.')
            print('fitres should be a list consisting of tuples (string, number, tuple of numbers) corresponds to (distribution name, deviance, parameters)')
            print('fitres should be the output of find_dist.')
            print('No p value is calculated. The output is None.')
            return None
        else: 
            # get the best fit distribution by user input
            # for python the counting starts from 0
            bestdist = fitres[bestfit-1]
            distname = bestdist[0]
            params = bestdist[2]

            #Load distribution
            dist = getattr(stats, distname)
            distribution = dist(*params)

            # Initialize a matrix to hold the p-values
            p_values = np.zeros_like(sim)

            # Loop through each element in sim
            iterable = tqdm(range(sim.shape[0]), desc="Rows") if progress_bar else range(sim.shape[0])
            for i in iterable:
                for j in range(sim.shape[1]):
                    p_values[i, j] = 1 - distribution.cdf(sim[i, j])
            # convert to dataframe and adding row/column names
            pval_df = pd.DataFrame(p_values, index=header1, columns=header2)
            # save the dataframe to csv file
            if outputname:
                pval_df.to_csv(f'{outputname}_pval.csv')

            return pval_df



    elif isinstance(fitres, np.ndarray):
        # check if fitres is a 1D array
        if len(fitres.shape) == 1:
            # Calculate p-values for entire sim matrix
            total_length_fitres = len(fitres)

            # alternative method to the loop below, faster but requires more memory
            # Use broadcasting to compare each element in sim with all elements in fitres
            # The result is a 3D array where the third dimension has size equal to total_length_fitres
            # comparison_result = fitres > sim[:,:,np.newaxis]
            # Sum along the third dimension and divide by total_length_fitres to get p-values
            # p_values = np.sum(comparison_result, axis=-1) / total_length_fitres

            # Initialize a matrix to hold the p-values
            p_values = np.zeros_like(sim)

            # Loop through each element in sim
            iterable = tqdm(range(sim.shape[0]), desc="Rows") if progress_bar else range(sim.shape[0])
            for i in iterable:
                for j in range(sim.shape[1]):
                    p_values[i, j] = np.sum(fitres > sim[i, j]) / total_length_fitres

            # convert to dataframe and adding row/column names
            pval_df = pd.DataFrame(p_values, index=header1, columns=header2)
            # save the dataframe to csv file
            if outputname:
                pval_df.to_csv(f'{outputname}_pval.csv')
            return pval_df 
        
        else:
            print('The dimension of fitres as a numpy array is wrong. fitres should be a 1D numpy array.')
            print('fitres should be the output of find_dist.')
            print('No p value is calculated. The output is None.')
            return None


    else:
        print('fitres should be the output of find_dist. It should be either a list of distributions or a numpy array.')
        print('No p value is calculated. The output is None.')
        return None



