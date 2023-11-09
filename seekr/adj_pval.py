###################################################################################################
### Description: 
# adjust p-values for multiple comparisons

### Details:
# this function performs multiple comparison correction for p values calculated by find_pval
# multiple comparison correction is performed by statsmodels.stats.multitest.multipletests function
# the input is a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences
# the output is a dataframe of adjusted p values, with the same dimension as the input dataframe
# if the input dataframe is a symmetric matrix: the row and column (or input1 and input2 in find_pval) are the same
# the number of multiple comparisons is automatically corrected to be half the matrix size, excluding the diagonal
# that is only the upper triangle of the matrix (excluding diagonal) is used for multiple comparison correction
# in the output dataframe, the upper triangle of the matrix (excluding diagonal) is filled with the adjusted p values
# the lower triangle and the diagonal is filled with NaN
# if the input dataframe is not a symmetric matrix the total matrix is used for multiple comparison correction

### Input:
# pvals: a dataframe of p values, with the row names (input 1) and column names (input 2) as the headers of the input sequences from find_pval
# method: the method used for multiple comparison correction; can be 'bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky'
# please refer to statsmodels.stats.multitest.multipletests for details
# alpha: the desired alpha level (family-wise error rate, FWER); default is 0.05
# The FWER is the probability of making one or more false discoveries, or Type I errors, when performing multiple hypotheses tests
# alpha represents the desired significance level for the entire set of tests, not just for an individual hypothesis test. 
# For example, if you set alpha to 0.05, you are aiming to have no more than a 5% chance of making at least one Type I error across all the tests you are conducting
# outputname: full path of the output file, default is None, which means the output is not saved to a csv file
# outputname can be set to 'test' under current directory or another folder: '/Users/username/Desktop/test'
# the trailing part _{method}_alpha{alpha}_adjusted_pval.csv will be automatically added


### Output:
# a dataframe of adjusted p values, with the row names (input 1) and column names (input 2) same as the input dataframe
# if outputname is set to a string, the dataframe is saved to a csv file named as outputname_{method}_alpha{alpha}_adjusted_pval.csv
# if the input pvals is a symmetric matrix, meaning the row and column (or input1 and input2 in find_pval) are the same
# the output dataframe will only contain adjusted pvals for the upper triangle of the matrix (excluding diagonal)
# the rest of the matrix is set to NaN, as they are not used for multiple comparison correction

### Example:
# fitres = find_dist(inputseq='default', k_mer=4, log2='Log2.post', models='common10', subsetting=True, subset_size = 10000, fit_model=True, statsmethod='ks',progress_bar=True, plotfit=None, outputname='test')
# pvals=find_pval(seq1file='test1.fa', seq2file='test2.fa', mean_path='bkg_mean_4mers.npy', std_path='bkg_std_4mers.npy', k_mer=4, fitres=fitres, log2='Log2.post', bestfit=1, outputname='test', progress_bar=True)
# adjpvals=adj_pval(pvals, method='bonferroni', alpha=0.05, outputname='test')

########################################################################################################

import pandas as pd
import numpy as np
import statsmodels.stats.multitest as ssm

# pvals=pd.read_csv('test_pval_sym.csv',header=0, index_col=0)

# Function to check if the DataFrame is symmetric
def is_symmetric(df):
    # Exclude the diagonal elements by filling them with NaN
    df_no_diagonal = df.where(~np.eye(df.shape[0], dtype=bool))
    # Check if DataFrame is symmetric
    return df_no_diagonal.equals(df_no_diagonal.T)

def adj_pval(pvals, method, alpha=0.05, outputname=None):

    # check if pvals is a dataframe
    if not isinstance(pvals, pd.DataFrame):
        print('The input pvals is not a dataframe. Please check the input.')
        return None
    else:
        # check if column and row numbers are the same for pvals
        # as if they are not the same it is definitely not a symmetric matrix
        if pvals.shape[0] == pvals.shape[1]:
            # check if pvals is a symmetric matrix
            if is_symmetric(pvals):
                print('The input pvals is a symmetric matrix. Only the upper triangle of the matrix (excluding diagonal) is used for multiple comparison correction.')
                # pvals is a symmetric matrix
                # only get the upper triangle of the matrix, excluding the diagonal
                # and flatten it to numpy array
                upper_triangle_indices = np.triu_indices_from(pvals, k=1)
                # Extract the p-values from the upper triangle
                upper_tri_values = pvals.values[upper_triangle_indices]
                # perform multiple comparison correction
                adj_p_vals = ssm.multipletests(upper_tri_values, alpha=alpha, method=method)[1]
    
                # Create a new DataFrame to hold the adjusted p-values with the same shape as the original
                adj_df = pd.DataFrame(np.nan, index=pvals.index, columns=pvals.columns)
                
                # Fill the upper triangle of the adjusted DataFrame with the adjusted p-values
                adj_df.values[upper_triangle_indices] = adj_p_vals

                if outputname:
                    adj_df.to_csv(f'{outputname}_{method}_alpha{alpha}_adjusted_pval.csv')
                return adj_df
            
            else:
                print('The input pvals is not a symmetric matrix. The total matrix is used for multiple comparison correction.')
                # convert and faltten the dataframe
                pvals_np = pvals.to_numpy()
                p_flatten = pvals_np.flatten()

                # Apply the adjustment method of your choice, for example, Bonferroni
                adj_p_vals = ssm.multipletests(p_flatten, alpha=alpha, method=method)[1]

                # Reshape the adjusted p-values back into the original matrix form
                adj_df = adj_p_vals.reshape(pvals.shape)

                # convert the np array to dataframe and add back the index and column names
                adj_df = pd.DataFrame(adj_df, index=pvals.index, columns=pvals.columns)

                if outputname:
                    adj_df.to_csv(f'{outputname}_{method}_alpha{alpha}_adjusted_pval.csv')
                return adj_df
            
        else:
            print('The input pvals is not a symmetric matrix. The total matrix is used for multiple comparison correction.')
            # convert and faltten the dataframe
            pvals_np = pvals.to_numpy()
            p_flatten = pvals_np.flatten()

            # Apply the adjustment method of your choice, for example, Bonferroni
            adj_p_vals = ssm.multipletests(p_flatten, alpha=alpha, method=method)[1]

            # Reshape the adjusted p-values back into the original matrix form
            adj_df = adj_p_vals.reshape(pvals.shape)

            # convert the np array to dataframe and add back the index and column names
            adj_df = pd.DataFrame(adj_df, index=pvals.index, columns=pvals.columns)

            if outputname:
                adj_df.to_csv(f'{outputname}_{method}_alpha{alpha}_adjusted_pval.csv')
            return adj_df


                

        

        


