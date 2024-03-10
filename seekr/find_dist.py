####################################################################################################################
### Description: 
# Find the best fitted distribution to the input background sequences 

### Details:
# data to be fitted -- all possible pairwise kmer pearson correlation for the background sequences
# can also choose to return the actual data by itself
# the results (best fitted distribution or actual data) will be used to calculate p-values

### Input:
# inputseq: full path to the fasta file containing background sequences. By default (inputseq='default') use mouse vM25 Long non-coding RNA unqiue transcript sequences
# 'default' lncRNA sequences were downloaded from here: https://www.gencodegenes.org/mouse/release_M25.html and duplicated sequences were removed
# There are 18856 total transcripts in gencode lncRNA. After removing duplicated sequences, there are 12996 unique sequences.
# k_mer: k-mer size, by default k_mer=4
# log2: whether to do log2 transform, possible options are: 'Log2.post','Log2.pre' and 'Log2.none', default is 'Log2.post'
# Log2.post -- Log2 transformation post-standardization: self.counts += np.abs(np.min(self.counts)) self.counts += 1 self.counts = np.log2(self.counts)
# Log2.pre -- Log2 transformation pre-standardization: self.counts += np.abs(np.min(self.counts)) self.counts += 1
# Log2.none -- no log2 transformation
# models: groups of candidate models for fitting, options are: 'all','common10' (default), or a list of distributions by user input for example ['norm','expon','pareto']
# 'all' - fit all available distributions from scipy stats, both rv_continuous and rv_discrete: https://docs.scipy.org/doc/scipy/reference/stats.html#continuous-distributions
# 'common10' - fit 10 common distributions ('cauchy', 'chi2', 'expon', 'exponpow', 'gamma', 'lognorm', 'norm', 'pareto', 'rayleigh', 'uniform')
# common10 distributions are deteremined from get_common_distributions() from the Fitter library with 'powerlaw' replaced by 'pareto' as pareto belongs to powerlaw family and bears a biology relation 
# user can also choose to input their own list of distributions in a list format, but be sure all distributions are available in scipy.stats
# subsetting: True (default) or False. whether or not to use a subset of the data for fitting or output. For large datasets, subsetting is recommended to save time
# here subsetting is performed after getting all possible pairwise seekr pearson correlation values for the background sequences
# subsetting is not performed on the input background sequences
# subset_size: the size of the subset to use for fitting or output. Default is 100000. only be used when subsetting=True
# if subset_size is larger than the actual data size, the actual data size will be used instead
# fit_model: True (default) or False. whether or not to fit the data to the distributions sepcified in models. 
# For small dataset, you can set fit_model=False, this will return the actual data as a npy array, without fitting to any distributions
# in this way, the p values will be calculated based on the actual data instead of the fitted distribution
# statsmethod: the stats method used to quantify the goodness of fit of the distributions vs data
# options are: 'ks' (default, Kolmogorov-Smirnov test), 'mse' (Mean Squared Error), 'aic' (Akaike Information Criterion), 'bic' (Bayesian Information Criterion)
# smaller values returned by all four statsmends generally indicate a better fit
# KS: measures the largest difference between the observed cumulative distribution function (CDF) of the data and the expected CDF of the model
# it is sensitive to differences in both the center and tails of the distribution and it provides an absolute measure of goodness-of-fit, which is intuitive and easy to understand
# it can be sensitive to the sample size, and it may not have sufficient power to detect subtle differences in distributions when the sample size is small
# it does not penalize for model complexity, so overfitting could be a potential problem
# MSE: is an error-based metric, calculates the mean squared error between the actual data and the fitted data, it is intuitive and easy to understand
# but it is sensitive to outliers, and it does not take into account the number of parameters in the model
# Both AIC and BIC are likelihood based metrics which has some advantages over error-based metrics:
# they balance the goodness of fit against the complexity of the model; they have broader applicability especially for models where the mean and variance are not well-defined (such as cauchy)
# AIC: Tends to favor models that fit the data very closely. It can sometimes choose overly complex models
# BIC: Incorporates a penalty term for the number of parameters in the model and the sample size, and therefore can result in the selection of simpler models compared to AIC
# but AIC and BIC takes longer time to calculate and the output of AIC and BIC are unitless, they are good for comparing different models rather than providing an interpretable absolute measure of "goodness of fit"
# The KS test measures how well the model fits the entire data distribution. MSE is designed to measure the predictive accuracy of a model but doesn't account for model complexity. AIC and BIC are designed for model selection, balancing fit and complexity.
# MSE is sensitive to outliers, while AIC and BIC can also be influenced by them depending on the likelihood formulation. The KS test is less sensitive to outliers but more focused on the overall shape of the distribution
# If not sure which statsmethod works for the input data, try several of them to get a well-rounded view and set plotfit=True to see the fitted distributions vs the actual data
# progress_bar: True or False (default). whether or not to show the progress bar when fitting the data to the distributions.
# plotfit: None (default, meaning not save the plot) or full path of file to save the plot of the fitted distributions (red dash line) vs the actual data (blue solid line) for all distributions in models.
# plotfit example: plotfit = 'testplot' or '/Users/username/Desktop/testplot', a trailing part .pdf will be added to ensure pdf format
# outputname: full path of the output csv file, default is None, which means not save the output to a csv file
# outputname example: 'test_fitres' or '/Users/username/Desktop/test_fitres', a trailing part .csv will be added to ensure csv format

### Output:
# If fit_model is True (fit ditributions to the data), returns the distributions name, goodness of fit (like D stats in ks) and fitted params for all distributions in models
# the results here is ordered by the goodness of fit, the first distribution (results[0]) is the best fitted distribution (to be sure, use manual inspection by setting plotfit=True)
# If fit_model is False (do not fit distributions to the data), returns the actual data (with or without subsetting based on subsetting argument) as a npy array, without fitting to any distributions
# if plotfit is given, a pdf plot with the given name and path will be saved, showing the fitted distributions (red dash line) vs the actual data (blue solid line) for all distributions in models
# if outputname is given, a csv file with the given name and path will be saved for either a list of fitted distributions or a npy array

### Example:
# fitres = find_dist(inputseq='default', k_mer=4, log2='Log2.post', models='common10', subsetting=True, subset_size = 10000, fit_model=False, statsmethod='ks',progress_bar=True, plotfit=None, outputname='test')

####################################################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import kstest
import warnings
import os

from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson

from tqdm import tqdm



# Finds and returns distributions and respective parameters that best fit input sequences
def find_dist(inputseq='default', k_mer=4, log2='Log2.post', models='common10', subsetting=True, subset_size = 100000, fit_model=True, statsmethod='ks',progress_bar=False, plotfit=None, outputname=None):

    #Determine what background transcriptome to use
    if inputseq == 'default':
        print('Using default background sequences: mouse vM25 Long non-coding RNA transcript sequences from Gencode.')
        print('https://www.gencodegenes.org/mouse/release_M25.html')
        print('Only including unique sequences and have the full length Airn(chr17:12741398-12830151) appended.')
        # Get the directory where the current file is located
        current_dir = os.path.dirname(os.path.realpath(__file__))
        
        # Construct the path to the .fa file
        inputseq = os.path.join(current_dir, 'data', 'gencode.vM25.lncRNA_transcripts.unique.genesequence_withfullairn.fa')
    
    #preprare the list of distributions to be used
    if models == 'common10':
        distributions = ['cauchy', 'chi2', 'expon', 'exponpow', 'gamma', 
                         'lognorm', 'norm', 'pareto', 'rayleigh', 'uniform']
        
    else: 
        #Create all distributions (continous and discrete) that would be used
        continuous_distributions = [d for d in dir(stats) if
                        isinstance(getattr(stats, d), stats.rv_continuous)]

        discrete_distributions = [d for d in dir(stats) if
                                isinstance(getattr(stats, d), stats.rv_discrete)]

        distributions = continuous_distributions + discrete_distributions 

        # it is a convention to start the name of private or special methods with an underscore
        # this line is used to filter out private or special methods from the scipy.stats module
        distributions = [d for d in distributions if d[0] != '_']
        # remove levy_stable and studentized_range, it has problem converging
        # and the shape of the distribution is not what we want
        distributions = [d for d in distributions if d != 'levy_stable']
        distributions = [d for d in distributions if d != 'studentized_range']
        # if 'gilbrat' is in distributions, replace 'gilbrat' with 'gibrat' as 'gilbrat' is a misspelling of 'gibrat'
        # if 'gilbrat' in distributions:
        #     distributions[distributions.index('gilbrat')] = 'gibrat'
        # use what it is from the stats package

        if models != 'all':
            # check if the user input distributions are included in the distributions list
            orilen1=len(models)
            distributions = [d for d in models if d in distributions]
            if len(distributions) < orilen1:
                print("Please enter valid distribution names available in scipy.stats. refer to https://docs.scipy.org/doc/scipy/reference/stats.html#continuous-distributions")
                # find out which distributions are excluded
                excluded1 = [d for d in models if d not in distributions]
                print(f"Excluding invalid distributions for fitting: {excluded1}")
        

    # Convert distribution names to distribution objects
    tempdist = distributions.copy()
    orilen2=len(distributions)
    distributions = [getattr(stats, d) for d in distributions]

    # Filter distributions - some distributions can throw errors during the fit step
    # filter out the distributions that do not have a 'fit' method
    distributions = [d for d in distributions if 'fit' in dir(d)]
    if (len(distributions) < orilen2) and (models != 'all' and models != 'common10') :
        # find out which distributions are excluded
        excluded2 = [d for d in tempdist if d not in distributions]
        print(f"Excluding distributions do not have a 'fit' method: {excluded2}")


    #Create normalization mean and std vector for the inputseq
    bkg_norm_counter = seekrBasicCounter(inputseq,log2=log2,k=k_mer,silent=True)
    bkg_norm_counter.get_counts()
    mean_path = f'bkg_mean_{k_mer}mers.npy'
    std_path = f'bkg_std_{k_mer}mers.npy'
    np.save(mean_path, bkg_norm_counter.mean)
    np.save(std_path, bkg_norm_counter.std)

    # Count k-mers
    bkg_counter = seekrBasicCounter(inputseq,mean=mean_path,std=std_path,k=k_mer,silent=True)
    bkg_counter.make_count_file()

    # Find similarities
    sim_counts = seekrPearson(bkg_counter.counts,bkg_counter.counts)

    # Extract the upper triangle and flatten it
    sim_triu = sim_counts[np.triu_indices(sim_counts.shape[0], k=1)]

    #If the user wants to subset the data
    if subsetting: 
        if len(sim_triu)>subset_size:
            # randomly sample the user defined number of the sim_triu
            sim_triu = np.random.choice(sim_triu, size=subset_size, replace=False)
        else:
            print("subset_size is larger than the actual data size, use the actual data size instead")
    
    

    #if fitting chosen, fit the data to the distributions
    if fit_model:
        
        if len(distributions)>50 and len(sim_triu)>5000000 and subsetting==False:
            print("The input sequence count and distribution number for fitting are both large, subsetting is recommended to save time")
    
        results = []
        # set iterable based on progress_bar
        iterable = tqdm(distributions) if progress_bar else distributions
        # Fit distributions to data
        for distribution in iterable:
            # Try to fit the distribution
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                try: 
                    params = distribution.fit(sim_triu)

                    if statsmethod == 'ks':
                        # Get the KS test result on fitted data
                        D, p = kstest(sim_triu, distribution.name, args=params)

                    elif statsmethod == 'mse':
        
                        if isinstance(distribution, stats.rv_continuous):
                            synthetic_data = distribution.rvs(*params, size=len(sim_triu))
                        else:
                            synthetic_data = distribution.rvs(*params[:-2], loc=params[-2], scale=params[-1], size=len(sim_triu))
                            
                        residuals = sim_triu - synthetic_data
                        D = np.mean(residuals**2)


                    elif statsmethod == 'aic' or statsmethod == 'bic':
                        # Calculate log-likelihood
                        if isinstance(distribution, stats.rv_continuous):
                            log_likelihood = np.sum(distribution.logpdf(sim_triu, *params))
                        else:
                            log_likelihood = np.sum(distribution.logpmf(sim_triu, *params[:-2], loc=params[-2], scale=params[-1]))
                        
                        # Number of parameters
                        k = len(params)
                        # Number of data points
                        n = len(sim_triu)
                        
                        # Calculate AIC and BIC
                        if statsmethod == 'aic':
                            D = 2 * k - 2 * log_likelihood

                        else: 
                            D = np.log(n) * k - 2 * log_likelihood

                    else:
                        print("Please enter a valid statsmethod: 'ks', 'mse', 'aic', or 'bic'. Use default 'ks' now.")
                        D, p = kstest(sim_triu, distribution.name, args=params)

                # catch the error if the distribution cannot be fit, exclude it from the results
                except Exception as e:
                    print(f"Could not fit {distribution.name} because {e}, excluding it from the results")
                    continue

                # save the name and test result
                results.append((distribution.name, D, params))

        # Sort by minimum D statistics
        results.sort(key=lambda x: x[1])
        # sort by p val
        # results.sort(key=lambda x: x[2], reverse=True)

        if plotfit:
            n = len(results)
            # get the num of cols as the minimal of 5 and n
            n_cols = min(5, n)
            n_rows = n // n_cols + (n % n_cols > 0)   # Calculate the number of rows needed

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3, n_rows*3))  # Increase figure size as necessary
            axes = axes.ravel()  # Flatten the axes array

            # Generate PDF for each distribution and plot it
            for idx, (ax, result) in enumerate(zip(axes, results)):
                dist_name, Dval,params = result
                # print(dist_name)
                distribution = getattr(stats, dist_name)
                
                # Generate PDF
                x = np.linspace(min(sim_triu), max(sim_triu), 1000)
                pdf = distribution.pdf(x, *params)

                # Plot the original data as histogram and the fitted model
                ax.hist(sim_triu, bins=100, density=True, alpha=0.6, color='skyblue')
                ax.plot(x, pdf, 'r--', linewidth=2)
                ax.set_title(f'{idx+1}: {dist_name} (Dev={Dval:.3f})')

            # Remove unused subplots
            for i in range(len(results), len(axes)):
                fig.delaxes(axes[i])

            plt.tight_layout()
            # save plot
            plt.savefig(f'{plotfit}.pdf',dpi=300)

        # results is a list of tuples like (distribution_name, D_statistics, params)
        if outputname:
            # convert to dataframe and adding row/column names
            results_df = pd.DataFrame(results, columns=['distribution_name', 'D_statistics', 'params'])
            # save the dataframe to csv file
            results_df.to_csv(f'{outputname}.csv', index=False)

        return results
    else:
        if plotfit:
            print('No plot will be produced as fit_model is set to False, please set fit_model=True to plot the fitted distributions vs the actual data')
        
        if outputname:
            # save sim_triu as a csv file
            np.savetxt(f'{outputname}.csv', sim_triu, delimiter=",")

        return sim_triu

