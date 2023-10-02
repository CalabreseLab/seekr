#!/usr/bin/env bash

set -e

TEST_DIR="$HOME/Desktop/test_env/"
rm -rf $TEST_DIR
mkdir $TEST_DIR
cd $TEST_DIR
virtualenv venv
source venv/bin/activate
pip_loc=`which pip`
echo "pip: $pip_loc"
pip install seekr
pip show seekr
seekr_download_gencode lncRNA -r 22
seekr_canonical_gencode v22_lncRNA.fa v22-01.fa
seekr_norm_vectors v22-01.fa
head -n 200 v22-01.fa > v22_head200.fa
seekr_kmer_counts v22_head200.fa -o 6mers.npy -mv mean.npy -sv std.npy -b -rl
seekr_pearson 6mers.npy 6mers.npy -o example_vs_self.npy -bi -bo
seekr_visualize_distro example_vs_self.npy example_vs_self.pdf
echo "Finding communities."
seekr_graph example_vs_self.npy .1 -s 0 -g example_vs_self.gml -c communities.csv
cat communities.csv

seekr_find_dist default -fm
seekr_find_pval seqs1.fa seqs2.fa bkg_mean_4mer.npy bkg_std_4mer.npy 4 test_fitres.csv
seekr_kmer_heatmap test_pval.csv 0 1
seekr_kmer_dendrogram test_pval.csv
seekr_kmer_leiden ldseq.fa bkg_mean_4mer.npy bkg_std_4mer.npy 4 -s
seekr_kmer_count_barplot test.fa bkg_mean_4mer.npy bkg_std_4mer.npy 4
seekr_kmer_msd_barplot test.fa bkg_mean_4mer.npy bkg_std_4mer.npy 4 
seekr_kmer_textplot seq1.fa seq2.fa 'ATTA,AAAA,ACTC'

printf "\nCompleted successfully!"
