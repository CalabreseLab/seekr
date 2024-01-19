###########################################################################################
### Description:
# plot Leiden community network for input fasta seqeunces
# it works better with small number of nodes (less than 50)
# for better visualization and customization, please set csvfile to a customized name (string) and use Gephi with the saved nodes and edges files

### Details:
# take as input a fasta file with multiple sequences
# calculte sequences distance matrix as seekr pearson correlation of the fasta file to itself
# use Leiden community to call cluster
# edge is weighted with the pearson correlation results
# darker and thicker edges have higher weights
# layout the nodes with networkx spring layout
# the network is undirected

### input:
# inputfile: the route to your input fasta file, please ensure all input sequences have unique headers, this will be used for node label
# mean: the route to your normalization mean vector file
# std: the route to your normalization std vector file
# k: the kmer size used in seekr, which should be consistent with the mean and std file
# algo: the Leiden algorithm used, default is 'RBERVertexPartition'
# other option of algo: 'ModularityVertexPartition', 'RBConfigurationVertexPartition','RBERVertexPartition', 'CPMVertexPartition', 'SurpriseVertexPartition','SignificanceVertexPartition'
# rs: the resolution parameter used in Leiden algorithm, default is 1.0, larger rs leads to more clusters and smaller rs leads to less clusters
# pearsoncutoff: set seekr pearson r value (from seekr.pearson) to 0 if it is less than the pearsoncutoff
# pearsoncutoff aids in removing weak correlation when constructuing Leiden community from a large number of nodes
# pearsoncutoff default is 0, which only removes the negative correlation as a basic requirement for Leiden partition algorithm
# setseed: if set to True, community detection results would be reproducible. 
# setseed default is False, in which case, the results would not be reproducible when the same code run multiple times. This can be used to test the robustness of the detected communities.
# edgecolormethod: the method used to set edge color, default is 'gradient' -- gradient scale of grey for egdes based on edge weights (pearson correlation), darker edges are higher weights
# edgecolormethod: the other option: 'threshold' -- edge color thresholded on edge weights (pearson correlation), smaller than threshold is grey, larger than threshold is black
# edgecolormethod: if other method is passed to edgecolormethod, it will do 'gradient' as the default
# edgethreshold: the threshold used in edgecolormethod 'threshold', default is 0.1
# labelfontsize: the font size of the node label, default is 12
# plotname: the name and path of your output plot, default is None, which suppress the plotting 
# plotname example: if given a string, such as 'test_kmer_leiden', the trailing part .pdf will be added automatically to ensure the pdf format
# plotname can be set to another folder: '/Users/username/Desktop/test_kmer_leiden'
# csvfile: the name and path of the nodes and edges csv files. Default is None, which suppress saving these files. 
# csvfile example: 'testdata', output csv files will be added a trailing part _edges_leiden.csv or _nodes_leiden.csv for edges and nodes, respectively
# The saved edges and nodes csv files are readily importable to Gephi for a better network visualization

### output:
# if plotname is given, a pdf file (plotname.pdf) of the Leiden community network
# if csvfile is set to a string, the corresponding edges and nodes csv files (csvfile_edges_leiden.csv and csvfile_nodes_leiden.csv) will also be saved, which can be readily imported to Gephi

### Example:
# kmer_leiden(inputfile='testld.fa', mean='bkg_mean_4mers.npy',
#             std='bkg_std_4mers.npy', k=4, algo='RBERVertexPartition', rs=1.0, 
#             pearsoncutoff=0.05, setseed=False,
#             edgecolormethod='gradient', edgethreshold=0.1, labelfontsize=12,
#             plotname='test_kmer_leiden', csvfile=None)

###################################################################################

import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg
import matplotlib.pyplot as plt


from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

def kmer_leiden (inputfile, mean, std, k, algo='RBERVertexPartition', rs=1.0, 
                 pearsoncutoff=0, setseed=False,
                 edgecolormethod='gradient', edgethreshold=0.1, labelfontsize=12,
                 plotname=None, csvfile=None):
    
    # first check whether input k_mer is compatible with the mean and std files
    meanfile = np.load(mean)
    stdfile = np.load(std)

    if len(meanfile) != 4**(k) | len(stdfile) != 4**(k):
        print('kmer size is not compatible with the normalization mean and/or std files.')
        print('Please make sure the normalization mean and std files are generated using the same kmer size as specified here in k.')
        print('No Leiden community is calculated or plotted. The output is None.')
        return None

    t1 = seekrBasicCounter(inputfile,
                        mean=mean,
                        std=std,
                        k=k,silent=True)
    t1.make_count_file()
    header1=seekrReader(inputfile).get_headers()
    # remove the '>' in the header
    header1 = [i[1:] for i in header1]

    # calculate pearson correlation
    ld_sim = seekrPearson(t1.counts,t1.counts)

    # set all pearson r value less than pearsoncutoff to 0 in ld_sim
    ld_sim[ld_sim<pearsoncutoff] = 0
    # set diagonal values to 0
    np.fill_diagonal(ld_sim, 0)

    # convert to pandas dataframe with row and column names as headers
    df = pd.DataFrame(ld_sim, columns=header1, index=header1)

    # Create a graph from the correlation matrix
    G = nx.from_pandas_adjacency(df)

    # Convert the networkx graph to an igraph graph
    # Note the "mode" argument is set to 'UNDIRECTED'
    ig_graph = ig.Graph.Adjacency((df.values > 0).tolist(),mode='UNDIRECTED')
    ig_graph.es['weight'] = df.values[df.values > 0].flatten()

    # print(ig_graph.es.attributes()) # check edge attributes- weight
    # print(ig_graph.ecount()) # check edge count

    # Define the resolution parameter: rs
    # default rs= 1.0, Adjust this to your needs
    # Map algorithm names to classes
    algo_dict = {
        'ModularityVertexPartition': leidenalg.ModularityVertexPartition,
        'RBConfigurationVertexPartition': leidenalg.RBConfigurationVertexPartition,
        'RBERVertexPartition': leidenalg.RBERVertexPartition,
        'CPMVertexPartition': leidenalg.CPMVertexPartition,
        'SurpriseVertexPartition': leidenalg.SurpriseVertexPartition,
        'SignificanceVertexPartition': leidenalg.SignificanceVertexPartition
        }

    # Run the Leiden algorithm 
    if setseed==True:
        kseed=1
    else:
        kseed=None


    if algo=='SignificanceVertexPartition':
      partition = leidenalg.find_partition(ig_graph,
                                             algo_dict[algo],
                                             seed=kseed)
      
    elif algo=='SurpriseVertexPartition' or algo=='ModularityVertexPartition':
      partition = leidenalg.find_partition(ig_graph,
                                           algo_dict[algo],
                                           weights='weight',
                                           seed=kseed)
    else:
        partition = leidenalg.find_partition(ig_graph,
                                             algo_dict[algo], 
                                             weights='weight',
                                             resolution_parameter=rs,
                                             seed=kseed)

    # Get a list of node names
    node_names = list(df.index)

    # check plotname to see if plotting is needed
    if plotname:

        if edgecolormethod == 'gradient':
            # gradient scale of grey for egdes based on edge weights (pearson correlation)
            # darker edges are higher weights
            # pearson correlation of 1 is colored black

            # Get upper triangular indices (excluding diagonal)
            row_indices, col_indices = np.triu_indices(df.shape[0], k=1)  # k=1 excludes diagonal

            # Extract values from the upper triangle of the matrix
            df_up = df.values[row_indices, col_indices]

            # Filter out zeros to get non-zero upper triangular values
            df_cl = df_up[df_up > 0]

            # Map normalized weights to [0.1, 1] to avoid totally transparent edges
            weights_normalized = (df_cl - df_cl.min()) / (df_cl.max() - df_cl.min())
            weights_mapped =0.1 + 0.9 * weights_normalized

            # Convert to shades of grey, where higher weights are darker
            edge_colors = [(1-w, 1-w, 1-w) for w in weights_mapped.flatten()]

            # Create edge widths based on the weights
            edge_widths = [(1 +  3*w) for w in weights_mapped.flatten()]

            # # Create dictionaries to map edges to their colors, widths, and original weights
            # edge_color_map = {}
            # edge_width_map = {}
            # edge_weight_map = {}
            # for i, edge in enumerate(G.edges()):
            #     edge_color_map[edge] = edge_colors[i]
            #     edge_width_map[edge] = edge_widths[i]
            #     edge_weight_map[edge] = G[edge[0]][edge[1]]['weight']
            # # Now, let's test:
            # edge_example = list(G.edges())[0]
            # print(f"Edge: {edge_example}")
            # print(f"Weight: {edge_weight_map[edge_example]}")
            # print(f"Color: {edge_color_map[edge_example]}")
            # print(f"Width: {edge_width_map[edge_example]}")

            # Assign community colors to nodes
            community_colors = plt.cm.rainbow(np.linspace(0, 1, max(partition.membership) + 1))
            node_colors = [community_colors[comm] for comm in partition.membership]

            # Use the spring layout (Fruchterman-Reingold force-directed algorithm)
            pos = nx.spring_layout(G, weight='weight')

            # Plot the graph using networkx and matplotlib
            plt.figure(figsize=(15, 15))
            # Turn off the axis
            plt.gca().axis('off')
            #plt.subplots_adjust(left=0.9, right=0.9, top=0.1, bottom=0.1)  # adjust margins
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500)
            nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths)
            nx.draw_networkx_labels(G, pos, font_size=labelfontsize, font_family="sans-serif")

            plt.tight_layout()

            # Save the figure
            plt.savefig(f'{plotname}.pdf')

            # # Plot the graph with the edge colors and edge widths based on weight
            # ig.plot(partition, vertex_label=node_names,
            #         layout=layout, edge_curved=False, margin=margin,
            #         edge_color=edge_colors,edge_width=edge_widths,
            #         target=f'{outputname}.pdf')
            
        elif edgecolormethod == 'threshold':
            # threshold for edge color based on edge weights (pearson correlation)
            # edge weight below threshold is colored grey
            # edge weight above threshold is colored black

            # Get upper triangular indices (excluding diagonal)
            row_indices, col_indices = np.triu_indices(df.shape[0], k=1)  # k=1 excludes diagonal

            # Extract values from the upper triangle of the matrix
            df_up = df.values[row_indices, col_indices]

            # Filter out zeros to get non-zero upper triangular values
            df_cl = df_up[df_up > 0]

            # Define the threshold
            threshold = edgethreshold

            # Assign colors based on whether the weight is above or below the threshold
            edge_colors = ["black" if w > threshold else "grey" for w in df_cl]

            # Create edge widths based on the weights
            edge_widths = [4 if w > threshold else 1 for w in df_cl]

            # Assign community colors to nodes
            community_colors = plt.cm.rainbow(np.linspace(0, 1, max(partition.membership) + 1))
            node_colors = [community_colors[comm] for comm in partition.membership]

            # Use the spring layout (Fruchterman-Reingold force-directed algorithm)
            pos = nx.spring_layout(G, weight='weight')

            # Plot the graph using networkx and matplotlib
            plt.figure(figsize=(15, 15))
            # Turn off the axis
            plt.gca().axis('off')
            #plt.subplots_adjust(left=0.9, right=0.9, top=0.1, bottom=0.1)  # adjust margins
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500)
            nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths)
            nx.draw_networkx_labels(G, pos, font_size=labelfontsize, font_family="sans-serif")

            plt.tight_layout()

            # Save the figure
            plt.savefig(f'{plotname}.pdf')

            # # Plot the graph with the edge colors based on the threshold
            # ig.plot(partition, vertex_label=node_names,
            #         layout=layout, edge_curved=False, margin=margin,
            #         edge_color=edge_colors,edge_width=edge_widths,
            #         target=f'{outputname}.pdf')
            
        else:
            # if input does not match 'gradient' or 'threshold', use default 'gradient'
            print("edgecolormethod must be either 'gradient' or 'threshold', use default 'gradient' now")
            
            # gradient scale of grey for egdes based on edge weights (pearson correlation)
            # darker edges are higher weights
            # pearson correlation of 1 is colored black

            # Get upper triangular indices (excluding diagonal)
            row_indices, col_indices = np.triu_indices(df.shape[0], k=1)  # k=1 excludes diagonal

            # Extract values from the upper triangle of the matrix
            df_up = df.values[row_indices, col_indices]

            # Filter out zeros to get non-zero upper triangular values
            df_cl = df_up[df_up > 0]

            # Map normalized weights to [0.1, 1] to avoid totally transparent edges
            weights_normalized = (df_cl - df_cl.min()) / (df_cl.max() - df_cl.min())
            weights_mapped =0.1 + 0.9 * weights_normalized

            # Convert to shades of grey, where higher weights are darker
            edge_colors = [(1-w, 1-w, 1-w) for w in weights_mapped.flatten()]

            # Create edge widths based on the weights
            edge_widths = [(1 +  3*w) for w in weights_mapped.flatten()]

            # Assign community colors to nodes
            community_colors = plt.cm.rainbow(np.linspace(0, 1, max(partition.membership) + 1))
            node_colors = [community_colors[comm] for comm in partition.membership]

            # Use the spring layout (Fruchterman-Reingold force-directed algorithm)
            pos = nx.spring_layout(G, weight='weight')

            # Plot the graph using networkx and matplotlib
            plt.figure(figsize=(15, 15))
            # Turn off the axis
            plt.gca().axis('off')
            #plt.subplots_adjust(left=0.9, right=0.9, top=0.1, bottom=0.1)  # adjust margins
            nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=500)
            nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths)
            nx.draw_networkx_labels(G, pos, font_size=labelfontsize, font_family="sans-serif")

            plt.tight_layout()

            # Save the figure
            plt.savefig(f'{plotname}.pdf')
            

    if csvfile: 
        # Initialize lists to hold labels and colors
        labels = []
        colors = []

        # Iterate over each community
        for i, community in enumerate(partition):
            # For each node in the community, append the node label and community number
            for node_index in community:
                labels.append(node_names[node_index])
                colors.append(i + 1)

        # Create a DataFrame
        df_out = pd.DataFrame({'Id': labels, 'Label': labels, 'Color': colors})

        # Write nodes to CSV
        df_out.to_csv(f'{csvfile}_nodes_leiden.csv', index=False)

        # Get edges
        # Extract the upper triangle (excluding the diagonal)
        mask = np.triu(np.ones(df.shape), k=1).astype(bool)

        # Filter dataframe values using the mask and then melt
        df_triu = df.where(mask).stack().reset_index()
        df_triu.columns = ['Source', 'Target', 'Weight']

        # Write edges to CSV
        df_triu.to_csv(f'{csvfile}_edges_leiden.csv', index=False)

