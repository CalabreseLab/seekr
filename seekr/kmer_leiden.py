###########################################################################################
### Description:
# plot Leiden community network for input fasta seqeunces
# it works better with small number of nodes (less than 50)
# for better visualization and customization, please set savecsv to True and use Gephi with the saved nodes and edges files

### Details:
# take as input a fasta file with multiple sequences
# calculte sequences distance matrix as seekr pearson correlation of the fasta file to itself
# use Leiden community to call cluster
# edge is weighted with the pearson correlation results
# layout the nodes with igraph
# the network is undirected

### input:
# inputfile: the route to your input fasta file, please ensure all input sequences have unique headers, this will be used for node label
# mean: the route to your normalization mean vector file
# std: the route to your normalization std vector file
# k: the kmer size used in seekr, which should be consistent with the mean and std file
# algo: the Leiden algorithm used, default is 'RBERVertexPartition'
# other option of algo: 'ModularityVertexPartition', 'RBConfigurationVertexPartition','RBERVertexPartition', 'CPMVertexPartition', 'SurpriseVertexPartition','SignificanceVertexPartition'
# rs: the resolution parameter used in Leiden algorithm, default is 1.0, larger rs leads to more clusters and smaller rs leads to less clusters
# edgecolormethod: the method used to set edge color, default is 'gradient' -- gradient scale of grey for egdes based on edge weights (pearson correlation)
# edgecolormethod: the other option: 'threshold' -- edge color thresholded on edge weights (pearson correlation), smaller than threshold is grey, larger than threshold is black
# edgecolormethod: if other method is passed to edgecolormethod, it will do 'gradient' as the default
# edgethreshold: the threshold used in edgecolormethod 'threshold', default is 0.1
# margin: the margin used in igraph plot, default is 100, set it to be large so that the node label will not be cut off
# outputname: the name and path of your output file, default is 'test', this will broadcast to all saved images (pdf) and csv files
# output images will be added a trailing part _gradient_leiden_network.pdf or _threshold_leiden_network.pdf based on the edgecolormethod
# output csv files will be added a trailing part _edges_leiden.csv or _nodes_leiden.csv for edges and nodes, respectively
# outputname other example: '/Users/username/Desktop/test'
# savecsv: whether to save the nodes and edges csv files, default is True. The saved edges and nodes csv files are readily importable to Gephi for a better network visualization

### output:
# a pdf file (outputname_gradient_leiden_network.pdf or outputname_threshold_leiden_network.pdf) of the Leiden community network
# if savecsv is True, the corresponding edges and nodes csv files (outputname_edges_leiden.csv and outputname_nodes_leiden.csv) will also be saved, which can be readily imported to Gephi

### Example:
# kmer_leiden(inputfile='testld.fa', mean='gencodevM25_unique_mean_4mers.npy',
#             std='gencodevM25_unique_std_4mers.npy', k=4, algo='RBERVertexPartition', rs=1.0,
#             edgecolormethod='gradient', edgethreshold=0.1, margin=100,
#             outputname='test', savecsv=True)

###################################################################################

import pandas as pd
import numpy as np
import networkx as nx
import igraph as ig
import leidenalg


from seekr.kmer_counts import BasicCounter as seekrBasicCounter
from seekr.pearson import pearson as seekrPearson
from seekr.fasta_reader import Reader as seekrReader

def kmer_leiden (inputfile, mean, std, k, algo='RBERVertexPartition', rs=1.0, 
                 edgecolormethod='gradient', edgethreshold=0.1, margin=100,
                 outputname='test', savecsv=True):

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

    # set all negative values to 0 in test_sim
    ld_sim[ld_sim<0] = 0
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
    if algo=='SurpriseVertexPartition':
        partition = leidenalg.find_partition(ig_graph,
                                             algo_dict[algo],
                                             weights='weight')
    else:
        partition = leidenalg.find_partition(ig_graph,
                                             algo_dict[algo], 
                                             weights='weight',
                                             resolution_parameter=rs)

    # Get a list of node names
    node_names = list(df.index)

    if edgecolormethod == 'gradient':
        # gradient scale of grey for egdes based on edge weights (pearson correlation)
        # darker edges are higher weights
        # pearson correlation of 1 is colored black

        # Normalize the weights to the range [0, 1]
        weights_normalized = (ig_graph.es['weight'] - min(ig_graph.es['weight'])) / (max(ig_graph.es['weight']) - min(ig_graph.es['weight']))

        # Map normalized weights to [0.1, 1] to avoid totally transparent edges
        weights_mapped = 0.1 + 0.9 * weights_normalized

        # Convert to shades of grey, where higher weights are darker
        edge_colors = [(1 - w, 1 - w, 1 - w) for w in weights_mapped]

        # Create edge widths based on the weights
        edge_widths = [1 + 3 * w for w in weights_mapped]

        # Use the Fruchterman-Reingold force-directed algorithm, considering edge weight
        layout = ig_graph.layout('fr', weights='weight')

        # Plot the graph with the edge colors and edge widths based on weight
        ig.plot(partition, vertex_label=node_names,
                layout=layout, edge_curved=False, margin=margin,
                edge_color=edge_colors,edge_width=edge_widths,
                target=f'{outputname}_gradient_leiden_network.pdf')
        
    elif edgecolormethod == 'threshold':
        # threshold for edge color based on edge weights (pearson correlation)
        # edge weight below threshold is colored grey
        # edge weight above threshold is colored black

        # Define the threshold
        threshold = edgethreshold

        # Assign colors based on whether the weight is above or below the threshold
        edge_colors = ["black" if w > threshold else "grey" for w in ig_graph.es['weight']]

        # Create edge widths based on the weights
        edge_widths = [4 if w > threshold else 1 for w in ig_graph.es['weight']]

        # Use the Fruchterman-Reingold force-directed algorithm, considering edge weight
        layout = ig_graph.layout('fr', weights='weight')

        # Plot the graph with the edge colors based on the threshold
        ig.plot(partition, vertex_label=node_names,
                layout=layout, edge_curved=False, margin=margin,
                edge_color=edge_colors,edge_width=edge_widths,
                target=f'{outputname}_threshold_leiden_network.pdf')
        
    else:
        # if input does not match 'gradient' or 'threshold', use default 'gradient'
        print("edgecolormethod must be either 'gradient' or 'threshold', use default 'gradient' now")
        
        # gradient scale of grey for egdes based on edge weights (pearson correlation)
        # darker edges are higher weights
        # pearson correlation of 1 is colored black

        # Normalize the weights to the range [0, 1]
        weights_normalized = (ig_graph.es['weight'] - min(ig_graph.es['weight'])) / (max(ig_graph.es['weight']) - min(ig_graph.es['weight']))

        # Map normalized weights to [0.1, 1] to avoid totally transparent edges
        weights_mapped = 0.1 + 0.9 * weights_normalized

        # Convert to shades of grey, where higher weights are darker
        edge_colors = [(1 - w, 1 - w, 1 - w) for w in weights_mapped]

        # Create edge widths based on the weights
        edge_widths = [1 + 3 * w for w in weights_mapped]

        # Use the Fruchterman-Reingold force-directed algorithm, considering edge weight
        layout = ig_graph.layout('fr', weights='weight')

        # Plot the graph with the edge colors and edge widths based on weight
        ig.plot(partition, vertex_label=node_names,
                layout=layout, edge_curved=False, margin=margin,
                edge_color=edge_colors,edge_width=edge_widths,
                target=f'{outputname}_gradient_leiden_network.pdf')
        

    if savecsv == True: 
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
        df_out.to_csv(f'{outputname}_nodes_leiden.csv', index=False)

        # Get edges
        # Extract the upper triangle (excluding the diagonal)
        mask = np.triu(np.ones(df.shape), k=1).astype(bool)

        # Filter dataframe values using the mask and then melt
        df_triu = df.where(mask).stack().reset_index()
        df_triu.columns = ['Source', 'Target', 'Weight']

        # Write edges to CSV
        df_triu.to_csv(f'{outputname}_edges_leiden.csv', index=False)



