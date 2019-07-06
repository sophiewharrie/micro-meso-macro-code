# Copyright (C) 2019 by
# Sophie Wharrie <sophie_wharrie@outlook.com>
# The University of Sydney, School of Mathematics and Statistics
"""
Implements the four community detection methods (Modularity, Infomap, Spectral and SBM).

Due to the methods originating from various authors and programming languages,
the purpose of this file is to provide a single cohesive Python interface to the four methods.
"""
import networkx as nx
import igraph as ig
from graph_tool import Graph
import graph_tool.inference as gt
import graph_tool.generation as gen
import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath('matlab')
eng.addpath('matlab/spectral_subroutines')
import random
import string

def get_modularity_communities(G):
    """
    Returns the community partition identified by the Modularity method,
    i.e. a list containing the community assignment of each node.

    Reference for Modularity method: A. Clauset, M.E.J. Newman, C. Moore, 
    Finding community structure in very large networks, Physical Review E, 
    70, 6 (2004).

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    Gi = convert_to_igraph(G)
    n = Gi.vcount()
    modularity_communities = Gi.community_fastgreedy().as_clustering()
    modularity_communities = reformat_igraph_output(modularity_communities,n)
    return(modularity_communities)

def get_infomap_communities(G):
    """
    Returns the community partition identified by the Infomap method,
    i.e. a list containing the community assignment of each node.

    Reference for Infomap method: M. Rosvall, D. Axelsson, C.T. Bergstrom, 
    The map equation, EPJ Special Topics, 178, 13-23 (2009).

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    Gi = convert_to_igraph(G)
    n = Gi.vcount()
    infomap_communities = Gi.community_infomap()
    infomap_communities = reformat_igraph_output(infomap_communities,n)
    return(infomap_communities)

def get_spectral_communities(G):
    """
    Returns the community partition identified by the Spectral method,
    i.e. a list containing the community assignment of each node.

    Note: the MATLAB code for the Spectral method takes a GML file as input, 
    so this function creates a file with the network saved in this format.

    Reference for Spectral method: A. Saade, F. Krzakala, L. Zdeborová, 
    Spectral clustering of graphs with the Bethe Hessian, Advances in Neural
    Information Processing Systems, 406-414 (2014).

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    filename = 'temp.gml'
    nx.write_gml(G, filename)
    spectral_communities = eng.spectral_method(filename, nargout=1)
    spectral_communities = [int(eng.single(x)) for x in spectral_communities] # convert to list format
    return(spectral_communities)

def get_sbm_communities(G):
    """
    Returns the community partition identified by the SBM method,
    i.e. a list containing the community assignment of each node.

    Note: the algorithm runs multiple trials due to the stochastic nature of the method
    (selecting the partition with the minimum description length).
    The number of trials can be varied in the code, if required. 

    Reference for SBM method: T.P. Peixoto, Efficient Monte Carlo and greedy heuristic 
    for the inference of stochastic block models, Physical Review E, 89 (2014).

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    trials = 3 # the number of trials to run of the SBM algorithm
    Gto = convert_to_gt(G)
    n = nx.number_of_nodes(G)
    desc_len = [0]*trials
    partitioning = [0]*trials

    # run multiple trials and select the partition with the minimum description length
    for trial in range(trials):
        partitioning_trial = gt.minimize_blockmodel_dl(Gto,B_min=1,B_max=n)
        desc_len[trial] = partitioning_trial.entropy() # description length of a fit (negative log-likelihood)
        partitioning[trial] = partitioning_trial
    
    mn,idx = min((desc_len[i],i) for i in range(0,len(desc_len)))
    sbm_communities = partitioning[idx] # the partitioning with the min description length
    sbm_communities = [x+1 for x in list(sbm_communities.get_blocks().get_array())] # convert to list format
    return(sbm_communities)

def convert_to_igraph(G):
    """
    Takes a NetworkX-formatted network and returns an equivalent igraph-formatted network.

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    n = nx.number_of_nodes(G)
    edges = G.edges()
    Gi = ig.Graph()
    Gi.add_vertices(n)
    Gi.add_edges(edges)
    return(Gi)

def convert_to_gt(G):
    """
    Takes a NetworkX-formatted network and returns an equivalent graphtool-formatted network.

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    n = nx.number_of_nodes(G)
    Gto = Graph(directed=False)
    Gto.add_vertex(n)
    edges = [edge for edge in G.edges()]
    for edge in edges:
        Gto.add_edge(Gto.vertex(edge[0]),Gto.vertex(edge[1]))
    return(Gto)

def reformat_igraph_output(igraph_output, n):
    """
    Takes output from igraph community detection methods and returns the data reformatted in a Python list,
    i.e. a list containing the community assignment of each node.

    Parameters
    ----------
    igraph_output: the community detection output from an igraph function
    n: number of nodes in the network
    """
    numcom = len(igraph_output)
    community_list = [0]*n
    for community in range(numcom):
        for entry in igraph_output[community]: 
            community_list[entry] = community+1
    return(community_list)

def get_number_communities(communities):
    """
    Returns the number of communities detected for a given community partition.

    Parameters
    ----------
    communities : the community partition (as a list)
    """
    return(int(max(communities)))

def get_similarity(communities1, communities2):
    """
    Returns the adjusted mutual information (AMI) between two community partitions,
    indicating the degree of similarity between them.
    
    Reference for AMI code: N.X. Vinh, J. Epps, J. Bailey, Information theoretic 
    measures for clusterings comparison: Variants, properties, normalization and 
    correction for chance, Journal of Machine Learning Research, 11, 2837-2854 (2010).

    Parameters
    ----------
    communities1 : the first set of communities (as a list)
    communities2: the second set of communities to compare (as a list)
    """
    # convert to MATLAB format
    matlab_communities1 = matlab.double(communities1)
    matlab_communities2 = matlab.double(communities2)

    # compute AMI
    ami = eng.ami(matlab_communities1,matlab_communities2)
    return(ami)
