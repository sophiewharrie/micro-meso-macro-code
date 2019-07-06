# Copyright (C) 2019 by
# Sophie Wharrie <sophie_wharrie@outlook.com>
# The University of Sydney, School of Mathematics and Statistics
"""
Network generators for the triadic closure and configuration models.

The Network class creates networks from these generator functions.
"""
import networkx as nx
import random as rand
from random import randint
import numpy as np
import math

def triadic_closure(p,N):
    """
    Generator function for triadic closure model networks.

    Networks are generated with average degree = 4.

    Model reference: adapted from the basic model described by G. Bianconi, R.K. Darst, 
    J. Iacovacci, S. Fortunato, Triadic closure as a basic generating mechanism of 
    communities in complex networks, Physical Review E, 90, 1 (2014).

    Parameters
    ----------
    p : the probability of triadic closure
    N : the network size
    """
    m = 2 # can be changed, but we always use this value in the paper

    # start with a small connected network of n0 >= m nodes and m0 >= m links
    n0 = 1
    while not(2*m <= n0*(n0-1) and n0 >= m): n0 = n0 + 1
    p0 = 0.3
    G = nx.gnp_random_graph(n0,p0)
    while(not (G.number_of_edges() >= m and nx.is_connected(G))): G = nx.gnp_random_graph(n0,p0)
    
    # at each time step a new node is added to the network with m links
    tmax = N-n0+1
    for t in range(tmax-1):
        # attach the first of the m links to a random node of the network
        new_node = n0+t
        link1 = randint(0,new_node-1)
        G.add_edge(new_node,link1)

        # attach remaining links
        for j in range(1,m):
            node_neighbours = [x for x in G.neighbors(new_node)]
            link_neighbours = list(set().union(*[[x for x in G.neighbors(i)] for i in node_neighbours])) # list of neighbours of new node neighbours
            link_neighbours.remove(new_node) # don't count new_node in neighbours list

            choice = np.random.binomial(1,p)
            # with probability p, attach the jth link to a node chosen randomly among the existing neighbours of the new node
            if choice and [i for i in link_neighbours if i not in node_neighbours]: # it is possible that new_node is already connected to all link_neighbours - skip to else case
                linkj = rand.sample([i for i in link_neighbours if i not in node_neighbours],1)[0]
            # with probability 1-p, attach the jth link to a random node in the network (not already linked)
            else:
                linkj = rand.sample([i for i in list(range(new_node)) if i not in node_neighbours],1)[0]

            G.add_edge(new_node,linkj)

    return(G)

def configuration(c,N):
    """
    Generator function for configuration model networks.

    Networks are generated with average degree = 4.
    The maximum attainable clustering c=0.2 is determined by the average degree.

    Model reference: adapted from the model described by M.E.J. Newman, Random graphs 
    with clustering, Physical Review Letters, 103, 058701 (2009).

    Parameters
    ----------
    c: the (global) clustering coefficient
    N : the network size
    """
    # sample joint degree sequence (s,t) from doubly Poisson distribution,
    # where t = number of triangles and s = number of independent links

    k = 4 # the average degree used in all trials
    s_avg = k*(c*k+c-1)/(c-1)
    t_avg = -c*(k**2)/(2*(c-1))
    s_max = 50
    t_max = 50 # max should be ok because probabilities are very small by this stage

    d = {e:(0,0) for e in range(s_max*t_max)}
    prob = [0]*(s_max*t_max) # the idea is to assign a probability to each (s,t) pair, then sample with these probabilities
    e = 0
    for s in range(s_max):
        for t in range(t_max):
            d[e] = (s, t)
            prob[e] = (math.exp(-s_avg)*((s_avg**s)/math.factorial(s)))*(math.exp(-t_avg)*((t_avg**t)/math.factorial(t)))
            e += 1

    # sample degree distribution from probabilities computed above
    options = list(range(len(prob)))
    indices = np.random.choice(options, N, p=prob)
    s = [d[i][0] for i in indices]
    t = [d[i][1] for i in indices]

    # ensure conditions for multiples of 2 and 3 are met
    while sum(s)%2 != 0:
        random_node = rand.sample(list(range(N)), 1)[0]
        if s[random_node] > 0:
            s[random_node] = s[random_node] - 1
    while sum(t)%3 != 0:
        random_node = rand.sample(list(range(N)), 1)[0]
        if t[random_node] > 0:
            t[random_node] = t[random_node] - 1

    # create network
    joint_degrees = zip(s,t)
    G = nx.random_clustered_graph(joint_degrees)
    G = nx.Graph(G) # remove parallel edges
    G.remove_edges_from(G.selfloop_edges()) # remove self loops
    #G = max(nx.connected_component_subgraphs(G), key=len) # if using largest connected component
    G = fix_graph(G)

    return(G)

def fix_graph(G):
    """
    Fixes a given network so that node labels are numbered consecutively.

    If nodes are not labelled consecutively, issues can arise when converting between
    the different network formats.

    Parameters
    ----------
    G : NetworkX-formatted network
    """
    gfix = nx.Graph()
    edges = G.edges()
    nodes = list(G.nodes())
    num_nodes = len(nodes)
    gfix.add_nodes_from(range(num_nodes))
    for edge in edges:
        gfix.add_edge(nodes.index(edge[0]), nodes.index(edge[1]))
    return(gfix)
