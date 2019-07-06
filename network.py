# Copyright (C) 2019 by
# Sophie Wharrie <sophie_wharrie@outlook.com>
# The University of Sydney, School of Mathematics and Statistics
"""
Base class for generated networks.

The Network class creates network objects from the given network generators,
and provides methods for examining the community structure.
"""
from network_generators import *
from community_detection import *
import matplotlib.pyplot as plt
import itertools

class Network:
    """
    Base class for generated networks.

    The Network class creates network objects from the given network generators,
    and provides methods for examining the community structure.

    Parameters (required)
    ----------
    model: network generation model (triadic_closure or configuration)
        The network generation model as a function, not as a string
        (e.g. as triadic_closure, not 'triadic_closure').
    N: network size (integer > 0)
        The number of nodes in the generated network as an integer greater
        than zero.
    t: clustering parameter (model dependent)
        For the triadic closure model, 0 <= t <= 1 (acts as the parameter p).
        For the configuration model, 0 <= t <= 0.2 (acts as the parameter c).
        These limits are determined by the model and the fact that we set 
        average degree = 4 across all networks.

    Examples
    ----------
    Create a triadic closure network with 1000 nodes and high clustering (p=1).

    >>> network = Network(model=triadic_closure, N=1000, t=1)

    **Plot:**

    View an image of the generated network. 
    
    The plot is saved in a location determined by the filename parameter.

    >>> filename = 'test_plot.png'
    >>> network.plot(filename)

    **NetworkX:**

    Retrieve the generated network as a NetworkX Graph object.

    This allows you to apply functions from the NetworkX package.

    >>> G = network.graph()

    For example, calculate the (global) clustering coefficient.

    >>> nx.transitivity(G)

    **Communities:**

    Determine communities using the Modularity method.

    >>> modularity_communities = network.get_communities(method='modularity')

    Get the number of communities.

    >>> get_number_communities(modularity_communities)

    Plot the network with nodes coloured by community.

    >>> filename = 'test_plot_colours.png'
    >>> network.plot(filename, communities=modularity_communities)

    Run the summary method that applies all community detection methods.

    >>> network.communities_summary()

    **Similarity:**

    Community similarity is quantified by the adjusted mutual information (AMI).

    Compare communities identified by the Modularity and SBM methods.

    >>> modularity_communities = network.get_communities(method='modularity')
    >>> SBM_communities = network.get_communities(method='sbm')
    >>> get_similarity(modularity_communities, SBM_communities)

    Note that the summary method reports the AMI between all pairs of methods.

    >>> network.communities_summary()
    """
    def __init__(self, model, N, t):
        """
        Initialise a network with given parameters.

        Parameters (required)
        ----------
        model: network generation model (triadic_closure or configuration)
            The network generation model as a function, not as a string
            (e.g. as triadic_closure, not 'triadic_closure').
        N: network size (integer > 0)
            The number of nodes in the generated network as an integer greater
            than zero.
        t: clustering parameter (model dependent)
            For the triadic closure model, 0 <= t <= 1 (acts as the parameter p).
            For the configuration model, 0 <= t <= 0.2 (acts as the parameter c).

        Examples
        ----------
        Create a triadic closure network with 1000 nodes and high clustering (p=1).

        >>> network = Network(model=triadic_closure, N=1000, t=1)

        Create a configuration network with 500 nodes and high clustering (c=0.2).

        >>> network = Network(model=configuration, N=500, t=0.2)

        Create a triadic closure network with 1000 nodes and no clustering (p=0, like 
        a random network).

        >>> network = Network(model=triadic_closure, N=1000, t=0)
        """
        # check input values
        models = [triadic_closure, configuration]
        if model not in models:
            raise ValueError("Invalid model. Expected one of: triadic_closure, configuration")
        if not (type(N)==int and N > 0):
            raise ValueError("Invalid N. Expected an integer greater than zero")
        if (model==triadic_closure and not 0 <= t <= 1) or (model==configuration and not 0 <= t <= 0.2):
            if model==triadic_closure:
                raise ValueError("Invalid t. Expected 0 <= t <= 1")
            elif model==configuration:
                raise ValueError("Invalid t. Expected 0 <= t <= 0.2")

        # if all good, create network
        self.model = model
        self.N = N
        self.t = t
        self.G = self.model(self.t, self.N)

    def graph(self):
        """
        Return the NetworkX graph object for the generated network.

        NetworkX functions can be applied to this version of the graph.
        """
        return self.G

    def plot(self, filename, communities=None):
        """
        Plot the network and save the image for viewing.

        If communities are given, the nodes will be colour-coded according to the community structure.
        Otherwise, all nodes are coloured black.

        Parameters
        ----------
        filename: (required) a filename for the network image
            The image is saved with the given filename, so that the user can open it for viewing.
            The filename is given as a string (e.g. 'network.png').
        communities: (optional) a list of community assignments for each node
            The communities are given in the form returned by the get_communities() function.
            The network nodes are then coloured according to their community assignment.
        """
        fig, ax = plt.subplots()
        pos = nx.spring_layout(self.G)
        nx.draw_networkx_edges(self.G, pos, alpha=0.4)
        
        if communities==None:
            nx.draw_networkx_nodes(self.G, pos, node_size=5, node_color='black')
        else:
            numcom = get_number_communities(communities)
            cm = plt.get_cmap('rainbow')
            colour_list = [cm(1.*i/numcom) for i in range(numcom)]
            colours = [colour_list[communities[i]-1] for i in range(self.N)]
            nx.draw_networkx_nodes(self.G, pos, node_size=5, node_color=colours)

        ax.axis('off')
        ax.figure.savefig(filename, dpi=250)
        print('Network plot saved as {}'.format(filename))

    def get_communities(self, method):
        """
        Runs a given community detection function on the network and returns the community partition.

        A community partition is a list of integers, giving the community number of each network node.

        For example, a return value of the form

        [1, 1, 2, 3, 1, 2, 3, 3 ...]

        indicates that the 1st node is in community 1, the 2nd node is in community 1, and so on.

        Parameter
        ----------
        method: the community detection method to be used (as a string)
            Options are 'modularity', 'infomap', 'spectral' and 'sbm'.
            These correspond to the four methods used in the paper.
        """
        methods = {
            'modularity': get_modularity_communities,
            'infomap': get_infomap_communities,
            'spectral': get_spectral_communities,
            'sbm': get_sbm_communities
        }
        # check input
        if method not in list(methods.keys()):
            raise ValueError("Invalid method. Expected one of: {}".format(list(methods.keys())))

        return(methods[method](self.G))

    def communities_summary(self):
        """
        Runs all four community detection methods on the network and prints a summary of the
        number of communities detected and community similarity between different methods.

        Note: this function may take a while to run for large networks.
        """
        methods = ['modularity', 'infomap', 'spectral', 'sbm']
        combinations = list(itertools.combinations(methods,2))
        communities = dict()
        number_communities = dict()
        similarity = dict()

        # get communities
        for method in methods:
            communities[method] = self.get_communities(method)
            number_communities[method] = get_number_communities(communities[method]) 

        # get similarity between communities
        for combination in combinations:
            method1 = combination[0]
            method2 = combination[1]
            similarity[combination] = get_similarity(communities[method1], communities[method2])

        print('\n########################\n')
        print('Network model: {}'.format(self.model.__name__))
        print('Network size: {}'.format(self.N))
        print('Clustering parameter: {}'.format(self.t))
        print('Clustering coefficient: {:.3f}'.format(nx.transitivity(self.G)))
        print('\n########################\n')

        print('NUMBER OF COMMUNITIES\n')
        for method in methods:
            print('{}: {}'.format(method,number_communities[method]))
        print('\n########################\n')
        
        print('COMMUNITY SIMILARITY (AMI)\n')
        for combination in combinations:
            print('{}: {:.3f}'.format(combination[0] + '-' + combination[1],similarity[combination]))
        print('\n########################\n')
