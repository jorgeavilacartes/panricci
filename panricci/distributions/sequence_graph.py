import networkx as nx

class DistributionNodes:
    """Distribution of nodes is weighted by
    the number of neighbors."""
    def __init__(self, G, alpha=0.5):
        self.G = G
        self.alpha = alpha
        self.distribution_nodes = self.compute_distribution_all_nodes()

    def __call__(self, node):
        "Return a dictionary with for each neighbor, including input 'node'"
        return self.distribution_nodes[node]

    def compute_distribution_all_nodes(self,):
        distribution_all_nodes = dict()
        for node in self.G.nodes():
            distribution_all_nodes[node] = self._compute_distribution(node)

        return distribution_all_nodes
    
    def _compute_distribution(self, node):
        """Returns a dictionary with the distribution of a node. 
        Keys are the neighbors of 'node', including the query 'node'.
        The values are the probability distribution defined by alpha and
        the paths"""
        neighborhood = list(self.G.neighbors(node))
        n_neighbors = len(neighborhood)
        # compute distribution for each neighbor of 'node'
        distribution = {node: self.alpha} if len(neighborhood)>0 else {node: 1.}
        for neighbor in neighborhood:
            distribution[neighbor] = (1-self.alpha) / n_neighbors
        
        return distribution