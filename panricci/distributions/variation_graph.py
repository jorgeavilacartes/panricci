import networkx as nx

class DistributionNodes:
    """Distribution of nodes is weighted by
    the node depth: number of paths using that node.
    node depth is in [0,1],. where 0 means no sequence
    is using that node, and 1 means all the sequences
    uses that node."""
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
        total_depth  = sum([self.G.nodes[n]["node_depth"] for n in neighborhood])

        # compute distribution for each node in the neighborhood of 'node', and for 'node'
        distribution = {node: self.alpha} if len(neighborhood)>0 else {node: 1.}
        for neighbor in neighborhood:
            neighbor_depth = self.G.nodes[neighbor]["node_depth"]
            distribution[neighbor] = (1-self.alpha) * neighbor_depth / total_depth
        
        return distribution