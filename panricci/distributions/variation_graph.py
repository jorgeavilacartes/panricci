import networkx as nx

class DistributionNodes:

    def __init__(self, G, alpha=0.5):
        self.G = G
        self.alpha = alpha

    def __call__(self, node):
        "Return a dictionary with for each neighbor, including input 'node'"
        
        neighborhood = list(self.G.neighbors(node))
        total_depth  = sum([self.G.nodes[n]["node_depth"] for n in neighborhood])

        # compute distribution for each node in the neighborhood of 'node', and for 'node'
        distribution = {node: self.alpha} if len(neighborhood)>0 else {node: 1.}
        for neighbor in neighborhood:
            neighbor_depth = self.G.nodes[neighbor]["node_depth"]
            distribution[neighbor] = (1-self.alpha) * neighbor_depth / total_depth
        
        return distribution
