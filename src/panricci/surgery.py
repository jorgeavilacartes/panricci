
class Surgery:

    def __call__(self, G, cutting_threshold: float = 1e5):
        """Given an input graph with weights, remove edges with weight larger than
        the cutting_threshold"""
        
        edges_to_remove = [ e for e in G.edges() if e["weight"] > cutting_threshold ]
        G.remove_edges_from(edges_to_remove)

        return G