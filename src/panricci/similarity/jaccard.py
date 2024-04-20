from typing import Any


class Jaccard:
    
    def __init__(self, threshold_cost: float = 0.5 ) -> None:
        self.threshold_cost = threshold_cost
        
    def __call__(self, alignment, graph1, graph2) -> float:

        # remove source and sink nodes from both graphs
        try:
            graph1.remove_nodes_from(["source","sink"])
            graph2.remove_nodes_from(["source","sink"])
        except:
            pass 
        
        opt_edges = [edge for edge, cost in alignment if cost <= self.threshold_cost]
        metric = 2*len(opt_edges)/(len(graph1)+len(graph2))
        
        return metric