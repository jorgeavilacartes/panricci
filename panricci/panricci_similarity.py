import numpy as np
import networkx as nx
from .ricci_flow import RicciFlow
from .distributions.variation_graph import DistributionNodes
from .utils.gfa_loader import GFALoader
from .utils.get_source_sink import get_sources_sinks
from .alignment.graph_alignment import GraphAlignment

class PanRicciSimilarity:
    DISTANCE_SOURCE=0
    DISTANCE_SINKS=0
    def __init__(self, alpha: float, iterations=10,undirected=False, weight_source_sink=0.1, threshold_alignment=0.5):
        self.alpha = alpha
        self.gfa_loader = GFALoader(undirected)
        self.iterations = iterations
        self.weight_source_sink = weight_source_sink
        self.threshold_alignment = threshold_alignment


    def __call__(self, path_gfa1, path_gfa2,):
        
        # load the graphs from gfa files
        _, _ , graph1 = self.gfa_loader(path_gfa1)
        _, _ , graph2 = self.gfa_loader(path_gfa2)
        
        # Ricci Flow
        ricci_graphs=[]
        for graph, path_graph in zip([graph1, graph2],[path_gfa1, path_gfa2]):
            
            # apply ricci flow
            distribution=DistributionNodes(graph, alpha=self.alpha)
            ricci_flow=RicciFlow(graph, distribution)
            ricci_graph=ricci_flow(iterations=self.iterations, save_last=False)
            
            # add dummy source and sink nodes (needed in the alignment step)
            sources, sinks = get_sources_sinks(path_graph)
            ricci_graph.add_edges_from([("source",node) for node in sources] , distance=self.DISTANCE_SOURCE, label="N")
            ricci_graph.add_edges_from([(node,"sink") for node in sinks] , distance=self.DISTANCE_SINKS, label="N")
            
            ricci_graphs.append(ricci_graph)

        # alignment of the graphs 
        aligner= GraphAlignment()
        alignment = aligner(ricci_graphs[0],ricci_graphs[1],path_gfa1, path_gfa2)
        
        metric = self.metric_from_alignment(alignment, ricci_graphs[0],ricci_graphs[1])
        return metric , alignment
    
    def metric_from_alignment(self, alignment, graph1, graph2):

        # remove source and sink nodes from both graphs
        try:
            graph1.remove_nodes_from(["source","sink"])
            graph2.remove_nodes_from(["source","sink"])
        except:
            pass 
        # _costs = [ 2-2*elem[1] for elem in alignment]
        opt_edges = [edge for edge, cost in alignment if cost <= self.threshold_alignment]
        metric = 2*len(alignment)/(len(graph1)+len(graph2))
        # metric = 2*sum([alignment[edge] for edge in opt_edges])/(len(graph1)+len(graph2))
        
        return metric