import numpy as np
import networkx as nx
from pathlib import Path

from .ricci_flow import RicciFlow
from .distributions.variation_graph import DistributionNodes
from .utils.gfa_loader import GFALoader
from .utils.get_source_sink import get_sources_sinks
from .alignment.graph_alignment import GraphAlignment

import logging 
logger = logging.getLogger(name = "PanRicciSimilarity")

class PanRicciSimilarity:
    WEIGHT_FROM_SINK=0
    WEIGHT_TO_SINKS=0
    def __init__(self, alpha: float, iterations=10, directed=True, threshold_alignment=1_000_000, dirsave=None):
        self.alpha = alpha
        self.gfa_loader = GFALoader(directed)
        self.iterations = iterations
        self.threshold_alignment = threshold_alignment
        self.dirsave = dirsave

    def __call__(self, path_gfa1, path_gfa2, name_gfa1, name_gfa2):
        
        # load the graphs from gfa files
        _, _ , graph1 = self.gfa_loader(path_gfa1)
        _, _ , graph2 = self.gfa_loader(path_gfa2)
        
        # Ricci Flow
        sources1, sinks1 = get_sources_sinks(path_gfa1)
        ricci_graph1 = self.apply_ricci_flow(graph1, sources=sources1, sinks=sinks1, name=name_gfa1)        
        
        sources2, sinks2 = get_sources_sinks(path_gfa2)
        ricci_graph2 = self.apply_ricci_flow(graph2, sources=sources2, sinks=sinks2, name=name_gfa2)

        # alignment of the graphs 
        aligner= GraphAlignment()
        alignment = aligner(ricci_graph1,ricci_graph2,path_gfa1, path_gfa2)
        
        metric = self.metric_from_alignment(alignment, ricci_graph1,ricci_graph2)
        return metric , alignment
    
    def apply_ricci_flow(self, graph, name=None, sources=None, sinks=None,):

        # apply ricci flow
        distribution=DistributionNodes(graph, alpha=self.alpha)
        ricci_flow=RicciFlow(graph, distribution, dirsave_graphs=self.dirsave)

        ricci_graph=ricci_flow.run(iterations=self.iterations, save_last=False, save_intermediate_graphs=True, name=name)

        # add dummy source and sink nodes (needed in the alignment step)
        if sources:
            ricci_graph.add_edges_from([("source",node) for node in sources] , weight=self.WEIGHT_FROM_SINK, label="N")
        if sinks:
            ricci_graph.add_edges_from([(node,"sink") for node in sinks] , weight=self.WEIGHT_TO_SINKS, label="N")
        
        return ricci_graph
    
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