from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import logging 
# import parasail         # smith-waterman
import numpy as np
import networkx as nx

# panricci
from .utils import get_sources_sinks
from .node_embeddings import (
    shortest_paths,
    compute_node_embeddings,
    compute_prefix_suffix_feature
)

from .node_embeddings import NodeEmbeddings

_Path = Optional[Union[Path, str]]


class GraphAlignment:
    """Given two graphs (after Ricci Flow), compute 
    an approximate alignment between the nodes of both graphs.
    The output is a dictionary with the nodes aligned (as edges)
    and the cost of each one of them. 
    """

    def __init__(self,  
                ricci_embedding: bool = True,
                seq_embedding: bool = False,
                kmer_size: Optional[int] = None,
                max_len: Optional[int] = None,
                path_save_bipartite: Optional[_Path] = None,
                threshold_alignment: float = 1e5,
                weight_node_labels: float = 0.0,
                log_level: str = "INFO", 
                ):
        
        logging.basicConfig(level=eval(f"logging.{log_level}"),
                    format='[Alignment] %(asctime)s.%(msecs)03d | %(message)s',
                    datefmt='%Y-%m-%d@%H:%M:%S')
        
        self.threshold_alignment = threshold_alignment # maximum cost allowed for the alignment (only used to filter final results)
        
        # node features
        self.ricci_embedding = ricci_embedding
        self.seq_embedding  = seq_embedding
        self.kmer_size = kmer_size
        self.max_len = max_len

        # cost function for alignment of nodes is = weight_node_embeddings*cost_embeddings + weight_node_labels*cost_labels
        self.weight_node_labels = weight_node_labels
        self.weight_node_embeddings = 1.0 - self.weight_node_labels

        # to save bipartite graph used for alignment
        self.path_save_bipartite = path_save_bipartite
        if path_save_bipartite:
            path_save_bipartite=Path(path_save_bipartite)
            path_save_bipartite.parent.mkdir(exist_ok=True, parents=True)
            self.path_save_bipartite = path_save_bipartite
        
    def __call__(self, ricci_graph1, ricci_graph2, name=None):
        "Alignment of two graphs with Ricci Metric"
        name = "bipartite-graph" if name is None else name

        # alignment
        logging.info("Creating bipartite graph")
        bipartite_graph = self.create_bipartite_graph(ricci_graph1, ricci_graph2,)
        
        if self.path_save_bipartite:
            logging.info("Saving bipartite graph")
            nx.write_edgelist(bipartite_graph, self.path_save_bipartite, data=True)

        # Compute alignment between the two graphs
        logging.info("Starting alignment on bipartite graph: minimum-weight-full-matching")            
        alignment = nx.bipartite.minimum_weight_full_matching(bipartite_graph, weight="weight")
        logging.info(f"filtering optimal alignment: avoiding edges with cost > {self.threshold_alignment}")
        opt_alignment=self.filter_optimal_alignment(alignment, bipartite_graph)
        logging.info("Done!")
        return opt_alignment

    def filter_optimal_alignment(self, alignment, bipartite_graph):
        """Get weight and edges in the optimal solution.
        Output is a dictionary (key:value) with edges (keys) sorted by
        the cost of the edge (values).
        """
        logging.info("start - parse_optimal_alignment")
        optimal_alignment=defaultdict(float)
        for node1, node2 in alignment.items():
            edge = (node1,node2)
            edge = tuple(sorted(edge))
            weight = bipartite_graph.edges[edge]["weight"]
            if weight <= self.threshold_alignment:
                optimal_alignment[edge] = weight

        logging.info("end - parse_optimal_alignment")
        return sorted(optimal_alignment.items(), key=lambda d: (d[0],d[1]),reverse=True)
            
    def create_bipartite_graph(self, ricci_graph1, ricci_graph2,):
        """Returns the bipartite graph between
        nodes of both input graphs, to be used for the alignment"""
        logging.info("start - create_bipartite_graph")
        ricci_graph1 = ricci_graph1.copy()
        ricci_graph2 = ricci_graph2.copy()

        # create node embeddings: relative representation of nodes w.r.t. source and sink nodes
        compute_node_embeddings = NodeEmbeddings(
                                    ricci_embedding=self.ricci_embedding,
                                    seq_embedding=self.seq_embedding, 
                                    kmer_size=self.kmer_size, 
                                    max_len=self.max_len      
                                    )

        nodes1_vector_rep = compute_node_embeddings(ricci_graph1,)
        nodes2_vector_rep = compute_node_embeddings(ricci_graph2,)

        # remove source and sink nodes from both graphs
        ricci_graph1.remove_nodes_from(["source","sink"])
        ricci_graph2.remove_nodes_from(["source","sink"])

        # create bipartite graph, create cost of alignment for each pair of nodes (node1,node2)
        bipartite_graph = nx.Graph()

        for node1 in ricci_graph1.nodes():
            for node2 in ricci_graph2.nodes(): 

                # define cost of aligning two nodes in the bipartite graph      
                cost_embeddings = self.compute_cost_embeddings(nodes1_vector_rep[node1], nodes2_vector_rep[node2])

                if self.weight_node_labels > 0:
                    cost_labels = self.compute_cost_labels(ricci_graph1.nodes[node1].get("label"), ricci_graph2.nodes[node2].get("label"))
                else:   
                    cost_labels = np.float64(0) 

                # cost of aligning two nodes is = weight_node_embeddings*cost_embeddings + weight_node_labels*cost_labels
                cost_alignment_nodes = self.weight_node_embeddings * cost_embeddings + self.weight_node_labels * cost_labels
                
                # update bipartite graph with the cost of aligning two nodes
                bipartite_graph.add_edge(node1+"-1", node2+"-2", weight=cost_alignment_nodes, cost_labels=cost_labels, cost_embeddings=cost_embeddings)

        logging.info("end - create_bipartite_graph")
        return bipartite_graph
         
    @staticmethod
    def compute_cost_labels(seq1, seq2,):
        "Weight to penalize cost of aligning two nodes"
        return np.float64(1) if seq1 != seq2 else np.float64(0)
    
    @staticmethod
    def compute_cost_embeddings(emb1, emb2):
        "Weight to penalize cost of aligning two nodes"
        return np.linalg.norm(emb1 - emb2)