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
from ..index.embeddings import Index

_Path = Optional[Union[Path, str]]

class GraphAlignment:
    """Given two graphs (after Ricci Flow), compute 
    an approximate alignment between the nodes of both graphs.
    The output is a dictionary with the nodes aligned (as edges)
    and the cost of each one of them. 
    """

    def __init__(self, threshold_alignment: float = 1e5, 
                 dirsave: Optional[_Path] = None,
                 ricci_embedding = True,
                 seq_embedding = False,
                 kmer_size = 1,
                 max_len: Optional[int] = None, 
                 ):
        self.threshold_alignment = threshold_alignment
        
        # node features
        self.ricci_embedding = ricci_embedding
        self.seq_embedding  = seq_embedding
        self.kmer_size = kmer_size
        self.max_len = max_len
        self.compute_node_embeddings = NodeEmbeddings(
                                            ricci_embedding,                                         # weights form ricci flow
                                            seq_embedding, kmer_size=kmer_size, max_len=max_len      # kmer distribution of prefix and suffix (split at node)
                                            )

        if dirsave:
            dirsave=Path(dirsave)
            dirsave.mkdir(exist_ok=True, parents=True)
            self.dirsave = dirsave
        else:
            self.dirsave=None
    
    def __call__(self, ricci_graph1, ricci_graph2, name=None):
        "Alignment of two graphs with Ricci Metric"
        name = "bipartite-graph" if name is None else name

        # alignment
        bipartite_graph = self.create_bipartite_graph(ricci_graph1, ricci_graph2,)
        
        if self.dirsave:
            nx.write_edgelist(bipartite_graph, self.dirsave.joinpath(f"{name}.edgelist"), data=True)

        # Compute alignment between the two graphs            
        alignment = nx.bipartite.minimum_weight_full_matching(bipartite_graph, weight="weight")
        opt_alignment=self.filter_optimal_alignment(alignment, bipartite_graph)
        
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

        # nodes are labeled as '<node>-1', if <node> it belongs to the first graph, and '<node>-2' from the second one
        # the cost/weight of an edge (u,v) correspong for the  
        bipartite_graph = nx.Graph()

        # sp_from_source1, sp_until_sink1 = shortest_paths(ricci_graph1)
        # sp_from_source2, sp_until_sink2 = shortest_paths(ricci_graph2)
        # nodes1_vector_rep = compute_node_embeddings(ricci_graph1, sp_from_source1, sp_until_sink1)
        # nodes2_vector_rep = compute_node_embeddings(ricci_graph2, sp_from_source2, sp_until_sink2)

        nodes1_vector_rep = self.compute_node_embeddings(ricci_graph1, )
        nodes2_vector_rep = self.compute_node_embeddings(ricci_graph2, )

    
        # remove source and sink nodes from both graphs
        ricci_graph1.remove_nodes_from(["source","sink"])
        ricci_graph2.remove_nodes_from(["source","sink"])

        for node1 in ricci_graph1.nodes():
            for node2 in ricci_graph2.nodes(): 

                # cost align two nodes      
                cost_embeddings = self.compute_cost_embeddings(nodes1_vector_rep[node1], nodes2_vector_rep[node2])
                # cost_labels     = self.compute_cost_labels(ricci_graph1.nodes[node1]["label"], ricci_graph2.nodes[node2]["label"])
                # print(f"(cost_embeddings={cost_embeddings}, cost_labels={cost_labels})")
                # weight = 0.5*cost_embeddings + 0.5*cost_labels
                weight = cost_embeddings
                bipartite_graph.add_edge(node1+"-1", node2+"-2", weight=weight) #, cost_embeddings=cost_embeddings, cost_labels=cost_labels)

        logging.info("end - create_bipartite_graph")
        return bipartite_graph
         
    @staticmethod
    def compute_cost_labels(seq1, seq2,):
        # TODO: use EDIT distance
        # https://pypi.org/project/editdistance/
        # https://pypi.org/project/python-Levenshtein/
        "Weight to penalize cost of aligning two nodes"
        # result = parasail.sw_stats(seq1,seq2, open=10, extend=2, matrix=parasail.dnafull)
        # m = result.matches
        # L = max(len(seq1),len(seq2))
        # w = (L-m)/L # weight based on smith-waterman
        w=0
        return w
    
    @staticmethod
    def compute_cost_embeddings(emb1, emb2):
        "Weight to penalize cost of aligning two nodes"
        return np.linalg.norm(emb1 - emb2)