import numpy as np
from collections import defaultdict
import parasail
import networkx as nx
from panricci.utils.get_source_sink import get_sources_sinks

class GraphAlignment:
    """Given two graphs (after Ricci Flow), compute 
    an approximate alignment between the nodes of both graphs.
    The output is a dictionary with the nodes aligned (as edges)
    and the cost of each one of them. 
    """

    def __init__(self, threshold_alignment: float = 0.5):
        self.threshold_alignment = threshold_alignment
    
    def __call__(self, ricci_graph1,ricci_graph2,path_gfa1, path_gfa2):
        "Alignment of two graphs with Ricci Metric"
        # alignment
        bipartite_graph = self.create_bipartite_graph(ricci_graph1, ricci_graph2, path_gfa1, path_gfa2)
        alignment = nx.bipartite.minimum_weight_full_matching(bipartite_graph, weight="weight")

        opt_alignment=self.parse_optimal_alignment(alignment, bipartite_graph)
        
        return opt_alignment

    def parse_optimal_alignment(self, alignment, bipartite_graph):
        """Get weight and edges in the optimal solution.
        Output is a dictionary (key:value) with edges (keys) sorted by
        the cost of the edge (values).
        """
        optimal_alignment=defaultdict(float)
        for node1, node2 in alignment.items():
            edge = (node1,node2)
            edge = tuple(sorted(edge))
            weight = bipartite_graph.edges[edge]["weight"]
            if weight <= self.threshold_alignment:
                optimal_alignment[edge] = weight
    
        return sorted(optimal_alignment.items(), key=lambda d: (d[0],d[1]),reverse=True)
            
    def create_bipartite_graph(self, ricci_graph1, ricci_graph2, path_gfa1, path_gfa2):
        """Returns the bipartite graph between
        nodes of both input graphs, to be used for the alignment"""
        ricci_graph1 = ricci_graph1.copy()
        ricci_graph2 = ricci_graph2.copy()
        
        # remove source and sink nodes from both graphs
        ricci_graph1.remove_nodes_from(["source","sink"])
        ricci_graph2.remove_nodes_from(["source","sink"])

        # nodes are labeled as '<node>-1', if <node> it belongs to the first graph, and '<node>-2' from the second one
        # the cost/weight of an edge (u,v) correspong for the  
        bipartite_graph = nx.Graph()

        ricci_graph1_nodes= [n for n in ricci_graph1.nodes()]
        ricci_graph2_nodes= [n for n in ricci_graph2.nodes()]
        for node1 in ricci_graph1_nodes:
            for node2 in ricci_graph2_nodes:
                graph1_vector_rep = self.compute_vector_representation(ricci_graph1, path_gfa1)
                graph2_vector_rep = self.compute_vector_representation(ricci_graph2, path_gfa2)
                
                weight_ricci_metric = np.linalg.norm(graph1_vector_rep[node1] - graph2_vector_rep[node2])
                weight_alignment    = self.compute_weight_alignment(ricci_graph1.nodes[node1]["label"], ricci_graph2.nodes[node2]["label"])
                
                weight = 0.5*weight_ricci_metric + 0.5*weight_alignment
                bipartite_graph.add_edge(node1+"-1", node2+"-2", weight=weight)

        return bipartite_graph
    
    @staticmethod
    def compute_vector_representation(ricci_graph, path_gfa):
        """Given a graph after Ricci Flow and with source and sink nodes,
        compute the vector representation of each node w.r.t. source and 
        sink nodes using the ricci metric.
        """  
        sources, sinks = get_sources_sinks(path_gfa)
        ricci_graph.add_edges_from([("source",node) for node in sources] , distance=1, label="N")
        ricci_graph.add_edges_from([(node,"sink") for node in sinks] , distance=1, label="N")
        
        sp_from_source = nx.shortest_path(ricci_graph, source="source", weight="distance", method="dijkstra")
        sp_until_sink = nx.shortest_path(ricci_graph, target="sink", weight="distance", method="dijkstra")

        del sp_from_source["source"]
        del sp_until_sink["sink"]
        
        costs_from_source = dict()
        costs_until_sink = dict()
        for start_node, path in sp_from_source.items():
            nodes = path
            edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1:])]
            cost  = np.sum([ricci_graph.edges[e]["distance"] for e in edges])
            costs_from_source[start_node] = cost

        for end_node, path in sp_until_sink.items():
            nodes = path
            edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1:])]
            cost  = np.sum([ricci_graph.edges[e]["distance"] for e in edges])
            costs_until_sink[end_node] = cost


        nodes_vector_representation = {node: np.array([costs_from_source[node], costs_until_sink[node]]) for node in ricci_graph.nodes() if node not in ["source", "sink"]}
        return nodes_vector_representation
    
    @staticmethod
    def compute_weight_alignment(seq1, seq2,):
        "Weight to penalize cost of aligning two nodes"
        result = parasail.sw_stats(seq1,seq2, open=10, extend=2, matrix=parasail.dnafull)
        m = result.matches
        L = max(len(seq1),len(seq2))
        w = (L-m)/L # weight based on smith-waterman
        return w