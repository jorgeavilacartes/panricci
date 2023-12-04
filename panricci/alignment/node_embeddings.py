import copy
import logging
import networkx as nx
import numpy as np
from typing import Optional 
from .utils import get_sources_sinks, count_kmers

def shortest_paths(G):
    "Return the shortest paths from source and sink"
    # G_copy = copy.deepcopy(G)
    G_copy=G
    sources, sinks = get_sources_sinks(G_copy)
    G_copy.add_edges_from([("source",node) for node in sources], weight=0, label="N")
    G_copy.add_edges_from([(node,"sink") for node in sinks], weight=0, label="N")
    
    sp_from_source = nx.shortest_path(G_copy, source="source", weight="weight", method="dijkstra")
    sp_until_sink = nx.shortest_path(G_copy, target="sink", weight="weight", method="dijkstra")

    del sp_from_source["source"]
    del sp_until_sink["sink"]
    
    return sp_from_source, sp_until_sink

def feature_from_seq(seq):
    "Return an 1d-array with the distribution of 1-mers in sorted order"
    kmers = count_kmers(seq, k=1)
    return np.array([kmers.get(c,0) for c in "ACGT"]) / len(seq)

def compute_prefix_suffix_feature(G, sp_from_source, sp_until_sink, max_len: Optional[int]=None):
    """compute vector with 1-mer normalized frequencies for prefix (source,node) and suffix (node,sink) sequences
    using the labels of the shortest paths from source and to sink
    """
    prefix_features = dict()
    for start_node, path in sp_from_source.items():
        if start_node != "sink":
            nodes = path[1:] # remove source node
            prefix_label = "".join(G.nodes[node]["label"] for node in nodes)
            
            if max_len:
                prefix_label = prefix_label[::-1][:max_len][::-1]

            prefix_feature =  feature_from_seq(prefix_label)
            prefix_features[start_node] = prefix_feature

    suffix_features = dict()
    for end_node, path in sp_until_sink.items():
        if end_node != "source": 
            nodes = path[:-1] # remove sink node
            suffix_label = "".join(G.nodes[node]["label"] for node in nodes)
            
            if max_len:
                suffix_label = suffix_label[:max_len]

            suffix_feature =  feature_from_seq(suffix_label)
            suffix_features[end_node] = suffix_feature

    nodes_features = {
        node: np.concatenate((prefix_features[node],suffix_features[node])) 
        for node in G.nodes() if node not in ["source", "sink"]
        }
    return nodes_features

def compute_node_embeddings(G, sp_from_source, sp_until_sink):
    """Given a graph after Ricci Flow, compute the vector representation based on the distances 
    from a dummy 'source' and dummy 'sink' nodes. Distances are computed using Dijkstra algorithm. 
    
    Return a dictionary with nodes id as keys, and a 2d-array with its vector representation
    """  

    costs_from_source = dict()
    costs_until_sink = dict()
    for start_node, path in sp_from_source.items():
        nodes = path
        edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1:])]
        cost  = np.sum([G.edges[e]["weight"] for e in edges])
        costs_from_source[start_node] = cost

    for end_node, path in sp_until_sink.items():
        nodes = path
        edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1:])]
        cost  = np.sum([G.edges[e]["weight"] for e in edges])
        costs_until_sink[end_node] = cost

    nodes_features = {
        node: np.array([costs_from_source[node], costs_until_sink[node]]) 
        for node in G.nodes() if node not in ["source", "sink"]
        }
    
    logging.info("end - compute_node_embeddings")
    return nodes_features

class NodeEmbeddings:

    def __init__(self, ricci_embedding=True, seq_embedding=False, max_len: Optional[int]= None):
        self.ricci_embedding = ricci_embedding
        self.seq_embedding = seq_embedding
        # self.kmer_size = kmer_size
        self.max_len = max_len

    def __call__(self, G):
        "Return a dictionary with node embeddings for each node"
        G_copy = copy.deepcopy(G)
        sp_from_source, sp_until_sink = shortest_paths(G_copy)

        if self.ricci_embedding:
            ricci_emb = compute_node_embeddings(G_copy, sp_from_source, sp_until_sink)

        if self.seq_embedding:
            seq_emb = compute_prefix_suffix_feature(G_copy, sp_from_source, sp_until_sink, max_len=self.max_len)
        
        if self.ricci_embedding and self.seq_embedding:
            emb = dict()
            G_copy.remove_nodes_from(["source","sink"])
            # for node in G_copy.nodes():
            #     emb[node] = np.concatenate([ricci_emb[node], seq_emb[node]]) 
            # return emb 
            return {node: np.concatenate([ricci_emb[node], seq_emb[node]]) for node in G_copy.nodes()}
        elif self.ricci_embedding:
            return ricci_emb
        elif self.seq_embedding:
            return seq_emb