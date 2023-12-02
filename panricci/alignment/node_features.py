import logging
import networkx as nx
import numpy as np

from .utils import get_sources_sinks, count_kmers

def shortest_paths(G):
    "Return the shortest paths from source and sink"
    sources, sinks = get_sources_sinks(G)
    G.add_edges_from([("source",node) for node in sources], weight=0, label="N")
    G.add_edges_from([(node,"sink") for node in sinks], weight=0, label="N")
    
    sp_from_source = nx.shortest_path(G, source="source", weight="weight", method="dijkstra")
    sp_until_sink = nx.shortest_path(G, target="sink", weight="weight", method="dijkstra")

    del sp_from_source["source"]
    del sp_until_sink["sink"]
    
    return sp_from_source, sp_until_sink

def feature_from_seq(seq):
    "Return a list with the distribution of 1-mers in sorted order"
    kmers = count_kmers(seq, k=1)
    return np.array([kmers.get(c,0) for c in "ACGT"]) / len(seq)

def compute_prefix_suffix_feature(G, sp_from_source, sp_until_sink):
    prefix_features = dict()
    suffix_features = dict()
    for start_node, path in sp_from_source.items():
        nodes = path
        prefix_label = "".join(G.nodes[node]["label"] for node in nodes)
        prefix_feature =  feature_from_seq(prefix_label)
        prefix_features[start_node] = prefix_feature

    for end_node, path in sp_until_sink.items():
        nodes = path
        suffix_label = "".join(G.nodes[node]["label"] for node in nodes)
        suffix_feature =  feature_from_seq(suffix_label)
        suffix_features[end_node] = suffix_feature

    nodes_features = {
        node: np.array(prefix_features[node]+suffix_features[node]) 
        for node in G.nodes() if node not in ["source", "sink"]
        }
    
    logging.info("end - compute_node_embeddings")
    return nodes_features

def compute_node_embeddings(G,):# sp_from_source, sp_until_sink):
    """Given a graph after Ricci Flow, compute the vector representation based on the distances 
    from a dummy 'source' and dummy 'sink' nodes. Distances are computed using Dijkstra algorithm. 
    
    Return a dictionary with nodes id as keys, and a 2d-array with its vector representation
    """  
    
    logging.info("start - compute_node_embeddings")
    sources, sinks = get_sources_sinks(G)
    G.add_edges_from([("source",node) for node in sources], weight=0, label="N")
    G.add_edges_from([(node,"sink") for node in sinks], weight=0, label="N")
    
    sp_from_source = nx.shortest_path(G, source="source", weight="weight", method="dijkstra")
    sp_until_sink = nx.shortest_path(G, target="sink", weight="weight", method="dijkstra")

    del sp_from_source["source"]
    del sp_until_sink["sink"]
    
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
