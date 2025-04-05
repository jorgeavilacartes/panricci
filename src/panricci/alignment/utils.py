from collections import defaultdict
import pandas as pd

def get_sources_sinks(G):
    "Return sources and sinks nodes in the graph G"
    sources = [x for x in G.nodes() if G.in_degree(x)==0]
    sinks = [x for x in G.nodes() if G.out_degree(x)==0]

    return sources, sinks

def parse_alignment(alignment, graph1, graph2):

    features = []
    for edge, weight in alignment:
        
        edge = sorted(edge, key=lambda x: int(x.split("-")[-1]) , reverse=False)
        node1,node2 = edge
        
        node1 = node1.replace("-1","")
        node2 = node2.replace("-2","")

        info1, info2 = graph1.nodes[node1], graph2.nodes[node2]
        features.append(
            dict(
                edge=edge,
                cost_alignment=weight,
                # cost_embeddings=cost_embeddings, 
                # cost_labels=cost_labels,
                node1=node1,
                node2=node2,
                label1=info1.get("label"),
                label2=info2.get("label"),
                node_depth1=info1.get("node_depth"),
                node_depth2=info2.get("node_depth"),
            )
        )

    return  pd.DataFrame(features)

def count_kmers(sequence: str, k: int): 
    freq_kmer = defaultdict(int)
    # representativity of kmers
    last_j = len(sequence) - k + 1   
    kmers  = (sequence[i:(i+k)] for i in range(last_j))
    # count kmers in a dictionary
    for kmer in kmers:
        if "N" not in kmer:
            freq_kmer[kmer] +=1
    
    return freq_kmer