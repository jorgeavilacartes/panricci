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
                node1=node1,
                node2=node2,
                label1=info1["label"],
                label2=info2["label"],
                node_depth1=info1["node_depth"],
                node_depth2=info1["node_depth"],
            )
        )

    return  pd.DataFrame(features)