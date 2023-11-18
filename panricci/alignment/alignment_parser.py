import pandas as pd

def feature_from_alignment(alignment, graph1, graph2):
    features = []
    for edge, weight in alignment:
        
        edge = sorted(edge, key=lambda x: int(x.split("-")[-1]) , reverse=False)
        node1,node2 = edge # TODO: function to parse solution to node id in the original graphs
        
        node1 = node1.replace("-1","")
        node2 = node2.replace("-2","")

        info1, info2 = graph1.nodes[node1], graph2.nodes[node2]
        features.append(
            dict(
                label1=info1["label"],
                label2=info2["label"],
                node_depth1=info1["node_depth"],
                node_depth2=info1["node_depth"],
            )
        )
    return pd.concat([pd.DataFrame(alignment, columns=["edge","weight"]), pd.DataFrame(features)], axis=1)