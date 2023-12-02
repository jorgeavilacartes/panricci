import networkx as nx

def edgelist_loader(path):
    assert str(path).endswith(".edgelist"), "Must be .edgelist"
    return nx.read_edgelist(
                        path,    # path checkpoint
                        nodetype=str, 
                        create_using=nx.DiGraph
                        )
