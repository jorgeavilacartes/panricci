from typing import Optional
from collections import defaultdict
import networkx as nx

NUC_COMPLEMENT = {n:c for n,c in zip ("ACGT","TGCA")}
def reverse_complement(seq):
    seq = seq[::-1]
    return "".join(NUC_COMPLEMENT[n] for n in seq)

class GFALoader:
    "Load GFA as a directed network with networkx"
    
    def __init__(self, undirected=False):
        self.undirected = undirected

    def __call__(self, path_gfa: str, path_chechkpoint_ricciflow: Optional[str]):
        
        G = nx.Graph() if self.undirected else nx.DiGraph()

        nodes, edges = self.load_nodes_edges(path_gfa)
        node_depth = self.node_depth(path_gfa, nodes)
        
        # add node_depth as attribute to nodes
        for node, node_depth in node_depth.items():
            nodes[node]["node_depth"] = node_depth
        
        G.add_nodes_from([(node, attrs) for node, attrs in nodes.items()])
        G.add_edges_from(edges,)

        if path_chechkpoint_ricciflow:
            data_edges = self.load_edge_attributes(path_chechkpoint_ricciflow)
            nx.set_edge_attributes(G, data_edges)
        
        return G

    def load_nodes_edges(self, path_gfa: str):
        nodes = dict()
        edges = []
        with open(path_gfa, "r") as fp:
            for line in fp.readlines():    
                # nodes
                if line.startswith("S"):
                    try:
                        _, nodeid, label = line.replace("\n","").split("\t")
                    except:
                        _, nodeid, label, _ = line.replace("\n","").split("\t")
                    nodes[nodeid] = {"label": label, "len": len(label)}

                # edges: L	4	+	86	+	0M
                if line.startswith("L"):
                    _, nodeid1, _, nodeid2, _, _ = line.replace("\n","").split("\t") 
                    edges.append((nodeid1, nodeid2))
        
        return nodes, edges
    
    def node_depth(self, path_gfa, nodes):
        paths_by_node=defaultdict(list)
        # once nodes are loaded, check the paths
        n_paths = 0
        with open(path_gfa, "r") as fp:
            for line in fp.readlines():

                # paths
                nodes_path=[]
                if line.startswith("P"):
                    n_paths += 1 
                    _, seq_id, path, *_ = line.replace("\n","").split("\t")

                    labels_path=[]
                    for node in path.split(","):
                        if "-" in node: # reverse label '-' next to node
                            nodeid = node.replace("-","")
                            label  = reverse_complement(nodes[nodeid]["label"])
                            labels_path.append(label.upper())
                        else: # forward label '+' next to node
                            nodeid = node.replace("+","")
                            label  = nodes[nodeid]["label"]
                            labels_path.append(label.upper())
                        
                        paths_by_node[nodeid].append(seq_id)
    
        return {node: len(set(paths))/n_paths for node, paths in paths_by_node.items()}

    def load_edge_attributes(self, path_chechkpoint_ricciflow):
        "Load RicciFlow checkpoint from .edgelist"
        G = self.edgelist_loader(path_chechkpoint_ricciflow)
        data_edges = dict()
        for u,v,data in G.edges(data=True):
            data_edges[(u,v)] = data
        return data_edges

    def edgelist_loader(self, path):
        assert str(path).endswith(".edgelist"), "Must be .edgelist"
        return nx.read_edgelist(
                            path,    # path checkpoint
                            nodetype=str, 
                            create_using=nx.DiGraph
                            )

