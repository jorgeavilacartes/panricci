import numpy as np
import networkx as nx
import ot

from .distributions.variation_graph import DistributionNodes
from typing import Callable
from tqdm import tqdm
from pathlib import Path

class RicciFlow:

    def __init__(self, G, dist_nodes: Callable, dirsave_graphs: str): 
        # TODO: include threshold_curvature: float 
        # to stop ricci flow if all curvatures are at most at threshold_curvature from the average curvature
        # it means that is constant
        
        self.G = G
        self.dist_nodes = dist_nodes
        self.dirsave = Path(dirsave_graphs)
        self.dirsave.mkdir(exist_ok=True, parents=True)
        # self.threshold_curvature = threshold_curvature

        # initialize curvature as 0 if it is not an attribute of the edges
        for edge in G.edges():
            self.G.edges[edge]["curvature"] = self.G.edges[edge].get("curvature",0)

    def __call__(self, iterations: int,):
        # compute curvature for all edges in the graph        
        for it in tqdm(range(iterations), total=iterations, desc="Ricci-Flow"):
            # distances are computed without modify the graph until all edges are used
            distance_by_edge = self.one_iteration()

            # update distances (Ricci-metric) for each edge
            nx.set_edge_attributes(self.G, distance_by_edge)

            # save graph with updated values
            path_save = self.dirsave.joinpath(f"graph-iter-{it}.edgelist")
            # nx.write_weighted_edgelist(self.G, path_save)
            nx.write_edgelist(self.G, path_save, data=True)
            
    def one_iteration(self,):
        "Compute new curvatures, and distances with Ricci Flow"
        distance_by_edge = {}

        for edge in self.G.edges():
            current_distance = self.G.edges[edge]["distance"]
            
            # compute curvature and update it 
            current_curvature = self.curvature(edge)
            self.G.edges[edge]["curvature"] = current_curvature

            # compute distance with Ricci-Flow
            new_distance = current_distance - current_curvature*current_distance
            distance_by_edge[edge] = {"distance": new_distance}

        return distance_by_edge
            
    def curvature(self, edge):
        "Compute curvature of an edge"
        node1, node2 = edge
        dist1 = self.dist_nodes(node1)
        dist2 = self.dist_nodes(node2)
        
        W = self.wasserstein(dist1,dist2)
        d = self.G.edges[edge]["distance"]
        
        return 1 - W/d 

    def wasserstein(self, dist1, dist2):
        "compute wasserstein distance between two distributions"

        # 1. extract subgraph containing only the involved nodes in dist1 and dist2 
        nodes_subgraph = list(dist1.keys()) + list(dist2.keys())
        subgraph = self.G.subgraph(nodes_subgraph)
        
        # 2. compute all-vs-all shortest paths in the subgraph
        distances_subgraph = nx.shortest_path(subgraph, weight="distance", method="dijkstra")#, subgraph.edges()
        
        # 3. create distance matrix for Optimal Transport
        M = np.zeros((len(dist1),len(dist2)))

        for i,source in enumerate(dist1.keys()):
            for j,target in enumerate(dist2.keys()):
                try:
                    nodes = distances_subgraph[source][target]
                    edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1])]
                    M[i,j] = np.sum([self.G.edges[e]["distance"] for e in edges])
                except:
                    continue
        M /= M.max() # rescale costs of matrix to [0,1]

        # vector with probabilities
        a,b=list(dist1.values()), list(dist2.values())
        
        return ot.emd2(a,b,M)