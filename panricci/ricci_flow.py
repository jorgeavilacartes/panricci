import numpy as np
import networkx as nx
import ot

from typing import Callable, Optional
from tqdm import tqdm
from pathlib import Path

import logging
class RicciFlow:

    def __init__(self, G, distribution: Callable, dirsave_graphs: Optional[str] = None): 
        # TODO: include threshold_curvature: float 
        # to stop ricci flow if all curvatures are at most at threshold_curvature from the average curvature
        # it means that is constant

        self.G = G
        self.distribution_nodes = distribution # mapping from nodes to its distributions
        self.dirsave = dirsave_graphs
        self.distance_by_edge = None

        self._counter_iters = 0
        if self.dirsave:
            self.dirsave = Path(self.dirsave)
            self.dirsave.mkdir(exist_ok=True, parents=True) # TODO: check is a directory
        # self.threshold_curvature = threshold_curvature

        # initialize curvature as 0 if it is not an attribute of the edges
        n_edges=len(G.edges())
        uniform_weight=1/n_edges
        for edge in G.edges():
            self.G.edges[edge]["curvature"] = self.G.edges[edge].get("curvature",0)
            self.G.edges[edge]["weight"] = self.G.edges[edge].get("weight",uniform_weight)

    def run(self, iterations: int, save_last: bool = True, save_intermediate_graphs: bool=False, name=None):
        # TODO: add callbacks
        # save_last and save_intermediate_graphs should be a callback
        name = name if name else "graph"
        # compute curvature for all edges in the graph        
        for it in tqdm(range(iterations), total=iterations, desc="RicciFlow"):
            self._counter_iters += 1
            logging.info(f"iteration {self._counter_iters}")
            # distances are computed without modify the graph until all edges are used
            new_weights, new_curvatures = self.iter_ricci_flow()
            
            # TODO: EarlyStopping-Callback: check if curvatures changed
            
            # update weights (Ricci-metric) for each edge
            nx.set_edge_attributes(self.G, new_weights)
            nx.set_edge_attributes(self.G, new_curvatures)

            if save_intermediate_graphs:
                name = name if name else "graph-iter"
                # save graph with updated values
                path_save = self.dirsave.joinpath(f"{name}-ricciflow-{it}.edgelist")
                nx.write_edgelist(self.G, path_save, data=True)

            if save_last:
                # save last graph with attributes on its edges
                path_save = self.dirsave.joinpath("{name}-ricciflow.edgelist")
                nx.write_edgelist(self.G, path_save, data=True)

        return self.G
            
    def iter_ricci_flow(self, eps=1):
        "Compute new curvatures, and distances with Ricci Flow"
        logging.info(f"Ricci-Flow iteration {self._counter_iters}")

        new_weights = {}
        new_curvatures = {}
        for edge in self.G.edges():

            # compute curvature and update it
            current_curvature = self.compute_curvature(edge)

            # update weights with Ricci-Flow
            current_weight = self.G.edges[edge]["weight"]
        
            new_weight = current_weight - eps*current_curvature*current_weight
            
            new_weights[edge] = {"weight": new_weight}
            new_curvatures[edge] = {"curvature": current_curvature}

        return new_weights, new_curvatures
            
    def compute_curvature(self, edge):
        "Compute curvature of an edge"
        node1, node2 = edge
        dist1 = self.distribution_nodes(node1)
        dist2 = self.distribution_nodes(node2)
        
        W = self.wasserstein(dist1,dist2)
        d = self.G.edges[edge]["weight"] #FIXME: should be distance computed with dijkstra
        
        return 1 - W/d 

    def wasserstein(self, dist1, dist2):
        "compute wasserstein distance between two distributions"

        # 1. extract subgraph containing only the involved nodes in dist1 and dist2 
        nodes_subgraph = list(dist1.keys()) + list(dist2.keys())
        subgraph = self.G.subgraph(nodes_subgraph)
        
        # 2. compute all-vs-all shortest paths in the subgraph
        distances_subgraph = nx.shortest_path(subgraph, weight="weight", method="dijkstra")
        
        # 3. create distance matrix for Optimal Transport
        M = np.zeros((len(dist1),len(dist2)))

        for i,source in enumerate(dist1.keys()):
            for j,target in enumerate(dist2.keys()):
                try:
                    nodes = distances_subgraph[source][target]
                    edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1:])]
                    M[i,j] = np.sum([self.G.edges[e]["weight"] for e in edges])
                except:
                    continue
        M /= M.max() # rescale costs of matrix to [0,1]

        # vector with probabilities
        a,b=list(dist1.values()), list(dist2.values())
        
        return ot.emd2(a,b,M)