import numpy as np
import networkx as nx
import ot

from typing import Callable, Optional
from rich.progress import track
from pathlib import Path

import logging

class RicciFlow:

    def __init__(self, 
                G, 
                distribution: Callable, 
                dirsave_graphs: Optional[str] = None, 
                overwrite: bool = True, 
                save_last: bool = True, 
                save_intermediate_graphs: bool=False,
                tol_curvature = 1e-11,
                log_level: str = "INFO",
                ): 
        # TODO: include threshold_curvature: float 
        # to stop ricci flow if all curvatures are at most at threshold_curvature from the average curvature
        # it means that is constant

        logging.basicConfig(level=eval(f"logging.{log_level}"),
                            format='[Ricci-Flow] %(asctime)s.%(msecs)03d | %(message)s',
                            datefmt='%Y-%m-%d@%H:%M:%S')
        
        self.G = G
        self.distribution_nodes = distribution # mapping from nodes to its distributions
        self.dirsave = dirsave_graphs
        self.dirsave_graphs = dirsave_graphs
        self.save_last = save_last
        self.save_intermediate_graphs = save_intermediate_graphs
        self.tol = tol_curvature # tolerance of minimum curvature to stop Ricci-Flow when |K(u,v)| < tol for all (u,v) edge of G
 
        self._counter_iters = 0

        if overwrite is False and self.dirsave:
            if self.dirsave.exists():
                raise(f"{str(self.dirsave)} already exists. To overwrite results set 'overwrite=True' or use another directory to save results.")

        if self.dirsave:
            self.dirsave = Path(self.dirsave)
            self.dirsave.mkdir(exist_ok=True, parents=True)

        self.initialize_edges()


    def run(self, iterations: int, name=None):

        name = name if name else "graph"
        
        # compute curvature for all edges in the graph       
        # while not self.is_curvature_below_tol() or self._counter_iters > iterations: 
        for it in track(range(iterations), total=iterations, description="Ricci-Flow", transient=False):
            
            self._counter_iters += 1
            logging.debug(f"iteration {self._counter_iters}")
            
            # distances are computed without modify the graph until all edges are used
            new_weights = self.iter_ricci_flow()
            
            # update weights (Ricci-metric) for each edge
            nx.set_edge_attributes(self.G, new_weights)
            
            # update curvature with the current weight 
            for edge in self.G.edges:
                self.G.edges[edge]["curvature"] = self.compute_curvature(edge)
                       
            self.checkpoints(it, name)

            if self.is_curvature_below_tol():
                logging.info(f"Stopping Ricci-Flow in iteration {it}, all curvatures below tol={self.tol}")
                break
            
        return self.G
            
    def iter_ricci_flow(self, eps=1):
        "Compute new curvatures, and distances with Ricci Flow"
        logging.debug(f"Ricci-Flow iteration {self._counter_iters}")

        new_weights = {}
        new_curvatures = {}
        for edge in self.G.edges():

            # compute curvature and update it
            current_curvature = self.G.edges[edge]["curvature"] #self.compute_curvature(edge)

            # update weights with Ricci-Flow
            current_weight = self.G.edges[edge]["weight"]
        
            new_weight = current_weight - eps*current_curvature*current_weight
            
            new_weights[edge] = {"weight": new_weight}
            # new_curvatures[edge] = {"curvature": current_curvature}

        return new_weights
            
    def compute_curvature(self, edge):
        "Compute curvature of an edge"
        node1, node2 = edge
        d_1 = self.distribution_nodes(node1)
        d_2 = self.distribution_nodes(node2)
        
        W = self.wasserstein(d_1,d_2)
        d = self.G.edges[edge]["weight"]
        
        return 1 - W/d

    def wasserstein(self, distribution_node1, distribution_node2):
        "compute wasserstein distance between two distributions of nodes"

        # 1. extract subgraph containing only the involved nodes in dist1 and dist2 
        nodes_subgraph = list(distribution_node1.keys()) + list(distribution_node2.keys())
        subgraph = self.G.subgraph(nodes_subgraph)
        
        # 2. compute all-vs-all shortest paths in the subgraph
        distances_subgraph = nx.shortest_path(subgraph, weight="weight", method="dijkstra")
        
        # 3. create distance matrix for Optimal Transport
        M = np.zeros((len(distribution_node1),len(distribution_node2)))

        for i,source in enumerate(distribution_node1.keys()):
            for j,target in enumerate(distribution_node2.keys()):
                try:
                    nodes = distances_subgraph[source][target]
                    edges = [(n1,n2) for n1,n2 in zip(nodes[:-1], nodes[1:])]
                    M[i,j] = np.sum([self.G.edges[e]["weight"] for e in edges])
                except:
                    continue
        M /= M.max() # rescale costs of matrix to [0,1]

        # vector with probabilities
        a,b=list(distribution_node1.values()), list(distribution_node2.values())
        
        # return wassersetein distance between the two distributions
        return ot.emd2(a,b,M)
    
    def is_curvature_below_tol(self,):
        "Check if all curvatures are below tol"
        for u,v, data in self.G.edges(data=True):
            if np.abs(data["curvature"]) > self.tol:
                logging.debug(f"curvature of edge K({u},{v}){data['curvature']} is > tol={self.tol}")
                return False
        return True
    
    def checkpoints(self, it, name):
        if self.save_intermediate_graphs:
            name = name if name else "graph-iter"
            # save graph with updated values
            path_save = self.dirsave.joinpath(f"{name}-ricciflow-{it+1}.edgelist")
            nx.write_edgelist(self.G, path_save, data=True)

        if self.save_last:
            # save last graph with attributes on its edges
            path_save = self.dirsave.joinpath("{name}-ricciflow.edgelist")
            nx.write_edgelist(self.G, path_save, data=True)

    def initialize_edges(self,):
        # initialize curvature as 0 if it is not an attribute of the edges
        n_edges=len(self.G.edges())
        weight_init = 1/n_edges

        for edge in self.G.edges():
            self.G.edges[edge]["weight"] = self.G.edges[edge].get("weight",weight_init)

        for edge in self.G.edges():
            self.G.edges[edge]["curvature"] = self.compute_curvature(edge)