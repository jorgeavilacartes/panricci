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
                step_size: float = 0.1,  # NEW: Add step size control
                min_weight: float = 1e-6,  # NEW: Minimum weight threshold
                ): 
        logging.basicConfig(level=eval(f"logging.{log_level}"),
                            format='[Ricci-Flow] %(asctime)s.%(msecs)03d | %(message)s',
                            datefmt='%Y-%m-%d@%H:%M:%S')
        
        self.G = G
        self.distribution_nodes = distribution
        self.dirsave = dirsave_graphs
        self.save_last = save_last
        self.save_intermediate_graphs = save_intermediate_graphs
        self.tol = tol_curvature
        self.step_size = step_size  # NEW
        self.min_weight = min_weight  # NEW
 
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
        
        for it in track(range(iterations), total=iterations, description="Ricci-Flow", transient=False):
            
            self._counter_iters += 1
            logging.debug(f"iteration {self._counter_iters}")
            
            new_weights = self.iter_ricci_flow()
            
            nx.set_edge_attributes(self.G, new_weights)
            
            for edge in self.G.edges:
                self.G.edges[edge]["curvature"] = self.compute_curvature(edge)
                       
            if self.dirsave: self.checkpoints(it, name)

            if self.is_curvature_below_tol():
                logging.info(f"Stopping Ricci-Flow in iteration {it}, all curvatures below tol={self.tol}")
                break
            
        return self.G
            
    def iter_ricci_flow(self, eps=None):
        """Compute new weights with Ricci Flow - FIXED VERSION"""
        if eps is None:
            eps = self.step_size  # Use instance step_size
            
        logging.debug(f"Ricci-Flow iteration {self._counter_iters}")

        new_weights = {}
        for edge in self.G.edges():

            current_curvature = self.G.edges[edge]["curvature"]
            current_weight = self.G.edges[edge]["weight"]
        
            # FIX 1: Use step_size to prevent overshooting
            new_weight = current_weight - eps * current_curvature * current_weight
            
            # FIX 2: Clamp to minimum weight to prevent zero/negative weights
            new_weight = max(new_weight, self.min_weight)
            
            # FIX 3: Log warning if weight is clamped
            if new_weight == self.min_weight:
                logging.warning(f"Edge {edge} weight clamped to min_weight={self.min_weight}")
            
            new_weights[edge] = {"weight": new_weight}

        return new_weights
            
    def compute_curvature(self, edge):
        """Compute curvature of an edge"""
        node1, node2 = edge
        d_1 = self.distribution_nodes(node1)
        d_2 = self.distribution_nodes(node2)
        
        W = self.wasserstein(d_1, d_2)
        d = self.G.edges[edge]["weight"]
        
        # FIX 4: Add safety check for division
        if d < self.min_weight:
            logging.warning(f"Edge {edge} has weight {d} < min_weight, using min_weight")
            d = self.min_weight
            
        print("W", W, "d", d)        
        return 1 - W/d

    def wasserstein(self, distribution_node1, distribution_node2):
        """Compute wasserstein distance between two distributions - FIXED VERSION"""

        # 1. extract subgraph
        nodes_subgraph = list(distribution_node1.keys()) + list(distribution_node2.keys())
        subgraph = self.G.subgraph(nodes_subgraph).to_undirected()
        
        # 2. compute all-vs-all shortest paths
        distances_subgraph = nx.shortest_path(subgraph, weight="weight", method="dijkstra")
        
        # 3. create distance matrix
        M = np.zeros((len(distribution_node1), len(distribution_node2)))

        for i, source in enumerate(distribution_node1.keys()):
            for j, target in enumerate(distribution_node2.keys()):
                try:
                    nodes = distances_subgraph[source][target]
                    edges = [(n1, n2) for n1, n2 in zip(nodes[:-1], nodes[1:])]
                    M[i, j] = np.sum([self.G.edges[e]["weight"] for e in edges])
                except:
                    continue
        
        # FIX 5: Safe normalization - avoid division by zero
        max_M = M.max()
        if max_M > 0:
            M = M / max_M  # rescale costs to [0,1]
        else:
            logging.warning("Distance matrix is all zeros, using uniform costs")
            M = np.ones_like(M) * self.min_weight

        # vector with probabilities
        a, b = list(distribution_node1.values()), list(distribution_node2.values())
        
        # FIX 6: Normalize probabilities to ensure they sum to 1
        a = np.array(a)
        b = np.array(b)
        a = a / a.sum() if a.sum() > 0 else np.ones_like(a) / len(a)
        b = b / b.sum() if b.sum() > 0 else np.ones_like(b) / len(b)
        
        # return wasserstein distance
        W = ot.emd2(a, b, M)
        
        # FIX 7: Check for NaN/Inf
        if not np.isfinite(W):
            logging.error(f"Wasserstein distance is not finite: {W}")
            return 0.0
            
        return W
    
    def is_curvature_below_tol(self):
        """Check if all curvatures are below tol"""
        for u, v, data in self.G.edges(data=True):
            if np.abs(data["curvature"]) > self.tol:
                logging.debug(f"curvature of edge K({u},{v}){data['curvature']} is > tol={self.tol}")
                return False
        return True
    
    def checkpoints(self, it, name):
        if self.save_intermediate_graphs:
            name = name if name else "graph-iter"
            path_save = self.dirsave.joinpath(f"{name}-ricciflow-{it+1}.edgelist")
            nx.write_edgelist(self.G, path_save, data=True)

        if self.save_last:
            path_save = self.dirsave.joinpath(f"{name}-ricciflow.edgelist")
            nx.write_edgelist(self.G, path_save, data=True)

    def initialize_edges(self):
        """Initialize edge weights and curvatures"""
        n_edges = len(self.G.edges())
        weight_init = 1 / n_edges

        for edge in self.G.edges():
            # Initialize weight if not present
            self.G.edges[edge]["weight"] = self.G.edges[edge].get("weight", weight_init)
            # Ensure weight is at least min_weight
            self.G.edges[edge]["weight"] = max(self.G.edges[edge]["weight"], self.min_weight)
            
        for edge in self.G.edges():
            self.G.edges[edge]["curvature"] = self.compute_curvature(edge)