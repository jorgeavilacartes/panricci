from typing import Callable, Optional
from .ricci_flow import RicciFlow

import logging 

class NormalizedRicciFlow(RicciFlow):

    def __init__(self, G, distribution: Callable, sigma: float = 1, dirsave_graphs: str | None = None, overwrite: bool = True, save_last: bool = True, save_intermediate_graphs: bool = False, tol=1e-11):
        self.sigma = sigma
        super().__init__(G, distribution, dirsave_graphs, overwrite, save_last, save_intermediate_graphs, tol)
        
    def iter_ricci_flow(self, eps=1):
        "Compute new curvatures, and distances with Ricci Flow"
        logging.info(f"Normalized Ricci-Flow iteration {self._counter_iters}")

        # TODO: (1) compute curvature with the current weights (2) compute normalization (3) compute weights for the next timestep 
        # 1. reveive a graph with weights and curvatures precomputed
        # 2. compute normalization
        # 3. update weights with normalized-ricci-flow
        
        normalization = self.compute_normalization()

        new_weights = {}
        for edge in self.G.edges():

            # update weights with Ricci-Flow
            curvature = self.G.edges[edge]["curvature"]
            weight = self.G.edges[edge]["weight"]

            new_weight = weight - eps*(curvature - 1/self.sigma * normalization) * weight
            new_weights[edge] = {"weight": new_weight}

        return new_weights

    def compute_normalization(self,):

        normalization = 0
        for edge in self.G.edges:
            normalization += self.G.edges[edge]["curvature"]*self.G.edges[edge]["weight"]

        return normalization
    
    def initialize_edges(self,):
        # initialize curvature as 0 if it is not an attribute of the edges
        n_edges=len(self.G.edges())
        weight_init = self.sigma/n_edges

        for edge in self.G.edges():
            # self.G.edges[edge]["curvature"] = self.G.edges[edge].get("curvature",1)
            self.G.edges[edge]["weight"] = self.G.edges[edge].get("weight",weight_init)

        for edge in self.G.edges():
            self.G.edges[edge]["curvature"] = self.compute_curvature(edge)