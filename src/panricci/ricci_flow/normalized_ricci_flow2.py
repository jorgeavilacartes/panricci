import numpy as np
import networkx as nx
from typing import Callable, Optional
import logging

# Assumes RicciFlow is imported from the fixed version
from .ricci_flow2 import RicciFlow

class NormalizedRicciFlow(RicciFlow):
    """
    Normalized Ricci Flow that preserves the total sum of edge weights.
    
    The normalization term ensures that sum(weights) = sigma remains constant,
    preventing the graph from collapsing to zero weights.
    
    Update rule: d_new = d - eps * (K - normalization/sigma) * d
    where normalization = sum(K_e * d_e) over all edges
    """
    
    def __init__(self, 
                 G, 
                 distribution: Callable, 
                 sigma: float = 1.0, 
                 dirsave_graphs: Optional[str] = None, 
                 overwrite: bool = True, 
                 save_last: bool = True, 
                 save_intermediate_graphs: bool = False, 
                 tol_curvature: float = 1e-11,
                 log_level: str = "INFO",
                 step_size: float = 0.1,
                 min_weight: float = 1e-6,
                 ):
        """
        Parameters:
        -----------
        sigma : float
            Target sum of all edge weights (default: 1.0)
            The normalization ensures sum(weights) = sigma is preserved
        """
        self.sigma = sigma
        
        # Call parent constructor with all parameters
        super().__init__(
            G, 
            distribution, 
            dirsave_graphs=dirsave_graphs, 
            overwrite=overwrite, 
            save_last=save_last, 
            save_intermediate_graphs=save_intermediate_graphs, 
            tol_curvature=tol_curvature,
            log_level=log_level,
            step_size=step_size,
            min_weight=min_weight
        )
    
    def iter_ricci_flow(self, eps=None):
        """
        Compute new weights with Normalized Ricci Flow.
        
        The normalization term (1/sigma * sum(K_e * d_e)) ensures that
        the total sum of weights remains constant at sigma.
        """
        if eps is None:
            eps = self.step_size
            
        logging.debug(f"Normalized Ricci-Flow iteration {self._counter_iters}")
        
        # Compute normalization term: sum(K_e * d_e) / sigma
        normalization = self.compute_normalization()
        
        logging.debug(f"Normalization term: {normalization / self.sigma:.6f}")
        logging.debug(f"Current sum of weights: {self.compute_total_weight():.6f}")
        
        new_weights = {}
        for edge in self.G.edges():
            curvature = self.G.edges[edge]["curvature"]
            weight = self.G.edges[edge]["weight"]
            
            # Normalized Ricci Flow update rule
            normalized_curvature = curvature - (normalization / self.sigma)
            new_weight = weight - eps * normalized_curvature * weight
            
            # Clamp to minimum weight to prevent numerical issues
            new_weight = max(new_weight, self.min_weight)
            
            if new_weight == self.min_weight:
                logging.warning(f"Edge {edge} weight clamped to min_weight={self.min_weight}")
            
            new_weights[edge] = {"weight": new_weight}
        
        # Log the sum of new weights for verification
        new_total = sum(w["weight"] for w in new_weights.values())
        logging.debug(f"New sum of weights: {new_total:.6f} (target: {self.sigma:.6f})")
        
        return new_weights
    
    def compute_normalization(self):
        """
        Compute the normalization term: sum(K_e * d_e) over all edges.
        
        This term is subtracted from each curvature to ensure the sum of
        weights remains constant.
        
        Returns:
        --------
        float : sum of (curvature * weight) for all edges
        """
        normalization = 0.0
        for edge in self.G.edges():
            curvature = self.G.edges[edge]["curvature"]
            weight = self.G.edges[edge]["weight"]
            normalization += curvature * weight
        
        return normalization
    
    def compute_total_weight(self):
        """
        Compute the current sum of all edge weights.
        
        This should remain approximately equal to sigma throughout iterations.
        """
        return sum(self.G.edges[edge]["weight"] for edge in self.G.edges())
    
    def initialize_edges(self):
        """
        Initialize edge weights such that their sum equals sigma.
        
        Each edge gets weight = sigma / n_edges
        """
        n_edges = len(self.G.edges())
        
        if n_edges == 0:
            logging.warning("Graph has no edges!")
            return
        
        # Initialize weights to ensure sum = sigma
        weight_init = self.sigma / n_edges
        
        for edge in self.G.edges():
            self.G.edges[edge]["weight"] = self.G.edges[edge].get("weight", weight_init)
            # Ensure weight is at least min_weight
            self.G.edges[edge]["weight"] = max(self.G.edges[edge]["weight"], self.min_weight)
        
        # Normalize to exactly sigma if clamping changed the sum
        total_weight = self.compute_total_weight()
        if abs(total_weight - self.sigma) > 1e-10:
            scale_factor = self.sigma / total_weight
            for edge in self.G.edges():
                self.G.edges[edge]["weight"] *= scale_factor
        
        # Compute initial curvatures
        for edge in self.G.edges():
            self.G.edges[edge]["curvature"] = self.compute_curvature(edge)
        
        logging.info(f"Initialized {n_edges} edges with total weight = {self.compute_total_weight():.6f}")
    
    def run(self, iterations: int, name=None):
        """
        Run Normalized Ricci Flow with verification of weight conservation.
        """
        logging.info(f"Starting Normalized Ricci Flow with sigma={self.sigma}")
        logging.info(f"Initial total weight: {self.compute_total_weight():.6f}")
        
        result = super().run(iterations, name)
        
        final_weight = self.compute_total_weight()
        logging.info(f"Final total weight: {final_weight:.6f}")
        logging.info(f"Weight preservation error: {abs(final_weight - self.sigma):.2e}")
        
        return result
