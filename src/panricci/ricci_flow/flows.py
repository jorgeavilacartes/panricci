import numpy as np
import networkx as nx
import ot
from abc import ABC, abstractmethod
from typing import Callable, Optional, Dict, Tuple, Any, Literal
from rich.progress import track
from pathlib import Path
import logging


class LearningRateScheduler:
    """
    Learning rate scheduler for Ricci Flow, similar to optimizers in deep learning.
    """
    
    def __init__(
        self, 
        initial_lr: float = 0.1,
        schedule: Literal["constant", "step", "exponential", "adaptive"] = "constant",
        step_size: int = 10,
        gamma: float = 0.5,
        patience: int = 5,
        factor: float = 0.5,
        min_lr: float = 1e-6,
    ):
        """
        Initialize learning rate scheduler.
        
        Parameters
        ----------
        initial_lr : float
            Initial learning rate (epsilon)
        schedule : str
            Type of schedule: "constant", "step", "exponential", or "adaptive"
        step_size : int
            For "step": decrease LR every step_size iterations
        gamma : float
            Multiplicative factor for LR decay
        patience : int
            For "adaptive": number of iterations to wait before reducing LR
        factor : float
            For "adaptive": factor to reduce LR by
        min_lr : float
            Minimum learning rate
        """
        self.initial_lr = initial_lr
        self.current_lr = initial_lr
        self.schedule = schedule
        self.step_size = step_size
        self.gamma = gamma
        self.patience = patience
        self.factor = factor
        self.min_lr = min_lr
        
        # For adaptive scheduling
        self.best_metric = float('inf')
        self.iterations_since_improvement = 0
        
    def step(self, iteration: int, metric: Optional[float] = None) -> float:
        """
        Update learning rate based on schedule.
        
        Parameters
        ----------
        iteration : int
            Current iteration number
        metric : Optional[float]
            Metric to monitor (e.g., max curvature) for adaptive scheduling
            
        Returns
        -------
        float
            Updated learning rate
        """
        if self.schedule == "constant":
            return self.current_lr
            
        elif self.schedule == "step":
            if iteration > 0 and iteration % self.step_size == 0:
                self.current_lr = max(self.current_lr * self.gamma, self.min_lr)
                logging.info(f"LR decreased to {self.current_lr:.6e}")
            return self.current_lr
            
        elif self.schedule == "exponential":
            self.current_lr = max(
                self.initial_lr * (self.gamma ** iteration), 
                self.min_lr
            )
            return self.current_lr
            
        elif self.schedule == "adaptive":
            if metric is None:
                return self.current_lr
                
            # Check if metric improved
            if metric < self.best_metric * 0.999:  # 0.1% improvement threshold
                self.best_metric = metric
                self.iterations_since_improvement = 0
            else:
                self.iterations_since_improvement += 1
                
            # Reduce LR if no improvement for patience iterations
            if self.iterations_since_improvement >= self.patience:
                old_lr = self.current_lr
                self.current_lr = max(self.current_lr * self.factor, self.min_lr)
                if old_lr != self.current_lr:
                    logging.info(
                        f"LR reduced from {old_lr:.6e} to {self.current_lr:.6e} "
                        f"(no improvement for {self.patience} iterations)"
                    )
                self.iterations_since_improvement = 0
                
            return self.current_lr
        
        return self.current_lr
    
    def reset(self):
        """Reset scheduler to initial state."""
        self.current_lr = self.initial_lr
        self.best_metric = float('inf')
        self.iterations_since_improvement = 0


class BaseRicciFlow(ABC):
    """
    Base class for Ricci Flow implementations on graphs.
    
    Implements discrete Ollivier-Ricci curvature for graph edges and
    provides infrastructure for running Ricci Flow iterations.
    """

    def __init__(
        self,
        G: nx.Graph,
        distribution: Callable[[Any], Dict[Any, float]],
        dirsave_graphs: Optional[str] = None,
        overwrite: bool = True,
        save_last: bool = True,
        save_intermediate_graphs: bool = False,
        tol_curvature: float = 1e-11,
        log_level: str = "INFO",
        eps: float = 0.1,
        lr_schedule: Literal["constant", "step", "exponential", "adaptive"] = "constant",
        lr_step_size: int = 10,
        lr_gamma: float = 0.5,
        lr_patience: int = 5,
        lr_factor: float = 0.5,
        min_lr: float = 1e-6,
    ):
        """
        Initialize Ricci Flow.

        Parameters
        ----------
        G : nx.Graph
            Input graph with edge weights
        distribution : Callable
            Function mapping nodes to probability distributions over neighbors
        dirsave_graphs : Optional[str]
            Directory to save intermediate graphs
        overwrite : bool
            Whether to overwrite existing save directory
        save_last : bool
            Whether to save the final graph
        save_intermediate_graphs : bool
            Whether to save graphs at each iteration
        tol_curvature : float
            Tolerance for stopping criterion (when all |K(u,v)| < tol)
        log_level : str
            Logging level (DEBUG, INFO, WARNING, ERROR)
        eps : float
            Initial step size (learning rate) for Ricci Flow updates
        lr_schedule : str
            Learning rate schedule: "constant", "step", "exponential", "adaptive"
        lr_step_size : int
            For step schedule: reduce LR every lr_step_size iterations
        lr_gamma : float
            Multiplicative factor for LR decay
        lr_patience : int
            For adaptive: iterations to wait before reducing LR
        lr_factor : float
            For adaptive: factor to reduce LR by
        min_lr : float
            Minimum learning rate
        """
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format="[Ricci-Flow] %(asctime)s.%(msecs)03d | %(message)s",
            datefmt="%Y-%m-%d@%H:%M:%S",
        )

        self.G = G
        self.distribution_nodes = distribution
        self.dirsave = Path(dirsave_graphs) if dirsave_graphs else None
        self.save_last = save_last
        self.save_intermediate_graphs = save_intermediate_graphs
        self.tol = tol_curvature
        self._counter_iters = 0
        
        # Initialize learning rate scheduler
        self.lr_scheduler = LearningRateScheduler(
            initial_lr=eps,
            schedule=lr_schedule,
            step_size=lr_step_size,
            gamma=lr_gamma,
            patience=lr_patience,
            factor=lr_factor,
            min_lr=min_lr,
        )

        # Validate save directory
        if not overwrite and self.dirsave and self.dirsave.exists():
            raise ValueError(
                f"{self.dirsave} already exists. "
                "Set 'overwrite=True' or use another directory."
            )

        # Create save directory
        if self.dirsave:
            self.dirsave.mkdir(exist_ok=True, parents=True)

        self.initialize_edges()

    def run(self, iterations: int, name: Optional[str] = None) -> nx.Graph:
        """
        Run Ricci Flow for specified number of iterations.

        Parameters
        ----------
        iterations : int
            Number of iterations to run
        name : Optional[str]
            Name prefix for saved graphs

        Returns
        -------
        nx.Graph
            Graph after Ricci Flow
        """
        name = name or "graph"

        for it in track(
            range(iterations),
            total=iterations,
            description="Ricci-Flow",
            transient=False,
        ):
            self._counter_iters += 1
            
            # Get current learning rate
            max_curvature = self.get_max_curvature()
            eps = self.lr_scheduler.step(it, metric=max_curvature)
            
            logging.debug(
                f"Iteration {self._counter_iters}, "
                f"LR={eps:.6e}, max|κ|={max_curvature:.6e}"
            )

            # Compute new weights without modifying graph
            new_weights = self.iter_ricci_flow(eps)

            # Update edge weights
            nx.set_edge_attributes(self.G, new_weights)

            # Update curvatures with new weights
            for edge in self.G.edges:
                self.G.edges[edge]["curvature"] = self.compute_curvature(edge)

            # Save checkpoints
            if self.dirsave:
                self.checkpoints(it, name)

            # Check stopping criterion
            if self.is_curvature_below_tol():
                logging.info(
                    f"Stopping Ricci-Flow at iteration {it+1}: "
                    f"all curvatures below tol={self.tol}"
                )
                break

        return self.G

    def get_max_curvature(self) -> float:
        """
        Get maximum absolute curvature across all edges.
        
        Returns
        -------
        float
            max |κ(u,v)|
        """
        curvatures = [abs(data["curvature"]) for _, _, data in self.G.edges(data=True)]
        return max(curvatures) if curvatures else 0.0

    @abstractmethod
    def iter_ricci_flow(self, eps: float) -> Dict[Tuple[Any, Any], Dict[str, float]]:
        """
        Compute one iteration of Ricci Flow.
        
        Parameters
        ----------
        eps : float
            Current learning rate (step size)
        
        Returns
        -------
        Dict[Tuple[Any, Any], Dict[str, float]]
            Dictionary mapping edges to new weight attributes
        """
        pass

    def compute_curvature(self, edge: Tuple[Any, Any]) -> float:
        """
        Compute Ollivier-Ricci curvature of an edge.

        The curvature is defined as: κ(u,v) = 1 - W(m_u, m_v) / d(u,v)
        where W is the Wasserstein distance and d is the edge weight.

        Parameters
        ----------
        edge : Tuple[Any, Any]
            Edge (u, v) to compute curvature for

        Returns
        -------
        float
            Curvature value
        """
        node1, node2 = edge
        d_1 = self.distribution_nodes(node1)
        d_2 = self.distribution_nodes(node2)

        W = self.wasserstein(d_1, d_2)
        d = self.G.edges[edge]["weight"]

        logging.debug(f"Edge {edge}: W={W:.6f}, d={d:.6f}")

        return 1 - W / d

    def wasserstein(
        self,
        distribution_node1: Dict[Any, float],
        distribution_node2: Dict[Any, float],
    ) -> float:
        """
        Compute Wasserstein distance between two node distributions.

        Parameters
        ----------
        distribution_node1 : Dict[Any, float]
            Probability distribution over neighbors of first node
        distribution_node2 : Dict[Any, float]
            Probability distribution over neighbors of second node

        Returns
        -------
        float
            Wasserstein-1 distance between distributions
        """
        # Extract subgraph with nodes from both distributions
        nodes_subgraph = list(distribution_node1.keys()) + list(
            distribution_node2.keys()
        )
        subgraph = self.G.subgraph(nodes_subgraph).to_undirected()

        # Compute shortest paths
        try:
            distances_subgraph = dict(
                nx.all_pairs_dijkstra_path(subgraph, weight="weight")
            )
        except nx.NetworkXNoPath as e:
            logging.warning(f"No path found in subgraph: {e}")
            return float("inf")

        # Build cost matrix for optimal transport
        M = np.zeros((len(distribution_node1), len(distribution_node2)))

        for i, source in enumerate(distribution_node1.keys()):
            for j, target in enumerate(distribution_node2.keys()):
                if source == target:
                    M[i, j] = 0.0
                    continue

                try:
                    path = distances_subgraph[source][target]
                    edges = [(n1, n2) for n1, n2 in zip(path[:-1], path[1:])]
                    M[i, j] = sum(self.G.edges[e]["weight"] for e in edges)
                except (KeyError, nx.NetworkXNoPath):
                    # No path exists, set to large value
                    M[i, j] = np.inf

        # Normalize cost matrix to [0, 1] if not all infinite
        max_val = M[np.isfinite(M)].max() if np.any(np.isfinite(M)) else 1.0
        if max_val > 0:
            M = np.where(np.isfinite(M), M / max_val, 1.0)

        # Prepare probability vectors
        a = np.array(list(distribution_node1.values()))
        b = np.array(list(distribution_node2.values()))

        # Normalize to ensure they sum to 1
        a = a / a.sum()
        b = b / b.sum()

        # Compute Wasserstein distance
        return ot.emd2(a, b, M)

    def is_curvature_below_tol(self) -> bool:
        """
        Check if all edge curvatures are below tolerance.

        Returns
        -------
        bool
            True if all |κ(u,v)| < tol
        """
        curvatures = [data["curvature"] for _, _, data in self.G.edges(data=True)]
        mean_curv = np.mean(curvatures)
        std_curv = np.std(curvatures)
        max_curv = np.max(np.abs(curvatures))
        
        logging.debug(
            f"Curvature stats: mean={mean_curv:.6e}, "
            f"std={std_curv:.6e}, max|κ|={max_curv:.6e}"
        )
        
        for u, v, data in self.G.edges(data=True):
            if np.abs(data["curvature"]) > self.tol:
                return False
        return True

    def checkpoints(self, iteration: int, name: str) -> None:
        """
        Save graph checkpoints.

        Parameters
        ----------
        iteration : int
            Current iteration number
        name : str
            Name prefix for saved files
        """
        if self.save_intermediate_graphs:
            path_save = self.dirsave / f"{name}-ricciflow-iter{iteration+1}.edgelist"
            nx.write_edgelist(self.G, path_save, data=True)
            logging.debug(f"Saved intermediate graph to {path_save}")

        if self.save_last:
            path_save = self.dirsave / f"{name}-ricciflow-final.edgelist"
            nx.write_edgelist(self.G, path_save, data=True)

    def initialize_edges(self) -> None:
        """
        Initialize edge weights and curvatures.
        
        Sets default weights (uniform) if not present, and computes
        initial curvatures for all edges.
        """
        n_edges = len(self.G.edges())
        weight_init = 1.0 / n_edges

        # Initialize weights
        for edge in self.G.edges():
            if "weight" not in self.G.edges[edge]:
                self.G.edges[edge]["weight"] = weight_init

        # Compute initial curvatures
        logging.info("Computing initial curvatures...")
        for edge in self.G.edges():
            self.G.edges[edge]["curvature"] = self.compute_curvature(edge)


class RicciFlow(BaseRicciFlow):
    """
    Standard Ricci Flow implementation.

    Updates edge weights according to:
    w^{(t+1)}_{(u,v)} = w^{(t)}_{(u,v)} - ε * κ^{(t)}_{(u,v)} * w^{(t)}_{(u,v)}
    
    Converges to uniform curvature κ = 0 (flat geometry).
    """

    def iter_ricci_flow(self, eps: float) -> Dict[Tuple[Any, Any], Dict[str, float]]:
        """
        Compute one iteration of standard Ricci Flow.

        Parameters
        ----------
        eps : float
            Current learning rate (step size)

        Returns
        -------
        Dict[Tuple[Any, Any], Dict[str, float]]
            Dictionary mapping edges to new weights
        """
        logging.debug(f"Standard Ricci-Flow iteration {self._counter_iters}")

        new_weights = {}
        for edge in self.G.edges():
            current_curvature = self.G.edges[edge]["curvature"]
            current_weight = self.G.edges[edge]["weight"]

            # Standard Ricci Flow update
            new_weight = current_weight - eps * current_curvature * current_weight

            # Ensure weight stays positive
            new_weight = max(new_weight, 1e-10)

            new_weights[edge] = {"weight": new_weight}

        return new_weights


class NormalizedRicciFlow(BaseRicciFlow):
    """
    Normalized Ricci Flow implementation with volume preservation.

    IMPORTANT: The formula w^(t+1) = w^(t) - ε*(κ - (1/σ)*Σκw)*w converges to 
    UNIFORM curvature (all κ equal), NOT κ=0. The equilibrium is when all curvatures 
    equal the weighted average: κ(u,v) = (1/σ)*Σκw for all edges.
    
    To preserve volume AND converge to κ=0, we use a two-step approach:
    1. Apply standard Ricci flow: w' = w - ε*κ*w
    2. Rescale to preserve total weight: w_new = w' * (σ_old / σ_new)
    
    This combines the standard flow's κ→0 convergence with volume preservation.
    """
    
    def __init__(self, *args, use_rescaling: bool = True, **kwargs):
        """
        Initialize Normalized Ricci Flow.
        
        Parameters
        ----------
        use_rescaling : bool
            If True, use rescaling method (converges to κ=0 with volume preservation)
            If False, use subtraction method (converges to uniform curvature ≠ 0)
        """
        self.use_rescaling = use_rescaling
        super().__init__(*args, **kwargs)

    def iter_ricci_flow(self, eps: float) -> Dict[Tuple[Any, Any], Dict[str, float]]:
        """
        Compute one iteration of normalized Ricci Flow.

        Parameters
        ----------
        eps : float
            Current learning rate (step size)

        Returns
        -------
        Dict[Tuple[Any, Any], Dict[str, float]]
            Dictionary mapping edges to new weights
        """
        logging.debug(f"Normalized Ricci-Flow iteration {self._counter_iters}")

        # Compute total weight σ = Σ w(i,j)
        sigma_old = sum(self.G.edges[edge]["weight"] for edge in self.G.edges())
        
        if self.use_rescaling:
            # METHOD 1: Rescaling approach (converges to κ=0)
            # Step 1: Apply standard Ricci flow
            new_weights = {}
            for edge in self.G.edges():
                current_curvature = self.G.edges[edge]["curvature"]
                current_weight = self.G.edges[edge]["weight"]
                
                # Standard update: w' = w - ε*κ*w = w*(1 - ε*κ)
                new_weight = current_weight * (1 - eps * current_curvature)
                new_weight = max(new_weight, 1e-10)
                new_weights[edge] = {"weight": new_weight}
            
            # Step 2: Rescale to preserve total weight
            sigma_new = sum(w["weight"] for w in new_weights.values())
            rescale_factor = sigma_old / sigma_new if sigma_new > 0 else 1.0
            
            for edge in new_weights:
                new_weights[edge]["weight"] *= rescale_factor
            
            logging.debug(
                f"Rescaling method: σ_old={sigma_old:.6f}, σ_new_before={sigma_new:.6f}, "
                f"rescale={rescale_factor:.6f}"
            )
            
        else:
            # METHOD 2: Subtraction approach (converges to uniform κ ≠ 0)
            # This is what you provided but it WON'T converge to κ=0!
            
            # Compute weighted curvature sum: Σ κ(i,j) * w(i,j)
            weighted_curvature_sum = sum(
                self.G.edges[edge]["curvature"] * self.G.edges[edge]["weight"]
                for edge in self.G.edges()
            )
            
            # Compute normalization term: (1/σ) * Σ κ(i,j) * w(i,j)
            normalization_term = weighted_curvature_sum / sigma_old if sigma_old > 0 else 0.0
            
            logging.debug(
                f"Subtraction method: σ={sigma_old:.6f}, "
                f"Σ(κ*w)={weighted_curvature_sum:.6f}, "
                f"norm_term={normalization_term:.6f}"
            )

            # Update all edges
            new_weights = {}
            for edge in self.G.edges():
                current_curvature = self.G.edges[edge]["curvature"]
                current_weight = self.G.edges[edge]["weight"]

                # w^(t+1) = w^(t) - ε * (κ - norm_term) * w^(t)
                delta = (current_curvature - normalization_term) * current_weight
                new_weight = current_weight - eps * delta
                new_weight = max(new_weight, 1e-10)
                new_weights[edge] = {"weight": new_weight}
        
        # Verify total weight
        new_total = sum(w["weight"] for w in new_weights.values())
        logging.debug(f"Weight conservation: σ_old={sigma_old:.6f}, σ_new={new_total:.6f}")

        return new_weights
