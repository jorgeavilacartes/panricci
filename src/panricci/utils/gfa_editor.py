from typing import Optional, List, Tuple, Dict, Set
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt


class GFAGraph:
    """
    Comprehensive utility for loading, modifying, and saving GFA (Graph Fragment Assembly) files.
    
    Supports:
    - Loading GFA files into memory
    - Removing nodes (with predecessor-successor reconnection)
    - Removing edges (with node collapsing and label concatenation)
    - Removing paths
    - Saving modified GFA files
    """
    
    def __init__(self, gfa_file: Optional[str] = None):
        """
        Initialize GFA graph.
        
        Parameters
        ----------
        gfa_file : Optional[str]
            Path to GFA file to load
        """
        self.header = "H\tVN:Z:1.1"
        self.nodes: Dict[str, str] = {}  # nodeid -> sequence
        self.edges: List[Tuple[str, str]] = []  # list of (from_node, to_node)
        self.paths: Dict[str, List[Tuple[str, str]]] = {}  # pathid -> [(nodeid, orientation)]
        
        if gfa_file:
            self.load(gfa_file)
    
    def load(self, filename: str) -> None:
        """
        Load a GFA file.
        
        Parameters
        ----------
        filename : str
            Path to GFA file
        """
        self.nodes.clear()
        self.edges.clear()
        self.paths.clear()
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t')
                record_type = parts[0]
                
                if record_type == 'H':
                    self.header = line
                    
                elif record_type == 'S':
                    # S	nodeid	sequence	[optional fields]
                    nodeid = parts[1]
                    sequence = parts[2]
                    self.nodes[nodeid] = sequence
                    
                elif record_type == 'L':
                    # L	from	from_orient	to	to_orient	overlap
                    from_node = parts[1]
                    to_node = parts[3]
                    self.edges.append((from_node, to_node))
                    
                elif record_type == 'P':
                    # P	path_name	node_list	overlaps
                    path_name = parts[1]
                    path_string = parts[2]
                    self.paths[path_name] = self._parse_path(path_string)
    
    def _parse_path(self, path_str: str) -> List[Tuple[str, str]]:
        """
        Parse path string into list of (nodeid, orientation) tuples.
        
        Parameters
        ----------
        path_str : str
            Path string like "1+,2+,4+" or "1-,2+,3-"
        
        Returns
        -------
        List[Tuple[str, str]]
            List of (nodeid, orientation) tuples
        """
        if not path_str or path_str == '*':
            return []
        
        nodes = []
        for node in path_str.split(','):
            node = node.strip()
            if node.endswith('+'):
                nodes.append((node[:-1], '+'))
            elif node.endswith('-'):
                nodes.append((node[:-1], '-'))
            else:
                # Default to forward orientation if not specified
                nodes.append((node, '+'))
        return nodes
    
    def _format_path(self, path: List[Tuple[str, str]]) -> str:
        """
        Format path list back to string.
        
        Parameters
        ----------
        path : List[Tuple[str, str]]
            List of (nodeid, orientation) tuples
        
        Returns
        -------
        str
            Formatted path string
        """
        if not path:
            return '*'
        return ','.join([f"{nodeid}{orient}" for nodeid, orient in path])
    
    def save(self, filename: str) -> None:
        """
        Save the graph to a GFA file.
        
        Parameters
        ----------
        filename : str
            Output file path
        """
        with open(filename, 'w') as f:
            # Write header
            f.write(self.header + '\n')
            
            # Write nodes (segments)
            for nodeid in sorted(self.nodes.keys(), key=lambda x: int(x) if x.isdigit() else x):
                seq = self.nodes[nodeid]
                f.write(f"S\t{nodeid}\t{seq}\n")
            
            # Write edges (links)
            for from_node, to_node in self.edges:
                f.write(f"L\t{from_node}\t+\t{to_node}\t+\t0M\n")
            
            # Write paths
            for path_name in sorted(self.paths.keys()):
                path_str = self._format_path(self.paths[path_name])
                f.write(f"P\t{path_name}\t{path_str}\t*\n")
    
    def get_predecessors(self, nodeid: str) -> List[str]:
        """Get all predecessor nodes."""
        return [from_node for from_node, to_node in self.edges if to_node == nodeid]
    
    def get_successors(self, nodeid: str) -> List[str]:
        """Get all successor nodes."""
        return [to_node for from_node, to_node in self.edges if from_node == nodeid]
    
    def remove_node(self, nodeid: str, reconnect: bool = True) -> None:
        """
        Remove a node from the graph.
        
        Parameters
        ----------
        nodeid : str
            Node identifier to remove
        reconnect : bool
            If True, connect all predecessors to all successors
        """
        if nodeid not in self.nodes:
            print(f"Warning: Node {nodeid} not found")
            return
        
        # Get predecessors and successors
        predecessors = self.get_predecessors(nodeid)
        successors = self.get_successors(nodeid)
        
        # Remove edges involving this node
        self.edges = [(f, t) for f, t in self.edges if f != nodeid and t != nodeid]
        
        # Reconnect if requested
        if reconnect:
            for pred in predecessors:
                for succ in successors:
                    if (pred, succ) not in self.edges:
                        self.edges.append((pred, succ))
        
        # Remove the node
        del self.nodes[nodeid]
        
        # Update all paths to remove this node
        for path_name in self.paths:
            self.paths[path_name] = [
                (nid, orient) for nid, orient in self.paths[path_name] 
                if nid != nodeid
            ]
        
        print(f"Removed node {nodeid}")
    
    def remove_edge(self, from_node: str, to_node: str, collapse_nodes: bool = True) -> None:
        """
        Remove an edge and optionally collapse the two nodes into one.
        
        When collapsing:
        - Creates new node with concatenated sequence
        - All edges to from_node become edges to new node
        - All edges from to_node become edges from new node
        - Paths are updated to use the new node
        
        Parameters
        ----------
        from_node : str
            Source node of edge
        to_node : str
            Target node of edge
        collapse_nodes : bool
            If True, collapse both nodes into a single node with concatenated label
        """
        if (from_node, to_node) not in self.edges:
            print(f"Warning: Edge {from_node} -> {to_node} not found")
            return
        
        # Remove the edge
        self.edges.remove((from_node, to_node))
        
        if collapse_nodes:
            # Create new node with concatenated sequence
            new_nodeid = f"{from_node}_{to_node}"
            new_sequence = self.nodes[from_node] + self.nodes[to_node]
            self.nodes[new_nodeid] = new_sequence
            
            # Redirect all edges
            # Edges coming into from_node now come into new_node
            # Edges going out of to_node now go out of new_node
            new_edges = []
            for f, t in self.edges:
                if f == from_node and t == to_node:
                    continue  # Skip the removed edge
                elif t == from_node:
                    new_edges.append((f, new_nodeid))
                elif f == from_node:
                    new_edges.append((new_nodeid, t))
                elif t == to_node:
                    new_edges.append((f, new_nodeid))
                elif f == to_node:
                    new_edges.append((new_nodeid, t))
                else:
                    new_edges.append((f, t))
            
            self.edges = new_edges
            
            # Update paths
            for path_name in self.paths:
                new_path = []
                i = 0
                path = self.paths[path_name]
                
                while i < len(path):
                    curr_nodeid, curr_orient = path[i]
                    
                    # Check if this is the start of the collapsed sequence
                    if (i < len(path) - 1 and 
                        curr_nodeid == from_node and 
                        path[i + 1][0] == to_node):
                        # Replace both nodes with the new collapsed node
                        new_path.append((new_nodeid, curr_orient))
                        i += 2  # Skip both nodes
                    elif curr_nodeid == from_node or curr_nodeid == to_node:
                        # Single occurrence, replace with new node
                        new_path.append((new_nodeid, curr_orient))
                        i += 1
                    else:
                        new_path.append((curr_nodeid, curr_orient))
                        i += 1
                
                self.paths[path_name] = new_path
            
            # Remove old nodes
            if from_node in self.nodes:
                del self.nodes[from_node]
            if to_node in self.nodes:
                del self.nodes[to_node]
            
            print(f"Removed edge {from_node} -> {to_node} and collapsed into node {new_nodeid}")
        else:
            print(f"Removed edge {from_node} -> {to_node} (no collapsing)")
            
            # Update paths that used this edge - remove the edge portion
            for path_name in list(self.paths.keys()):
                path = self.paths[path_name]
                new_path = []
                
                for i in range(len(path)):
                    curr_nodeid, curr_orient = path[i]
                    
                    # Check if this edge was part of the path
                    if (i < len(path) - 1 and 
                        curr_nodeid == from_node and 
                        path[i + 1][0] == to_node):
                        # Edge was in path - keep only the from_node
                        new_path.append((curr_nodeid, curr_orient))
                        # Skip to_node
                        continue
                    
                    new_path.append((curr_nodeid, curr_orient))
                
                self.paths[path_name] = new_path
    
    def remove_path(self, path_name: str, remove_unused_nodes: bool = False) -> None:
        """
        Remove a path from the graph.
        
        Parameters
        ----------
        path_name : str
            Path identifier to remove
        remove_unused_nodes : bool
            If True, also remove nodes that are not used by any remaining path
        """
        if path_name not in self.paths:
            print(f"Warning: Path {path_name} not found")
            return
        
        # Remove the path
        del self.paths[path_name]
        print(f"Removed path {path_name}")
        
        if remove_unused_nodes:
            # Find nodes used by remaining paths
            used_nodes = set()
            for path in self.paths.values():
                for nodeid, _ in path:
                    used_nodes.add(nodeid)
            
            # Find and remove unused nodes
            unused_nodes = [nid for nid in self.nodes.keys() if nid not in used_nodes]
            
            for nodeid in unused_nodes:
                self.remove_node(nodeid, reconnect=False)
            
            if unused_nodes:
                print(f"Removed {len(unused_nodes)} unused nodes: {unused_nodes}")
    
    def get_sources_and_sinks(self) -> Tuple[Set[str], Set[str]]:
        """
        Get source and sink nodes from paths.
        
        Returns
        -------
        Tuple[Set[str], Set[str]]
            (sources, sinks) - sets of node IDs
        """
        sources = set()
        sinks = set()
        
        for path in self.paths.values():
            if path:
                sources.add(path[0][0])   # First node in path
                sinks.add(path[-1][0])    # Last node in path
        
        return sources, sinks
    
    def to_networkx(self, directed: bool = True) -> nx.Graph:
        """
        Convert GFA graph to NetworkX graph.
        
        Parameters
        ----------
        directed : bool
            If True, create directed graph; otherwise undirected
        
        Returns
        -------
        nx.Graph
            NetworkX graph representation
        """
        G = nx.DiGraph() if directed else nx.Graph()
        
        # Add nodes with attributes
        for nodeid, sequence in self.nodes.items():
            G.add_node(nodeid, label=sequence, length=len(sequence))
        
        # Add edges
        G.add_edges_from(self.edges)
        
        return G
    
    def compute_node_depth(self) -> Dict[str, float]:
        """
        Compute node depth: fraction of paths that pass through each node.
        
        Returns
        -------
        Dict[str, float]
            nodeid -> depth (0.0 to 1.0)
        """
        if not self.paths:
            return {}
        
        paths_by_node = defaultdict(set)
        
        for path_name, path in self.paths.items():
            for nodeid, _ in path:
                paths_by_node[nodeid].add(path_name)
        
        n_paths = len(self.paths)
        return {
            nodeid: len(path_set) / n_paths 
            for nodeid, path_set in paths_by_node.items()
        }
    
    def draw(
        self, 
        layout: str = 'hierarchical',
        figsize: Tuple[int, int] = (14, 8),
        node_size: int = 3000,
        font_size: int = 10,
        show_sequence: bool = True,
        max_seq_length: int = 8,
        color_by_depth: bool = True,
        save_path: Optional[str] = None,
        dpi: int = 150,
    ) -> None:
        """
        Draw the GFA graph using networkx and matplotlib.
        
        Parameters
        ----------
        layout : str
            Layout algorithm: 'hierarchical', 'spring', 'circular', 'kamada_kawai', 'planar'
        figsize : Tuple[int, int]
            Figure size (width, height)
        node_size : int
            Size of nodes
        font_size : int
            Font size for labels
        show_sequence : bool
            If True, show sequences in node labels; otherwise just node IDs
        max_seq_length : int
            Maximum sequence length to display (truncate longer sequences)
        color_by_depth : bool
            If True, color nodes by their depth (fraction of paths using them)
        save_path : Optional[str]
            If provided, save figure to this path
        dpi : int
            DPI for saved figure
        """
        # Create NetworkX graph
        G = self.to_networkx(directed=True)
        
        # Compute layout
        if layout == 'hierarchical':
            pos = self._hierarchical_layout(G)
                
        elif layout == 'spring':
            pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        elif layout == 'planar':
            try:
                pos = nx.planar_layout(G)
            except:
                print("Graph is not planar, using spring layout")
                pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
        else:
            raise ValueError(f"Unknown layout: {layout}")
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Prepare node colors based on depth
        if color_by_depth:
            node_depths = self.compute_node_depth()
            node_colors = [node_depths.get(node, 0.0) for node in G.nodes()]
            cmap = plt.cm.YlOrRd
            vmin, vmax = 0, 1
        else:
            node_colors = 'lightblue'
            cmap = None
            vmin, vmax = None, None
        
        # Draw nodes first
        nx.draw_networkx_nodes(
            G, pos,
            node_color=node_colors,
            node_size=node_size,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            edgecolors='black',
            linewidths=2,
            ax=ax
        )
        
        # Draw edges with proper arrows - use FancyArrowPatch for better control
        from matplotlib.patches import FancyArrowPatch, ArrowStyle
        
        for edge in G.edges():
            x1, y1 = pos[edge[0]]
            x2, y2 = pos[edge[1]]
            
            arrow = FancyArrowPatch(
                (x1, y1), (x2, y2),
                arrowstyle=ArrowStyle.Simple(head_length=15, head_width=15, tail_width=2),
                color='gray',
                alpha=0.6,
                connectionstyle="arc3,rad=0.1",
                zorder=1,
                linewidth=2
            )
            ax.add_patch(arrow)
        
        # Prepare node labels
        labels = {}
        for node in G.nodes():
            if show_sequence:
                seq = self.nodes[node]
                if len(seq) > max_seq_length:
                    seq = seq[:max_seq_length] + '...'
                labels[node] = f"{node}\n{seq}"
            else:
                labels[node] = node
        
        # Draw labels
        nx.draw_networkx_labels(
            G, pos,
            labels=labels,
            font_size=font_size,
            font_weight='bold',
            ax=ax
        )
        
        # Add colorbar if coloring by depth
        if color_by_depth and len(G.nodes()) > 0:
            sm = plt.cm.ScalarMappable(
                cmap=cmap,
                norm=plt.Normalize(vmin=0, vmax=1)
            )
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('Node Depth (Fraction of Paths)', rotation=270, labelpad=20)
        
        # Styling
        ax.set_title(
            f"GFA Graph: {len(self.nodes)} nodes, {len(self.edges)} edges, {len(self.paths)} paths",
            fontsize=14,
            fontweight='bold',
            pad=20
        )
        ax.axis('off')
        plt.tight_layout()
        
        # Save or show
        if save_path:
            plt.savefig(save_path, dpi=dpi, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
        
        plt.show()
    
    def _hierarchical_layout(self, G: nx.DiGraph) -> Dict[str, Tuple[float, float]]:
        """
        Create hierarchical layout that shows bubble structure.
        Uses BFS to assign levels and spreads nodes vertically within levels.
        
        Parameters
        ----------
        G : nx.DiGraph
            NetworkX directed graph
        
        Returns
        -------
        Dict[str, Tuple[float, float]]
            Node positions {node: (x, y)}
        """
        from collections import deque
        
        # Find source nodes (nodes with no predecessors)
        sources = [n for n in G.nodes() if G.in_degree(n) == 0]
        
        if not sources:
            # If no clear sources (cyclic graph), use nodes from paths
            path_sources, _ = self.get_sources_and_sinks()
            sources = list(path_sources) if path_sources else [list(G.nodes())[0]]
        
        # BFS to assign levels (x-coordinates)
        levels = {}
        queue = deque([(s, 0) for s in sources])
        visited = set(sources)
        
        for s in sources:
            levels[s] = 0
        
        while queue:
            node, level = queue.popleft()
            
            for successor in G.successors(node):
                if successor not in visited:
                    visited.add(successor)
                    levels[successor] = level + 1
                    queue.append((successor, level + 1))
                else:
                    # Node already visited, but might need higher level
                    levels[successor] = max(levels.get(successor, 0), level + 1)
        
        # Handle any nodes not reached by BFS
        for node in G.nodes():
            if node not in levels:
                levels[node] = max(levels.values()) + 1 if levels else 0
        
        # Group nodes by level
        level_groups = defaultdict(list)
        for node, level in levels.items():
            level_groups[level].append(node)
        
        # Assign positions
        pos = {}
        for level, nodes in level_groups.items():
            # Sort nodes within level for consistency
            sorted_nodes = sorted(nodes, key=lambda x: int(x) if x.isdigit() else x)
            n_nodes = len(sorted_nodes)
            
            # Spread nodes vertically, centered around y=0
            for i, node in enumerate(sorted_nodes):
                x = level
                y = (i - (n_nodes - 1) / 2) * 1.5  # Spread with spacing factor
                pos[node] = (x, y)
        
        return pos
    
    def __str__(self) -> str:
        """String representation of the graph."""
        sources, sinks = self.get_sources_and_sinks()
        return (
            f"GFA Graph:\n"
            f"  Nodes: {len(self.nodes)}\n"
            f"  Edges: {len(self.edges)}\n"
            f"  Paths: {len(self.paths)}\n"
            f"  Sources: {sources}\n"
            f"  Sinks: {sinks}"
        )
    
    def __repr__(self) -> str:
        return self.__str__()


# Example usage and testing
if __name__ == "__main__":
    # Create example GFA content
    example_gfa = """H	VN:Z:1.1
S	1	AAAA
S	2	AAC
S	3	C
S	4	T
S	5	A
S	6	CCTAAA
S	7	AC
S	8	CC
S	9	T
S	10	C
S	11	AA
S	12	CCT
S	13	AC
S	14	T
S	15	C
S	16	AA
S	17	CCT
S	18	AC
L	1	+	2	+	0M
L	1	+	3	+	0M
L	2	+	4	+	0M
L	3	+	4	+	0M
L	4	+	5	+	0M
L	5	+	6	+	0M
L	5	+	7	+	0M
L	5	+	8	+	0M
L	6	+	9	+	0M
L	7	+	9	+	0M
L	8	+	9	+	0M
L	9	+	10	+	0M
L	10	+	11	+	0M
L	10	+	12	+	0M
L	11	+	13	+	0M
L	12	+	13	+	0M
L	13	+	14	+	0M
L	14	+	15	+	0M
L	14	+	16	+	0M
L	14	+	17	+	0M
L	15	+	18	+	0M
L	16	+	18	+	0M
L	17	+	18	+	0M
P	seq1	1+,2+,4+,5+,6+,9+,10+,11+,13+,14+,15+,18+	*
P	seq2	1+,2+,4+,5+,7+,9+,10+,12+,13+,14+,16+,18+	*
P	seq3	1+,3+,4+,5+,8+,9+,10+,11+,13+,14+,17+,18+	*
P	seq4	1+,3+,4+,5+,6+,9+,10+,11+,13+,14+,17+,18+	*"""
    
    # Write example file
    with open('example.gfa', 'w') as f:
        f.write(example_gfa)
    
    print("=" * 70)
    print("EXAMPLE 1: Remove a node")
    print("=" * 70)
    gfa1 = GFAGraph('example.gfa')
    print("Original:")
    print(gfa1)
    print(f"\nNode 6 sequence: {gfa1.nodes['6']}")
    print(f"Node 6 predecessors: {gfa1.get_predecessors('6')}")
    print(f"Node 6 successors: {gfa1.get_successors('6')}")
    
    gfa1.remove_node('6', reconnect=True)
    gfa1.save('example_node_removed.gfa')
    print("\nAfter removing node 6:")
    print(gfa1)
    print("Saved to: example_node_removed.gfa\n")
    
    print("=" * 70)
    print("EXAMPLE 2: Remove an edge (with node collapsing)")
    print("=" * 70)
    gfa2 = GFAGraph('example.gfa')
    print("Original:")
    print(f"Node 4 sequence: {gfa2.nodes['4']}")
    print(f"Node 5 sequence: {gfa2.nodes['5']}")
    print(f"Edge 4->5 exists: {(('4', '5') in gfa2.edges)}")
    
    gfa2.remove_edge('4', '5', collapse_nodes=True)
    print(f"\nAfter collapsing edge 4->5:")
    print(f"New node 4_5 sequence: {gfa2.nodes.get('4_5', 'Not found')}")
    print(f"Path seq1: {gfa2._format_path(gfa2.paths['seq1'])}")
    gfa2.save('example_edge_collapsed.gfa')
    print("Saved to: example_edge_collapsed.gfa\n")
    
    print("=" * 70)
    print("EXAMPLE 3: Remove a path")
    print("=" * 70)
    gfa3 = GFAGraph('example.gfa')
    print("Original paths:", list(gfa3.paths.keys()))
    print("Original nodes:", len(gfa3.nodes))
    
    gfa3.remove_path('seq1', remove_unused_nodes=False)
    print(f"\nAfter removing seq1:")
    print("Remaining paths:", list(gfa3.paths.keys()))
    print("Nodes:", len(gfa3.nodes))
    gfa3.save('example_path_removed.gfa')
    print("Saved to: example_path_removed.gfa\n")
    
    print("=" * 70)
    print("EXAMPLE 4: Node depth computation")
    print("=" * 70)
    gfa4 = GFAGraph('example.gfa')
    node_depths = gfa4.compute_node_depth()
    print("Node depths (fraction of paths passing through):")
    for nodeid in sorted(node_depths.keys(), key=lambda x: int(x) if x.isdigit() else x):
        print(f"  Node {nodeid}: {node_depths[nodeid]:.2f}")
    
    print("\n" + "=" * 70)
    print("EXAMPLE 5: Visualize the GFA graph")
    print("=" * 70)
    gfa_viz = GFAGraph('example.gfa')
    print("Drawing original graph with hierarchical layout...")
    gfa_viz.draw(
        layout='hierarchical',
        figsize=(16, 10),
        node_size=3000,
        font_size=9,
        show_sequence=True,
        max_seq_length=6,
        color_by_depth=True,
        save_path='gfa_visualization.png'
    )
    
    print("\n" + "=" * 70)
    print("EXAMPLE 6: Visualize after modifications")
    print("=" * 70)
    gfa_modified = GFAGraph('example.gfa')
    print("Removing edge 10->11 with collapsing...")
    gfa_modified.remove_edge('10', '11', collapse_nodes=True)
    print("Drawing modified graph...")
    gfa_modified.draw(
        layout='hierarchical',
        figsize=(16, 10),
        node_size=3000,
        font_size=9,
        show_sequence=True,
        color_by_depth=True,
        save_path='gfa_modified_visualization.png'
    )