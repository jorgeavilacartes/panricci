MARKDOWN = """
# `panricci` WHAT-TO-DO?   
- `panricci` is a tool for aligning pangenome graphs (sequence and variation graphs). 

- Graphs are manifolds, and they are evolved until the notion of curvature is constant. 
- This state gave us a metric (weights for the arcs in the graph) which depends on the node information
which can be, the paths, the neighbors, the label, etc, everything you can model as a probabity distribution over the
neighbors of a node. 

- `panricci` uses the Ricci-Flow (and Normalized Ricci-Flow) to evolve the graphs. 

- The alignment process is independent, it only requires the final state of the graph
to create vector representations of each node.

A step-by-step guide to help you

1. Evolve your graph using Ricci-Flow | `panspace ricci-flow --help`
2. Do the same for another (or the same) graph | `panspace ricci-flow --help`
3. Align them: are they isomorphic? | `panspace alignment --help`
"""

import typer
import logging
from pathlib import Path
from typing_extensions import Annotated
# types for typer
from rich.progress import track
from rich import print 
from rich.console import Console
from rich.markdown import Markdown

console = Console()
app = typer.Typer(name="PanRicci", rich_markup_mode="rich",
                  help="""
                  :cat: Welcome to [white bold]PanRicci[/white bold] :  Alignment of pangenome graphs with Ricci-Flow
                  """
                  )
# @app.command("docs", help="Open documentation webpage.")
# def github() -> None:
#     typer.launch("https://github.com/jorgeavilacartes/panricci")

@app.command("ricci-flow", help="apply ricci flow to a graph.")
def ricci_flow(
    gfa: Annotated[Path, typer.Option("--gfa", "-g", help="Path to the GFA file.")], 
    iterations: Annotated[int, typer.Option("--iterations", "-i", help="Maximum number of iterations to run Ricci-Flow.")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory to save .")] = "output-ricci-flow/ricci-graph", 
    tol_curvature: Annotated[float, typer.Option("--tol-curvature", "-t", help="Tolerance for curvature. If all curvatures are smaller than this, then the algorithm stop.")] = 1e-11,
    sequence_graph: Annotated[bool, typer.Option("--sequence-graph", "-s", help="If set, define node distributions as a sequence graph, otherwise use the one for variation graphs (considering paths).")] = False,
    log_level: Annotated[str, typer.Option("--log-level", "-l", help="Log level. Default: INFO.")] = "INFO",
    save_intermediate_graphs: Annotated[bool, typer.Option("--save-intermediate-graphs", "-si", help="Save intermediate graphs.")] = False,
    undirected: Annotated[bool, typer.Option("--undirected", "-u", help="Load GFA as undirected graph")] = False,
    ):
    from pathlib import Path
    from panricci.ricci_flow import RicciFlow
    from panricci.utils import GFALoader
    
    if sequence_graph:
        from panricci.node_distributions.sequence_graph import DistributionNodes
    else:
        from panricci.node_distributions.variation_graph import DistributionNodes
    
    dirsave    = Path(outdir); dirsave.mkdir(exist_ok=True, parents=True)
    gfa_loader = GFALoader(undirected=undirected)

    # load graph
    G = gfa_loader(gfa)

    # compute distribution of nodes of the graph
    distribution = DistributionNodes(G, alpha=0.5)

    # run Ricci-Flow
    ricci_flow = RicciFlow(G, 
                    distribution, 
                    save_last=True, 
                    save_intermediate_graphs=save_intermediate_graphs, 
                    dirsave_graphs=dirsave,
                    tol_curvature=tol_curvature,
                    log_level=log_level,
                    )

    G_ricci = ricci_flow.run(iterations=iterations, name=Path(outdir).stem)

    return G_ricci


@app.command("ricci-flow2", help="apply ricci flow to a graph.")
def ricci_flow2(
    gfa: Annotated[Path, typer.Option("--gfa", "-g", help="Path to the GFA file.")], 
    iterations: Annotated[int, typer.Option("--iterations", "-i", help="Maximum number of iterations to run Ricci-Flow.")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory to save .")] = "output-ricci-flow/ricci-graph", 
    tol_curvature: Annotated[float, typer.Option("--tol-curvature", "-t", help="Tolerance for curvature. If all curvatures are smaller than this, then the algorithm stop.")] = 1e-11,
    sequence_graph: Annotated[bool, typer.Option("--sequence-graph", "-s", help="If set, define node distributions as a sequence graph, otherwise use the one for variation graphs (considering paths).")] = False,
    log_level: Annotated[str, typer.Option("--log-level", "-l", help="Log level. Default: INFO/DEBUG")] = "INFO",
    undirected: Annotated[bool, typer.Option("--undirected", "-u", help="Load GFA as undirected graph")] = False,
    eps: Annotated[float, typer.Option("--eps", help="Step size for curvature to update weights")] = 0.1,
    min_weight: Annotated[float, typer.Option("--min-weight", help="minimum weight for edge")] = 0.000001,
    ):
    from pathlib import Path
    from panricci.ricci_flow import RicciFlow2 as RicciFlow
    from panricci.utils import GFALoader
    
    if sequence_graph:
        from panricci.node_distributions.sequence_graph import DistributionNodes
    else:
        from panricci.node_distributions.variation_graph import DistributionNodes
    
    dirsave    = Path(outdir); dirsave.mkdir(exist_ok=True, parents=True)
    gfa_loader = GFALoader(undirected=undirected)

    # load graph
    G = gfa_loader(gfa)

    # compute distribution of nodes of the graph
    distribution = DistributionNodes(G, alpha=0.5)

    # run RicciFlow
    ricci = RicciFlow(
        G, 
        distribution,
        dirsave_graphs=outdir,
        step_size=eps,      # Start with 0.1 instead of 1.0
        min_weight=min_weight,    # Prevent weights from going to zero
        log_level=log_level,   # See warnings
        tol_curvature=tol_curvature,
    )

    G_ricci = ricci.run(iterations=iterations)
    return G_ricci

@app.command("normalized-ricci-flow", help="Normalized Ricci Flow")
def normalized_ricci_flow(
    gfa: Annotated[Path, typer.Option("--gfa", "-g", help="Path to the GFA file.")], 
    iterations: Annotated[int, typer.Option("--iterations", "-i", help="Maximum number of iterations to run Ricci-Flow.")],
    sigma: Annotated[float, typer.Option("--sigma", help="Sum of weights after each iteration")] = 5,
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory to save .")] = "output-ricci-flow/ricci-graph", 
    tol_curvature: Annotated[float, typer.Option("--tol-curvature", "-t", help="Tolerance for curvature. If all curvatures are smaller than this, then the algorithm stop.")] = 1e-11,
    sequence_graph: Annotated[bool, typer.Option("--sequence-graph", "-s", help="If set, define node distributions as a sequence graph, otherwise use the one for variation graphs (considering paths).")] = False,
    log_level: Annotated[str, typer.Option("--log-level", "-l", help="Log level. Default: INFO/DEBUG")] = "INFO",
    undirected: Annotated[bool, typer.Option("--undirected", "-u", help="Load GFA as undirected graph")] = False,
    eps: Annotated[float, typer.Option("--eps", help="Step size for curvature to update weights")] = 0.1,
    min_weight: Annotated[float, typer.Option("--min-weight", help="minimum weight for edge")] = 0.000001,
    ):

    from pathlib import Path
    from panricci.ricci_flow import NormalizedRicciFlow #NormalizedRicciFlow2 as NormalizedRicciFlow
    from panricci.utils import GFALoader
    
    if sequence_graph:
        from panricci.node_distributions.sequence_graph import DistributionNodes
    else:
        from panricci.node_distributions.variation_graph import DistributionNodes
    
    dirsave    = Path(outdir); dirsave.mkdir(exist_ok=True, parents=True)
    gfa_loader = GFALoader(undirected=undirected)

    # load graph
    G = gfa_loader(gfa)

    # compute distribution of nodes of the graph
    distribution = DistributionNodes(G, alpha=0.5)

    # Create Normalized Ricci Flow instance
    nrf = NormalizedRicciFlow(
        G,
        distribution=distribution,
        log_level=log_level,
        use_rescaling=True,
        eps=eps,
        tol_curvature=tol_curvature,
        dirsave_graphs=outdir,
        lr_schedule="adaptive",
    )

    # Run for specified iterations
    G_final = nrf.run(iterations=iterations, name="test_graph")
    
    print(f"\nFinal edge weights:")
    for u, v, data in G_final.edges(data=True):
        print(f"Edge ({u}, {v}): weight={data['weight']:.6f}, curvature={data['curvature']:.6f}")

@app.command("align", help="Alignment of ricci graphs.")
def align(
    ricci_graph1: Annotated[Path, typer.Option("--ricci-graph1", "-r1", help="Path to the first Ricci graph file.")],
    ricci_graph2: Annotated[Path, typer.Option("--ricci-graph2", "-r2", help="Path to the second Ricci graph file.")],
    path_save: Annotated[Path, typer.Option("--path-save", "-p", help="Path to save the alignment results. Default: 'output-ricci-flow/align/ricci1-ricci2.tsv'.")] = "output-ricci-flow/align/ricci1-ricci2.tsv",
    weight_node_labels: Annotated[float, typer.Option("--weight-node-labels", "-w", min=0, max=0.99, help="Weight for node labels in the alignment cost function. Default: 0. If > 0 then --gfa1 and --gfa2 must be provided")] = 0.0,
    gfa1: Annotated[str, typer.Option("--gfa1", "-g1", help="Path to the first GFA file. Required if metadata-nodes is set.")] = None,
    gfa2: Annotated[str, typer.Option("--gfa2", "-g2", help="Path to the second GFA file. Required if metadata-nodes is set.")] = None, 
    log_level: Annotated[str, typer.Option("--log-level", "-l", help="Log level. Default: INFO.")] = "INFO",
    store_bipartite: Annotated[bool, typer.Option("--store-bipartite", "-sb", help="Store the bipartite graph used for alignment.")] = False,
):
    
    from pathlib import Path
    import networkx as nx

    from panricci.utils import GFALoader
    from panricci.alignment import GraphAlignment, parse_alignment

    dirsave=Path(path_save).parent
    dirsave.mkdir(exist_ok=True, parents=True)

    path_save_bipartite = dirsave.joinpath("bipartite-graph.edgelist") if store_bipartite else None
    aligner = GraphAlignment(
                ricci_embedding = True, 
                weight_node_labels = weight_node_labels,
                path_save_bipartite = path_save_bipartite,
                log_level = log_level,
                ) 

    # load pangenome graphs with node info, and add the weight
    g1 = nx.read_edgelist(ricci_graph1, data=True, create_using=nx.DiGraph) # DiGraph to find source and sinks nodes   
    g2 = nx.read_edgelist(ricci_graph2, data=True, create_using=nx.DiGraph) # Otherwise, provide source and sinks nodes (TODO)

    # add weights to graphs
    if gfa1 and gfa2:
        gfa_loader = GFALoader(undirected=False)
        graph1 = gfa_loader(gfa1)
        for edge, d in g1.edges.items():
            graph1.edges[edge]["weight"] = d["weight"]
            graph1.edges[edge]["curvature"] = d["curvature"]

        graph2 = gfa_loader(gfa2)
        for edge, d in g2.edges.items():
            graph2.edges[edge]["weight"] = d["weight"]
            graph2.edges[edge]["curvature"] = d["curvature"]
    else:
        graph1 = g1 
        graph2 = g2

    alignment = aligner(graph1, graph2, name="alignment") 

    parse_alignment(alignment, graph1, graph2).\
    sort_values(by="cost_alignment").\
    to_csv(path_save, sep="\t", float_format="%.16f")

@app.command("modify-gfa", help="Modify a GFA file by removing edges, nodes or paths.")
def modify_gfa(
    gfa: Annotated[Path, typer.Option("--gfa", "-g", help="Path to the GFA file.")],
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Output directory to save the modified GFA file.")] = "output-ricci-flow/modified-gfa",
    remove_edges: Annotated[str, typer.Option("--remove-edges", help="Comma-separated list of edges to remove (format: node1-node2). Nodes are collapsed by default.")] = None,
    remove_nodes: Annotated[str, typer.Option("--remove-nodes", help="Comma-separated list of nodes to remove. Predecessors and successors are reconnected.")] = None,
    remove_paths: Annotated[str, typer.Option("--remove-paths", help="Comma-separated list of paths to remove.")] = None,
    no_collapse: Annotated[bool, typer.Option("--no-collapse", help="When removing edges, do not collapse adjacent nodes.")] = False,
    no_reconnect: Annotated[bool, typer.Option("--no-reconnect", help="When removing nodes, do not reconnect predecessors to successors.")] = False,
    remove_unused_nodes: Annotated[bool, typer.Option("--remove-unused-nodes", help="When removing paths, also remove nodes not used by any remaining path.")] = False,
):
    from pathlib import Path
    from panricci.utils import GFAGraph

    dirsave = Path(outdir); dirsave.mkdir(exist_ok=True, parents=True)

    # load graph
    G = GFAGraph(gfa)
    print(f"Loaded: {G}")

    # Remove specified edges
    if remove_edges:
        for edge in remove_edges.split(","):
            node1, node2 = edge.strip().split("-")
            G.remove_edge(node1.strip(), node2.strip(), collapse_nodes=not no_collapse)

    # Remove specified nodes
    if remove_nodes:
        for node in remove_nodes.split(","):
            G.remove_node(node.strip(), reconnect=not no_reconnect)

    # Remove specified paths
    if remove_paths:
        for path_name in remove_paths.split(","):
            G.remove_path(path_name.strip(), remove_unused_nodes=remove_unused_nodes)

    # Save modified graph as GFA
    modified_gfa_path = dirsave.joinpath(f"{Path(gfa).stem}_modified.gfa")
    G.save(str(modified_gfa_path))
    print(f"\nSaved modified GFA to: {modified_gfa_path}")
    print(f"Result: {G}")


if __name__ == "__main__":
    app()