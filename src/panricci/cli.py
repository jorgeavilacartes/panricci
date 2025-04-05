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
    undirected: Annotated[bool, typer.Option("--undirected", "-u", help="Treat the graph as undirected for Wasserstein distance computation.")] = False,
    sequence_graph: Annotated[bool, typer.Option("--sequence-graph", "-s", help="If set, define node distributions as a sequence graph, otherwise use the one for variation graphs (considering paths).")] = False,
    log_level: Annotated[str, typer.Option("--log-level", "-l", help="Log level. Default: INFO.")] = "INFO",
    save_intermediate_graphs: Annotated[bool, typer.Option("--save-intermediate-graphs", "-si", help="Save intermediate graphs.")] = False,
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

@app.command("align", help="Alignment of ricci graphs.")
def align(
    ricci_graph1: Annotated[Path, typer.Option("--ricci-graph1", "-r1", help="Path to the first Ricci graph file.")],
    ricci_graph2: Annotated[Path, typer.Option("--ricci-graph2", "-r2", help="Path to the second Ricci graph file.")],
    path_save: Annotated[Path, typer.Option("--path-save", "-p", help="Path to save the alignment results. Default: 'output-ricci-flow/align/ricci1-ricci2.tsv'.")] = "output-ricci-flow/align/ricci1-ricci2.tsv",
    metadata_nodes: Annotated[bool, typer.Option("--metadata-nodes", "-m", help="Include metadata for nodes in the alignment, in which case --gfa1 and --gfa2 must be provided")] = False,
    gfa1: Annotated[str, typer.Option("--gfa1", "-g1", help="Path to the first GFA file. Required if metadata-nodes is set.")] = None,
    gfa2: Annotated[str, typer.Option("--gfa2", "-g2", help="Path to the second GFA file. Required if metadata-nodes is set.")] = None, 
    weight_node_labels: Annotated[float, typer.Option("--weight-node-labels", "-w", min=0, max=0.99, help="Weight for node labels in the alignment cost function. Default: 0.0.")] = 0.0,
    log_level: Annotated[str, typer.Option("--log-level", "-l", help="Log level. Default: INFO.")] = "INFO",
    store_bipartite: Annotated[bool, typer.Option("--store-bipartite", "-sb", help="Store the bipartite graph used for alignment.")] = False,
):
    
    from pathlib import Path
    import networkx as nx

    from panricci.utils import GFALoader
    from panricci.alignment import GraphAlignment, parse_alignment

    dirsave=Path(path_save).parent
    dirsave.parent.mkdir(exist_ok=True, parents=True)

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
    if metadata_nodes:
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
    to_csv(path_save,sep="\t")


if __name__ == "__main__":
    app()