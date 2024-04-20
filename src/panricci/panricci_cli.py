import typer

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

# types for typer
from rich.progress import track
from rich import print 
from rich.console import Console
from rich.markdown import Markdown

console = Console()
app = typer.Typer(name="PanRicci", rich_markup_mode="rich",
                  help="""
                  :cat: Welcome to [bold]Pan[/bold]Ricci Alignment of pangenome graphs with Ricci-Flow
                  """
                  )


@app.command("docs", help="Open documentation webpage.")
def github() -> None:
    typer.launch("https://github.com/jorgeavilacartes/panricci")


@app.command("ricci-flow", help="apply ricci flow to a graph.")
def ricci_flow(
        gfa: str, 
        iterations: int,
        undirected: bool = False,
        sequence_graph: bool = False,
        dirsave: str = "output-ricci-flow", 
    ):
    from pathlib import Path
    from panricci import RicciFlow
    from panricci.utils import GFALoader
    
    if sequence_graph:
        from panricci.distributions.sequence_graph import DistributionNodes
    else:
        from panricci.distributions.variation_graph import DistributionNodes
    
    dirsave = Path(dirsave)
    dirsave.mkdir(exist_ok=True, parents=True)

    # load graph
    gfa_loader = GFALoader(undirected=undirected)
    G = gfa_loader(gfa)

    # compute distribution of nodes
    distribution = DistributionNodes(G, alpha=0.5)

    # Initialize ricci-flow
    ricci_flow = RicciFlow(G, 
                               distribution, 
                               save_last=False, 
                               save_intermediate_graphs=True, 
                               dirsave_graphs=dirsave,
                               )
    # ricci_flow = RicciFlow(G, distribution, dirsave_graphs=dirsave, save_last=False, save_intermediate_graphs=True)
    G_ricci = ricci_flow.run(iterations=iterations, name=Path(gfa).stem)

    return G_ricci

@app.command("align", help="Alignment of ricci graphs.")
def alignment(
        gfa1: str,
        gfa2: str,
        iterations: int,
        dirsave: str = "output-ricci-flow/alignment",
):
    
    from pathlib import Path; dirsave=Path(dirsave)
    dirsave.mkdir(exist_ok=True, parents=True)
    from panricci.alignment import GraphAlignment, parse_alignment

    g1 = ricci_flow(gfa1, iterations)   
    g2 = ricci_flow(gfa2, iterations)

    aligner = GraphAlignment(
    ricci_embedding = True, 
    seq_embedding = False, ) # kmer_size=4)
    alignment = aligner(g1, g2, name="alignment") 

    parse_alignment(alignment, g1, g2).\
    sort_values(by="cost_alignment").\
    to_csv(dirsave.joinpath(f"alignment-{Path(gfa1).stem}-{Path(gfa2).stem}.tsv"),sep="\t")

if __name__ == "__main__":
    app()