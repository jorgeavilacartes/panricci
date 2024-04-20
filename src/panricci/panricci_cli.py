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
# def ricci_flow(path_graph: str,
               
               
#                ): 


if __name__ == "__main__":
    app()