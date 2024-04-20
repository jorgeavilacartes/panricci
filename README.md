<!-- ![pantera](images/PANRICCI-removebg-preview.png) -->
<img src="images/PANRICCI-removebg-preview.png" width="200" height="200">

# Alignment of Pangenome Graphs with Ricci Flow

`panricci` is a library that deals with pangenome graphs (variation graphs and sequence graphs)
with tools from Riemannian Geometry. 

A pangenome graph in the `panricci` universe is a manifold $\mathcal{X}$ provided of a metric $d$ 
(weights of the edges), where nodes encode information in probability distributions over its (1-hop) neighborhood. 

This manifold is evolved over time by using the Ricci-Flow, an algorithm that leverages the notion
of curvature of edges to modify its weights until the curvature is constant (equal to 0).

Once the graph reach the state of constant curvature, we can use it to perform alignment of two graphs
by definining coordinates with respect to source and sink nodes in our (directed) pangenome graphs.

**NOTE** Each pangenome graph (input file, and manifold) is assumed to be one single connected component.  

## 1. Install
Clone this repository, then from the main folder, follow the next steps

### Create environment
(use `conda`/`miniconda/``mamba`/`micromamba`)
```bash
mamba env create -n panricci -f panricci.yml
```

Activate environment
```bash
mamba activate panricci
```

Developer note: This will install `panricci` as a pip library in editable mode.
**NOTE** If changes are not recognized, reinstall library with the pevious command.

___
## CLI
`panricci` works with pangenome graphs (variation graphs or sequence graphs) in `.gfa` format.

1. You need to evolve the metric of the graph, this is randomnly initialized.
2. Once you have two graphs with the metrics after Ricci-Flow, we can align them.
3. Check the results. 

```{bash}
panricci --help
                                                                                           
 Usage: panricci [OPTIONS] COMMAND [ARGS]...                                               
                                                                                           
 🐱 Welcome to PanRicci Alignment of pangenome graphs with Ricci-Flow                      
                                                                                           
╭─ Options ───────────────────────────────────────────────────────────────────────────────╮
│ --install-completion          Install completion for the current shell.                 │
│ --show-completion             Show completion for the current shell, to copy it or      │
│                               customize the installation.                               │
│ --help                        Show this message and exit.                               │
╰─────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ──────────────────────────────────────────────────────────────────────────────╮
│ align        Alignment of ricci graphs.                                                 │
│ docs         Open documentation webpage.                                                │
│ ricci-flow   apply ricci flow to a graph.                                               │
╰─────────────────────────────────────────────────────────────────────────────────────────╯


```

___ 
## PanRicci examples

To see examples check
TODO: add links 
1. Compute Ricci-Flow on a pangenome graph `notebooks/ricci-flow.ipynb`
2. Alignment of pangenome graphs with Ricci-Flow `notebooks/graph-alignment.ipynb`

___ 
# TODO
- [ ] Plot distributions of nodes after each iteration of ricci flow