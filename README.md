<!-- ![pantera](images/PANRICCI-removebg-preview.png) -->
<img src="images/PANRICCI-removebg-preview.png" width="200" height="200">


# Alignment of Pangenome Graphs with Ricci Flow

## Installation
```bash
pip install git+https://github.com/jorgeavilacartes/panricci.git
```

or (recommended for developer) create a conda environment with the library in edition mode
```bash
git clone git@github.com:jorgeavilacartes/panricci.git
cd panricci
conda env create -f panricci.yml
conda activate panricci

panricci --help
```

`panricci` is a library that deals with pangenome graphs (variation graphs and sequence graphs)
with tools from Riemannian Geometry. 

A pangenome graph in the `panricci` universe is a manifold $\mathcal{X}$ provided of a metric $d$ 
(weights of the edges), where nodes encode information in probability distributions over its (1-hop) neighborhood. 

This manifold is evolved over time by using the **Ricci-Flow**, an algorithm that leverages the notion
of curvature of edges to modify its weights until the curvature is constant (equal to 0).

Once the graph reach the state of constant curvature, we can use it to perform alignment of two graphs
by definining coordinates with respect to source and sink nodes in our (directed) pangenome graphs.

**NOTE** Each pangenome graph (input file, and manifold) is assumed to be one single connected component.
___
## CLI
`panricci` works with pangenome graphs (variation graphs or sequence graphs) in `.gfa` format.

```bash
$ panricci --help
```

Available commands:

| Command | Description |
|---------|-------------|
| `ricci-flow` | Apply Ricci Flow to a graph |
| `ricci-flow2` | Apply Ricci Flow (v2) with configurable step size and minimum weight |
| `normalized-ricci-flow` | Apply Normalized Ricci Flow with weight rescaling |
| `align` | Alignment of ricci graphs |
| `modify-gfa` | Modify a GFA file by removing edges, nodes or paths |

---

### `ricci-flow`

Apply the discrete Ricci-Flow algorithm to evolve the metric of a pangenome graph.

```bash
$ panricci ricci-flow \
	--gfa data/test5.gfa \
	--iterations 1000 \
	--tol-curvature 1e-15 \
	--outdir output/test5
```

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--gfa` | `-g` | Path to the GFA file | *required* |
| `--iterations` | `-i` | Maximum number of iterations | *required* |
| `--outdir` | `-o` | Output directory | `output-ricci-flow/ricci-graph` |
| `--tol-curvature` | `-t` | Tolerance for curvature convergence | `1e-11` |
| `--sequence-graph` | `-s` | Use sequence graph node distributions (instead of variation graph) | `False` |
| `--log-level` | `-l` | Log level | `INFO` |
| `--save-intermediate-graphs` | `-si` | Save graphs at each iteration | `False` |
| `--undirected` | `-u` | Load GFA as undirected graph | `False` |

- The algorithm runs for at most `--iterations` iterations, or until all curvatures are smaller than `--tol-curvature` (the goal is to stop when all curvatures are 0).
- The `--undirected` option allows the computation of the Wasserstein distance over the local subgraph avoiding directions, which is needed to keep the Wasserstein distance as a metric.
- By default node distributions are defined for a variation graph (paths are considered). If you have a sequence graph where paths are not part of the `.gfa` file, use `--sequence-graph`.

---

### `ricci-flow2`

An alternative Ricci Flow implementation with configurable step size and minimum edge weight.

```bash
$ panricci ricci-flow2 \
	--gfa data/test5.gfa \
	--iterations 1000 \
	--eps 0.1 \
	--min-weight 0.000001 \
	--outdir output/test5
```

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--gfa` | `-g` | Path to the GFA file | *required* |
| `--iterations` | `-i` | Maximum number of iterations | *required* |
| `--outdir` | `-o` | Output directory | `output-ricci-flow/ricci-graph` |
| `--tol-curvature` | `-t` | Tolerance for curvature convergence | `1e-11` |
| `--sequence-graph` | `-s` | Use sequence graph node distributions | `False` |
| `--log-level` | `-l` | Log level (`INFO`/`DEBUG`) | `INFO` |
| `--undirected` | `-u` | Load GFA as undirected graph | `False` |
| `--eps` | | Step size for curvature to update weights | `0.1` |
| `--min-weight` | | Minimum weight for an edge (prevents weights going to zero) | `0.000001` |

---

### `normalized-ricci-flow`

Apply Normalized Ricci Flow, which rescales weights after each iteration to maintain a constant total weight.

```bash
$ panricci normalized-ricci-flow \
	--gfa data/test5.gfa \
	--iterations 1000 \
	--sigma 5 \
	--eps 0.1 \
	--outdir output/test5
```

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--gfa` | `-g` | Path to the GFA file | *required* |
| `--iterations` | `-i` | Maximum number of iterations | *required* |
| `--sigma` | | Target sum of weights after each iteration | `5` |
| `--outdir` | `-o` | Output directory | `output-ricci-flow/ricci-graph` |
| `--tol-curvature` | `-t` | Tolerance for curvature convergence | `1e-11` |
| `--sequence-graph` | `-s` | Use sequence graph node distributions | `False` |
| `--log-level` | `-l` | Log level (`INFO`/`DEBUG`) | `INFO` |
| `--undirected` | `-u` | Load GFA as undirected graph | `False` |
| `--eps` | | Step size for curvature to update weights | `0.1` |
| `--min-weight` | | Minimum weight for an edge | `0.000001` |

---

### `align`

Align two ricci graphs by finding a mapping between their nodes based on relative node representations.

```bash
$ panricci align \
	--ricci-graph1 output/test5/test5-ricciflow-52.edgelist \
	--ricci-graph2 output/test5/test5-ricciflow-52.edgelist \
	--path-save output/alignment.tsv
```

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--ricci-graph1` | `-r1` | Path to the first Ricci graph file | *required* |
| `--ricci-graph2` | `-r2` | Path to the second Ricci graph file | *required* |
| `--path-save` | `-p` | Path to save alignment results | `output-ricci-flow/align/ricci1-ricci2.tsv` |
| `--weight-node-labels` | `-w` | Weight for node labels in the cost function (0-0.99). If > 0, `--gfa1` and `--gfa2` must be provided | `0.0` |
| `--gfa1` | `-g1` | Path to the first GFA file (for node metadata) | `None` |
| `--gfa2` | `-g2` | Path to the second GFA file (for node metadata) | `None` |
| `--log-level` | `-l` | Log level | `INFO` |
| `--store-bipartite` | `-sb` | Store the bipartite graph used for alignment | `False` |

The output is a TSV file with the following columns:

- `node1`: the identifier of the node in ricci graph 1
- `node2`: the identifier of the node in ricci graph 2
- `edge`: tuples of strings of the form `<node_id>-<graph_id>`
- `cost_alignment`: the cost of aligning `node1` and `node2` (euclidean distance between their relative node representations)

Example output:

```
 	 edge          	 cost_alignment	 node1	 node2
0	 ['4-1', '4-2']	            0.0	     4	     4
1	 ['3-1', '3-2']	            0.0	     3	     3
2	 ['2-1', '2-2']	            0.0	     2	     2
3	 ['1-1', '1-2']	            0.0	     1	     1
```

To include node metadata (label and node depth), provide the original `.gfa` files:

```bash
$ panricci align \
	--ricci-graph1 output/test5/test5-ricciflow-52.edgelist \
	--ricci-graph2 output/test5/test5-ricciflow-52.edgelist \
	--path-save output/alignment.tsv \
	--gfa1 data/test5.gfa \
	--gfa2 data/test5.gfa
```

To penalize alignment of nodes with different labels, use `--weight-node-labels`:

```bash
$ panricci align \
	--ricci-graph1 output/test5/test5-ricciflow-52.edgelist \
	--ricci-graph2 output/test5/test5-ricciflow-52.edgelist \
	--path-save output/alignment.tsv \
	--weight-node-labels 0.5 \
	--gfa1 data/test5.gfa \
	--gfa2 data/test5.gfa
```

---

### `modify-gfa`

Modify a GFA file by removing edges, nodes, or paths. When removing edges, the adjacent nodes are collapsed into a single node with concatenated sequences by default. When removing nodes, predecessors and successors are reconnected by default.

```bash
$ panricci modify-gfa \
	--gfa data/test5.gfa \
	--outdir output/modified \
	--remove-nodes "6,7"
```

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--gfa` | `-g` | Path to the GFA file | *required* |
| `--outdir` | `-o` | Output directory for the modified GFA file | `output-ricci-flow/modified-gfa` |
| `--remove-edges` | | Comma-separated edges to remove (format: `node1-node2`) | `None` |
| `--remove-nodes` | | Comma-separated nodes to remove | `None` |
| `--remove-paths` | | Comma-separated paths to remove | `None` |
| `--no-collapse` | | Do not collapse adjacent nodes when removing edges | `False` |
| `--no-reconnect` | | Do not reconnect predecessors to successors when removing nodes | `False` |
| `--remove-unused-nodes` | | Remove nodes not used by any remaining path after path removal | `False` |

#### Examples

Remove edges (collapses nodes by default):
```bash
$ panricci modify-gfa --gfa input.gfa --remove-edges "4-5,10-11"
```

Remove edges without collapsing:
```bash
$ panricci modify-gfa --gfa input.gfa --remove-edges "4-5" --no-collapse
```

Remove nodes (reconnects predecessors to successors):
```bash
$ panricci modify-gfa --gfa input.gfa --remove-nodes "6,7,8"
```

Remove paths and clean up unused nodes:
```bash
$ panricci modify-gfa --gfa input.gfa --remove-paths "seq1,seq2" --remove-unused-nodes
```

The modified GFA file is saved as `<original_name>_modified.gfa` in the output directory.


## Python API

#### 1. Create ricci-graphs
```python
from panricci.ricci_flow import RicciFlow
from panricci.utils import GFALoader
from panricci.node_distributions.variation_graph import DistributionNodes

# load variation graph
gfa_loader = GFALoader(undirected=False)
G = gfa_loader(gfa)

# compute distribution for each node of the graph
distribution_nodes = DistributionNodes(G, alpha=0.5)

# Initialize Ricci-Flow
ricci_flow = RicciFlow(G, 
	distribution=distribution_nodes, # the distribution over the 1-hop neighborhood for each node in the graph
	save_last=False,                 # will overwrite the results in each iteration to keep the last one: outfile will be "{name}-ricciflow-{it}.edgelist"
	save_intermediate_graphs=True,   # will save results for all iteration: outfile will be "{name}-ricciflow.edgelist"
	dirsave_graphs="outdir",         # directory to save results
	tol_curvature=1e-15,             # tolerance of minimum curvature to stop Ricci-Flow
	overwrite=False,                 # If the dirsave_graphs directory exists, it will raise an Exception
	)

# apply Ricci-Flow
G_ricci = ricci_flow.run(
	iterations=1000,       # maximum number of iterations to run Ricci-Flow  
	name="id-graph"        # some identifier to store results
	)
``` 

#### 2. Align ricci-graphs

```python
from pathlib import Path
import networkx as nx

from panricci.utils import GFALoader
from panricci.alignment import GraphAlignment, parse_alignment

dirsave=Path(path_save)
dirsave.parent.mkdir(exist_ok=True, parents=True)

# initialize aligner
aligner = GraphAlignment(
	ricci_embedding = True, 
	weight_node_labels = 0.5,    # set to 0 if you only want to consider embedding cost
	path_save_bipartite = "bipartite-graph.edgelist", # store bipartite graph used for alignment
	) 

# load ricci graphs
g1 = nx.read_edgelist(ricci_graph1, data=True, create_using=nx.DiGraph)    
g2 = nx.read_edgelist(ricci_graph2, data=True, create_using=nx.DiGraph)

gfa_loader = GFALoader(undirected=False)
graph1 = gfa_loader(gfa1)
for edge, d in g1.edges.items():
	graph1.edges[edge]["weight"] = d["weight"]
	graph1.edges[edge]["curvature"] = d["curvature"]

graph2 = gfa_loader(gfa2)
for edge, d in g2.edges.items():
	graph2.edges[edge]["weight"] = d["weight"]
	graph2.edges[edge]["curvature"] = d["curvature"]

alignment = aligner(g1, g2, name="alignment") 
df_alignment = parse_alignment(alignment, graph1, graph2)
df_aligment.to_csv("alignment.tsv", sep="\t")
```