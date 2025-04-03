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

1. You need to evolve the metric of the graph, this is randomnly initialized.
2. Once you have two graphs with the metrics after Ricci-Flow, we can align them.
3. Check the results. 

```{bash}
$ panricci --help 
                                                                                              
 Usage: panricci [OPTIONS] COMMAND [ARGS]...                                                  
                                                                                              
 🐱 Welcome to PanRicci :  Alignment of pangenome graphs with Ricci-Flow                      
                                                                                              
╭─ Options ──────────────────────────────────────────────────────────────────────────────────╮
│ --install-completion          Install completion for the current shell.                    │
│ --show-completion             Show completion for the current shell, to copy it or         │
│                               customize the installation.                                  │
│ --help                        Show this message and exit.                                  │
╰────────────────────────────────────────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────────────────────────────────────────╮
│ ricci-flow   apply ricci flow to a graph.                                                  │
│ align        Alignment of ricci graphs.                                                    │
╰────────────────────────────────────────────────────────────────────────────────────────────╯
```

### Example 
We will use one graph from the `data` folder

#### 1. Create ricci-graphs
```bash
$ panricci ricci-flow --gfa data/test5.gfa --iterations 1000 --tol-curvature 1e-15 --undirected  --outdir output/test5 
``` 

- this tells `panricci` to apply the discrete Ricci-Flow algorithm for at most 1000 iterations, or until all curvatures are smaller than 1e-15 (our goal is to stop when all curvatures are 0)
- the results for each iteration (weights and curvature for each edge of the graph) will be saved in the `output/test5` directory
- the option `--undirected` allows the computation of the Wasserstein distance over the local subgraph of two nodes avoiding the directions, which is needed to keep the Wasserstein distance as a metric, so we recommend to keep this option.
- by default node distributions are defined for a variation graph, i.e. paths are considered. If you have a sequence graph, where paths are not part of the .gfa file, then you can add the option `--sequence-graph` 

#### 2. Align ricci-graphs
The **alignment** consists on finding a mapping **between nodes of two graphs**.

You can either create the ricci graph for another gfa in the folder, or we can use the same ricci graph to be align against itself, which we will do now.

Notice that in the `output/test5` folder you will find a lot (one for each iteration), we will take the last one, which should be the file `test5-ricciflow-52.edgelist`

```bash
$ panricci align \
	--ricci-graph1 output/test5/test5-ricciflow-52.edgelist \
	--ricci-graph2 output/test5/test5-ricciflow-52.edgelist \
	--path-save output/alignment.tsv
```

you should get a csv file with the following information inside

```bash
 	 edge          	 cost_alignment	 node1	 node2
0	 ['4-1', '4-2']	            0.0	     4	     4	       	       	            	 
1	 ['3-1', '3-2']	            0.0	     3	     3	       	       	            	 
2	 ['2-1', '2-2']	            0.0	     2	     2	       	       	            	 
3	 ['1-1', '1-2']	            0.0	     1	     1	       	       	            	 
```

the columns,

- `node1`: the identifier of the node in the ricci graph 1
- `node2`: the identifier of the node in the ricci graph 2
- `edge` : are tuples of strings of the form '<node_id>-<graph_id>' 
- `cost_alignment`: the cost of aligning `node1` and `node2` (euclidean distance between their relative node representations)


as we can see from this result, the first row says that node 4 of the graph1 was aligned to node 4 of graph2, with a cost alignment of 0, which means that they have the same relative node representation.


Additionally, you can ask to add node metadata (label and node depth) to the output by providing the original .gfa files as follows

```bash
$ panricci align \
	--ricci-graph1 output/test5/test5-ricciflow-52.edgelist \
	--ricci-graph2 output/test5/test5-ricciflow-52.edgelist \
	--path-save output/alignment.tsv \
	--node-metadata \
	--gfa1  data/test5.gfa \
	--gfa2 data/test5.gfa
```

The output should look like this

```bash
 	 edge          	 cost_alignment	 node1	 node2	 label1      	 label2      	 node_depth1	 node_depth2
0	 ['4-1', '4-2']	            0.0	     4	     4	 CAACGTTTTTTT	 CAACGTTTTTTT	         0.5	         0.5
1	 ['3-1', '3-2']	            0.0	     3	     3	 CTAGA       	 CTAGA       	         1.0	         1.0
2	 ['2-1', '2-2']	            0.0	     2	     2	 A           	 A           	         0.5	         0.5
3	 ['1-1', '1-2']	            0.0	     1	     1	 C           	 C           	         1.0	         1.0
```