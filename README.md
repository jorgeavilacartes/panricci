<!-- ![pantera](images/PANRICCI-removebg-preview.png) -->
<img src="images/PANRICCI-removebg-preview.png" width="200" height="200">

# Alignment of Pangenome Graphs with Ricci Flow

`panricci` is a library that deals with pangenome graphs (variation graphs and sequence graphs)
with tools from Riemannian Geometry. 

A pangenome graph in the `panricci` universe is a manifold $\mathcal{X}$ provided of a metric $d$ 
(weights of the edges), where nodes encode information in probability distributions. 

This manifold is evolved over time by using the Ricci-Flow, an algorithm that leverages the notion
of curvature of edges to modify its weights until the curvature is constant.

Once the graph reach the state of constant curvature, we can use it to perform alignment of two graphs.

## 1. Install
Clone this repository, then from the main folder, follow the next steps

### Create environment
(use `conda`/`miniconda/``mamba`/`micromamba`)
```bash
mamba env create -n panricci -f envs/panricci.yml
```

Activate environment
```bash
mamba activate panricci
```

Install `panricci` as library in the environment.

<!-- [check](https://goodresearch.dev/setup#pip-install-your-package) -->
```bash
pip install -e .
```
**NOTE** If changes are not recognized, reinstall library with the pevious command.

## 2. PanRicci

To see examples check
1. `notebooks/ricci-flow.ipynb`
2. `notebooks/graph-alignment.ipynb`

# TODO
- [ ] Plot distributions of nodes after each iteration of ricci flow