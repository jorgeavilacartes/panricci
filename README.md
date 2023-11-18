<!-- ![pantera](images/PANRICCI-removebg-preview.png) -->
<img src="images/PANRICCI-removebg-preview.png" width="200" height="200">

# Clustering Pangenome Graphs with Ricci Flow

## Create environment
```bash
mamba env create -n panricci -f envs/panricci.yml
```

Activate environment
```bash
mamba activate panricci
```

Install panricci as library in the environment
From the main folder (parent folder of panricci/) run the following
[check](https://goodresearch.dev/setup#pip-install-your-package)
```bash
pip install -e .
```
**NOTE** If changes are not recognized, reinstall library with the pevious command.


```python
from panricci import  (
    RicciFlow,
    Index,
    PanRicciSimilarity
)
```


# TODO
- [ ] Plot distributions of nodes after each iteration of ricci flow