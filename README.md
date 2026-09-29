# Ae4 salivary secretion

Code, numerical data and supplementary material for **Intracellular chloride availability limits Slc4a9-dependent salivary secretion**.

## Start here

**Run the manuscript model**

```bash
python supplementary_material/model/run_model.py
```

**Model implementation**  
[`supplementary_material/model/ae4_salivary_model.py`](supplementary_material/model/ae4_salivary_model.py)

**Supplementary material**  
[`supplementary_material/`](supplementary_material/README.md)

The supplementary directory contains the model used for the principal control and Ae4 knockout simulations, its preserved transport dependencies, full-precision parameter and onset records, stored simulation data, and the Python plotting scripts.

The principal public model entry point is `PairedModel` in `supplementary_material/model/ae4_salivary_model.py`. The executable entry point is `supplementary_material/model/run_model.py`.

## Repository structure

- `supplementary_material/` manuscript model, parameter records, numerical data and plotting scripts
- `manuscript/` manuscript source and figures
- `analysis/` numbered scientific analyses and development history
- `archive/` earlier material retained for provenance
- `docs/` project ledgers and supporting documentation

The historical analysis directories are retained to document the development of the model. For the simulations reported in the current manuscript, use the model identified above.
