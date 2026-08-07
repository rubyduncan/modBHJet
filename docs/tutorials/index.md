# Tutorials

The tutorials are written as lightweight documentation pages and are paired with
the notebooks in `examples/`. The Read the Docs build does not execute heavy
model notebooks by default.

````{grid} 1 1 2 2
:gutter: 2

```{grid-item-card} BLJet dynamics
Inspect the jet speed, radius, density, magnetic field, and nonthermal fraction.
```

```{grid-item-card} Full spectrum
Compute a multi-zone BHJet spectrum on an observed energy grid.
```

```{grid-item-card} Cutoff comparison
Compare exponential, super-exponential, and sech-squared particle cutoffs.
```

```{grid-item-card} Gammapy wrapper
Use the `BHJetSpectralModel` wrapper for data/model comparison workflows.
```
````

## Existing notebooks

- `tutorials/welcome.ipynb` — Binder/JupyterHub landing page
- `tutorials/1_pilot_explorer.ipynb` — widget-based parameter explorer
- `tutorials/2_source_SED_analysis.ipynb` — one-zone to multi-zone SED exercise
- `examples/python_classes/bljet_dynamics_example.ipynb`
- `examples/python_classes/single_zone_example.ipynb`
- `examples/python_classes/bhjet_example.ipynb`
- `examples/gammapy/gammapy_intro.ipynb` — advanced development example
- `examples/gammapy/gammapy_bhjet_comparison.ipynb` — advanced development example
