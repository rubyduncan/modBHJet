# modBHJet Documentation

modBHJet is a modular black-hole jet model for computing steady-state jet spectra.
It combines a high-level Python interface, C++ jet/radiation-zone classes, and the
Kariba radiation backend.

This site is organized around how users actually work with the model: install it,
choose a jet dynamics model, compute particles and radiation, inspect outputs, and
then connect those outputs to data-analysis workflows.

```{toctree}
:maxdepth: 2
:caption: Getting started

getting-started/installation
getting-started/local-development
getting-started/import-checks
```

```{toctree}
:maxdepth: 2
:caption: Model and physics

model/architecture
model/data-flow
physics/notes
```

```{toctree}
:maxdepth: 2
:caption: Tutorials and demos

tutorials/index
tutorials/bljet-dynamics
tutorials/full-spectrum
tutorials/cutoff-comparison
tutorials/gammapy
tutorials/interactive-demo
deployment/tutorial-hosting
```

```{toctree}
:maxdepth: 2
:caption: API reference

api/python
api/cpp
```

```{toctree}
:maxdepth: 2
:caption: Kariba backend

kariba/index
kariba/about
kariba/installation
kariba/using
```
