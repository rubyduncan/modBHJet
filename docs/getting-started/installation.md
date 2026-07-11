# Installation

The recommended user workflow is to install `bhjet` into a clean Python
environment. The package builds a compiled extension module and links against
GSL and the local Kariba backend.

## Environment

For a micromamba environment:

```bash
micromamba create -n bhjet_env python=3.12 compilers cmake make gsl nanobind litgen \
  pandas numpy scipy matplotlib jupyterlab ipywidgets ipympl
micromamba activate bhjet_env
```

For notebook and tutorial work, also install:

```bash
python -m pip install -r docs/requirements.txt
```

## Install from a local checkout

From the `modBHJet` directory:

```bash
python -m pip install .
```

This invokes `scikit-build-core`, configures CMake, compiles the C++ extension,
and installs the Python package into the active environment.

## Verify the install

```python
import bhjet
from bhjet import BLJet, BHJet

print(bhjet.__file__)
print("cutoff_type" in BLJet.__init__.__doc__)
```

If `bhjet.__file__` points to a stale `site-packages` install, reinstall or use
the local development workflow.
