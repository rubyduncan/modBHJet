# Installation

The recommended user workflow is to install `bhjet` into a clean Python
environment. The package builds a compiled extension module and links against
GSL and the local Kariba backend.

## Environment

The root `environment.yml` provides the complete
supported local environment, including the C++ compiler toolchain, GSL, OpenMP,
Nanobind, and optional notebook/Gammapy dependencies:

```bash
micromamba env create -f environment.yml
micromamba activate bhjet
```

`conda env create -f environment.yml` and `mamba env create -f environment.yml`
are equivalent alternatives.

## Kariba backend

Clone the Kariba source beside the modBHJet checkout and pin it to the tutorial
revision:

```bash
git clone https://github.com/rubyduncan/kariba.git ../kariba
git -C ../kariba checkout 6ca6968adf8568e99b2010212fd1419a79d94d35
```

## Install from a local checkout

From the `modBHJet` directory:

```bash
python -m pip install --no-build-isolation -v .
```

This invokes `scikit-build-core`, configures CMake, compiles the C++ extension,
and installs the Python package into the active environment. If Kariba lives
elsewhere, pass `CMAKE_ARGS="-Dkariba_SOURCE_DIR=/absolute/path/to/kariba"` to
the install command.

## Verify the install

```python
import bhjet
from bhjet import BLJet, BHJet

print(bhjet.__file__)
print("cutoff_type" in BLJet.__init__.__doc__)
```

If `bhjet.__file__` points to a stale `site-packages` install, reinstall or use
the local development workflow.
