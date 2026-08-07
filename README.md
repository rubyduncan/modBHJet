
## Link to the online notebook: 

[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rubyduncan/modBHJet/bhjet_tutorial?urlpath=lab/tree/tutorials/00_welcome.ipynb)

The public tutorial includes an overview of what BHJet calculates with an interactive
notebook to change parameters, and a test of fitting a BHJet spectral energy distribution to M87 data. Binder sessions are temporary and public, anything that you produce and want to save needs to be downloaded. 

## Installation

If you'd like to install the code locally on your computer instead of using Binder: BHJet has a compiled C++ extension and uses the Kariba radiation library (this is a separate code that works in tandem with BHJet and holds the radiation defintions). The [`environment.yml`](environment.yml) will take care of package dependencies. 

### 1. Clone BHJet and Kariba

The default build expects Kariba next to the BHJet directory. The commands below
reproduce the current Binder tutorial locally (from the branch bhjet_tutorial), including this environment file and the tutorial notebooks.

```bash
git clone --branch bhjet_tutorial https://github.com/rubyduncan/modBHJet.git
cd modBHJet
git clone https://github.com/rubyduncan/kariba.git ../kariba
git -C ../kariba checkout 6ca6968adf8568e99b2010212fd1419a79d94d35
```

The final command pins Kariba to the version used by the Binder tutorial, so
the local and public environments use the same radiation backend.

### 2. Create the environment

Using conda or micromamba, create the environment: 

```bash
conda env create -f environment.yml
conda activate bhjet
```

### 3. Build and verify BHJet

```bash
python -m pip install --no-build-isolation -v .
python tutorials/smoke_test.py
```

The smoke test imports the compiled extension, loads the bundled example data,
and calculates a small spectrum. Start JupyterLab with:

```bash
jupyter lab
```

### Kariba in a different directory

If you already have Kariba elsewhere, keep it there and modify the kariba_SOURCE_DIR path in the [CMake file](CMakeLists.txt) so that it can be found when installing BHJet:

To rebuild after changing C++ or Python code, repeat the installation command.
