# Local Development Build

For active C++ development, use a dedicated CMake build directory and import the
compiled extension directly from `modBHJet/python`. This avoids repeatedly
reinstalling with `pip`.

## Configure once

From `modBHJet`:

```bash
cmake -S . -B build-nb -G "Unix Makefiles" \
  -Dkariba_SOURCE_DIR=/Users/ruby/Desktop/University/phd/bhjet_mod/kariba
```

If CMake cannot find dependencies in your environment, pass their locations:

```bash
cmake -S . -B build-nb -G "Unix Makefiles" \
  -Dkariba_SOURCE_DIR=/Users/ruby/Desktop/University/phd/bhjet_mod/kariba \
  -Dnanobind_DIR=/Users/ruby/anaconda3/envs/mod_bhjet/lib/python3.12/site-packages/nanobind/cmake \
  -DPython_EXECUTABLE=/Users/ruby/anaconda3/envs/mod_bhjet/bin/python3.12 \
  -DGSL_ROOT_DIR=/Users/ruby/anaconda3/envs/mod_bhjet
```

## Rebuild after edits

```bash
cmake --build build-nb
```

The extension is written into:

```text
python/bhjet/bhjet.cpython-312-darwin.so
```

Restart Python or the Jupyter kernel after rebuilding the compiled extension.

## Import the local build in notebooks

```python
from pathlib import Path
import sys

ROOT = Path("/Users/ruby/Desktop/University/phd/bhjet_mod/modBHJet")
sys.path.insert(0, str(ROOT / "python"))

import bhjet
print(bhjet.__file__)
```
