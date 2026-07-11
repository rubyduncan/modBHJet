# Import and Kernel Checks

Most notebook import problems come from loading an old installed extension
instead of the locally rebuilt one.

Use this cell at the top of development notebooks:

```python
from pathlib import Path
import sys

candidates = []
for base in [Path.cwd(), *Path.cwd().parents]:
    candidates.extend([base / "python", base / "modBHJet" / "python"])

local_python = next((c for c in candidates if (c / "bhjet" / "__init__.py").exists()), None)
if local_python is None:
    raise RuntimeError("Could not find local modBHJet/python package directory")

sys.path.insert(0, str(local_python))

import bhjet
from bhjet import BLJet, BHJet

print("Using bhjet from:", bhjet.__file__)
print("BLJet supports cutoff_type:", "cutoff_type" in (BLJet.__init__.__doc__ or ""))
```

Expected output should point inside `modBHJet/python/bhjet`, not an unrelated
`site-packages` directory.
