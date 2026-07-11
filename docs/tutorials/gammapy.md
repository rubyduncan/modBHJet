# Gammapy Wrapper

The package includes a Gammapy-facing wrapper in `bhjet.gammapy`.

```python
import astropy.units as u
import numpy as np
from bhjet.gammapy import BHJetSpectralModel

model = BHJetSpectralModel()
energy = np.logspace(-8, 2, 20) * u.TeV
flux = model(energy)
```

The wrapper converts Gammapy energy inputs to erg, computes a BHJet spectrum,
and returns differential photon flux with units compatible with Gammapy.

Use the starting notebooks for concrete workflows:

- `examples/starting/gammapy_intro.ipynb`
- `examples/starting/gammapy_bhjet_comparison.ipynb`
