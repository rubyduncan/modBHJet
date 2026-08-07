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

- `tutorials/01_bhjet_gammapy.ipynb` for the Binder-safe guided exercise
- `examples/gammapy/gammapy_intro.ipynb` for the advanced introduction
- `examples/gammapy/gammapy_bhjet_comparison.ipynb` for the advanced comparison
