# Particle Cutoff Comparison

`cutoff_type` is passed through the dynamics model into each radiation zone and
then into Kariba particle distributions.

```python
from bhjet import BLJet

for cutoff_type in [0, 1, 2]:
    dynamics = BLJet(cutoff_type=cutoff_type)
    dynamics.compute_jet_dynamics()
    print(cutoff_type, sorted(set(dynamics.get_cutoff_type_grid())))
```

Cutoff types:

- `0`: $\exp(-x)$
- `1`: $\exp(-x^2)$
- `2`: $\mathrm{sech}^2(x)$

For a spectral comparison, run the same `BHJet` setup for each cutoff type and
plot the post-dissipation synchrotron or total flux arrays on the same axes.
