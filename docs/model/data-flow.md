# Data Flow

The full model calculation follows this sequence:

1. Instantiate a dynamics model, usually `BLJet`.
2. Instantiate `BHJet` with observer/source settings.
3. Attach the dynamics model with `BHJet.init_jet_dynamics(...)`.
4. Optionally add target photon fields.
5. Call `BHJet.compute_full_jet(photon_energy_grid)`.
6. Read total or component fluxes from `BHJet` and its `radiation_zones`.

```python
import numpy as np
from bhjet import BLJet, BHJet

erg_per_ev = 1.0 / 6.242e11
energy_ev = np.logspace(-8, 19, 200)

dynamics = BLJet(cutoff_type=0)
model = BHJet(theta_obs=17, distance=16e3, redshift=0.00428)
model.init_jet_dynamics(dynamics)
model.compute_full_jet((energy_ev * erg_per_ev).tolist())

flux = model.get_observed_photon_flux_total()
```

## Units

- Photon energies passed to `compute_full_jet` are in erg.
- Many input distances inside dynamics classes are in gravitational radii.
- `distance` in `BHJet` is in kpc.
- Observed photon flux getters return number-flux style arrays compatible with
  the internal Kariba luminosity conversion.
