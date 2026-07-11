# BLJet Dynamics

Use `BLJet` directly when you want to inspect the physical jet profile before
computing radiation.

```python
import numpy as np
import matplotlib.pyplot as plt
from bhjet import BLJet

jet = BLJet(
    mass_bh=5e8,
    jet_power_eddington=2e-3,
    z_jet_launching=2,
    r_initial=5,
    z_end_of_acceleration=100,
    z_dissipation=50,
    z_max_calculation=1e6,
    gamma_final=2,
    cutoff_type=0,
)

jet.compute_jet_dynamics()

z_rg = np.asarray(jet.get_z_center_grid()) / jet.r_g
beta_gamma = jet.get_beta_gamma_grid()

plt.loglog(z_rg, beta_gamma)
plt.xlabel(r"z [$r_g$]")
plt.ylabel(r"$\beta\gamma$")
```

Useful getters include:

- `get_radius_grid()`
- `get_magnetic_field_grid()`
- `get_electron_density_grid()`
- `get_fraction_nonthermal_electrons_grid()`
- `get_cutoff_type_grid()`
