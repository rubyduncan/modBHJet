# Full BHJet Spectrum

`BHJet` combines jet dynamics and radiation-zone calculations.

```python
import numpy as np
import matplotlib.pyplot as plt
from bhjet import BLJet, BHJet

erg_per_ev = 1.0 / 6.242e11
energy_ev = np.logspace(-8, 19, 200)

dynamics = BLJet(
    mass_bh=1e9,
    jet_power_eddington=1e-5,
    z_dissipation=1e2,
    z_max_calculation=1e6,
    gamma_final=15,
    fraction_nonthermal_electrons=0.1,
    cutoff_type=0,
)

model = BHJet(
    theta_obs=17,
    distance=16e3,
    redshift=0.00428,
    include_counterjet=True,
    compton_threshold=1e-5,
)
model.init_jet_dynamics(dynamics)
model.compute_full_jet((energy_ev * erg_per_ev).tolist())

flux = np.asarray(model.get_observed_photon_flux_total())

plt.loglog(energy_ev, energy_ev * erg_per_ev * flux)
plt.xlabel("Observed energy [eV]")
plt.ylabel(r"Observed energy flux [erg cm$^{-2}$ s$^{-1}$]")
```

Use the integrated component getters to separate pre- and post-dissipation
regions:

```python
post_syn = model.get_observed_photon_integrated_flux_electron_cyclosyn(
    dynamics.z_dissipation * dynamics.r_g,
    dynamics.z_max_calculation * dynamics.r_g,
)
```
