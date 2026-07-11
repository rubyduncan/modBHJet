# Physics Notes

This section summarizes the model assumptions that are most useful when
choosing parameters or interpreting spectra.

## Zone discretization

The jet is split into independent zones along the jet axis. The dynamics model
sets each zone's position, height, radius, bulk speed, magnetic field, particle
density, temperature, and nonthermal-particle parameters. Radiation is computed
zone-by-zone and then interpolated onto the requested observed energy grid.

## Particle distributions

The electron distribution in each zone depends on
`fraction_nonthermal_electrons`:

- `0`: thermal distribution
- `0 < f < 0.5`: mixed thermal plus power-law tail
- `0.5 <= f < 1`: broken power-law approximation
- `1`: pure power law

The high-energy cutoff shape is controlled by `cutoff_type`:

- `0`: exponential, $\exp(-x)$
- `1`: super-exponential, $\exp(-x^2)$
- `2`: $\mathrm{sech}^2(x)$

where $x = p / p_\mathrm{cut}$.

## Emission components

`RadiationZone` computes cyclosynchrotron emission and uses the resulting
synchrotron photons as seed photons for Compton scattering. Additional target
photon fields can be attached to `BHJet` and are transformed into each zone's
frame before the particle and radiation calculations.

## Observed flux

Each zone computes emission in the observer frame, accounting for bulk beaming
and the counterjet setting. `BHJet` sums all zones on the energy grid requested
by the caller.
