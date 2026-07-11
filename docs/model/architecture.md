# Architecture

modBHJet separates the model into three layers:

1. **Jet dynamics** compute spatial profiles along the jet.
2. **Radiation zones** compute local particle distributions and emission.
3. **BHJet orchestration** builds zones, adds target fields, and sums the
   observed spectrum.

## Main classes

`BLJet`
: Magnetically driven jet dynamics model. It computes the zone grid, bulk
  Lorentz factors, magnetic field, particle densities, temperatures, and
  nonthermal-particle controls.

`IsoJet`
: Quasi-isothermal jet dynamics model with the same `JetDynamics` interface.

`BHJet`
: High-level model wrapper. It accepts a `JetDynamics` instance, constructs
  `RadiationZone` objects, computes emission, and exposes total or
  region-integrated observed fluxes.

`RadiationZone`
: Local homogeneous emission zone. It constructs thermal, mixed, broken
  power-law, or power-law electron distributions and computes cyclosynchrotron
  and Compton components.

`TargetPhotonField`
: Abstract interface for external target photon fields.

`TargetBlackBody`
: Concrete blackbody target field used both as a Compton seed field and,
  optionally, as a direct component of the total observed flux.

## Backend

Kariba supplies particle distributions and radiative processes. modBHJet wraps
those lower-level classes into a multi-zone jet workflow.
