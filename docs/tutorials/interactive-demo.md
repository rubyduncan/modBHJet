# Interactive Demo

The first interactive layer should use notebooks rather than a custom web app.
This keeps the documentation build reliable while still giving users a way to
change parameters and plot spectra.

## Binder launch

[![Launch Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rubyduncan/modBHJet/main?labpath=examples/starting/bhjet_pilot.ipynb)

The default public demo opens the pilot notebook:

```text
https://mybinder.org/v2/gh/rubyduncan/modBHJet/main?labpath=examples/starting/bhjet_pilot.ipynb
```

Recommended demo notebooks:

- `examples/starting/bhjet_pilot.ipynb`
- `examples/starting/gammapy_intro.ipynb`

## Future web app boundary

A future app should treat BHJet as a model service with:

- a parameter schema matching `BLJet`, `BHJet`, and target photon fields
- an observed energy grid input
- optional uploaded flux points
- a spectrum output table containing energy, total flux, and selected components

The web app should be separate from the static documentation site.
