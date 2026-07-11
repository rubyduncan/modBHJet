# Python API

The user-facing Python module is `bhjet`. It is generated with nanobind and
exports the main C++ classes.

```python
from bhjet import BLJet, IsoJet, BHJet, RadiationZone, TargetBlackBody
```

## Core classes

`BLJet`
: Magnetically driven jet dynamics.

`IsoJet`
: Quasi-isothermal jet dynamics.

`BHJet`
: Multi-zone orchestration and observed spectrum summation.

`RadiationZone`
: Local particle and radiation calculation.

`TargetBlackBody`
: Blackbody target photon field.

## Generated type stub

The generated stub file is included below as a compact signature reference.

```{literalinclude} ../../_modbhjet.pyi
:language: python
:caption: Generated nanobind/litgen stub
```
