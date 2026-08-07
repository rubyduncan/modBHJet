"""Minimal post-build check for the public BHJet tutorial environment."""

from pathlib import Path

import numpy as np
from astropy.table import Table

from bhjet import BHJet, BLJet


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATA_FILE = REPOSITORY_ROOT / "examples" / "gammapy" / "data" / "ngc4261" / "ngc4261_dnde.csv"

if not DATA_FILE.is_file():
    raise FileNotFoundError(f"Bundled tutorial data is missing: {DATA_FILE}")

table = Table.read(DATA_FILE, format="ascii.csv", delimiter=" ")
if len(table) == 0 or not {"e_ref", "dnde", "dnde_err"}.issubset(table.colnames):
    raise ValueError("The NGC 4261 tutorial table is not usable.")

dynamics = BLJet(dlgz=0.2, verbosity_level=0)
jet = BHJet(verbosity_level=0)
jet.init_jet_dynamics(dynamics)

energy_eV = np.logspace(-6, 9, 24)
jet.compute_full_jet(photon_energy_grid=(energy_eV / 6.242e11).tolist())
flux = np.asarray(jet.get_observed_photon_flux_total(), dtype=float)

if flux.shape != energy_eV.shape or not np.isfinite(flux).all():
    raise RuntimeError("BHJet did not return a finite spectrum on the tutorial grid.")

print(f"BHJet tutorial smoke test passed with {len(table)} bundled data points.")
