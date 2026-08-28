"""
bhjet_explorer.py
-----------------
Unified interactive BHJet SED explorer.

Supports both the pure bhjet classes and the 3ML BHJetPlugin wrapper.
3ML is imported optionally — if not available, the mode toggle is hidden.

Usage in a notebook:
    from bhjet_explorer import BHJetExplorer
    explorer = BHJetExplorer()
    explorer.display()          # SED + jet shape + parameter tabs
    explorer.display_zones()    # zone inspector (separate figure + slider)
"""

from html import escape

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import ipywidgets as widgets
from ipywidgets import VBox, HBox, Tab
from IPython.display import display

from bhjet import BLJet, BHJet, RadiationZone, TargetBlackBody

# ── optional 3ML import ───────────────────────────────────────────────────────
try:
    from bhjet.threeML_plugins.BHJetPlugin import BHJetPlugin
    from bhjet.threeML_plugins.TargetBlackBody import TargetBlackBody as BHJetTargetBB
    _THREEML_AVAILABLE = True
except ImportError:
    _THREEML_AVAILABLE = False

# ── constants ────────────────────────────────────────────────────────────────
erg2eV  = 6.242e+11
eV2erg  = 1.0 / erg2eV
c_cm_s  = 3e10
h_eV_s  = 4.135667696e-15

_MOM_GRID = np.logspace(-8, 8, 200) / c_cm_s   # g cm/s

# Representative frequencies for the broad observing bands shown above the
# tutorial SED. They are guide labels, rather than sharp band boundaries.
OBSERVING_BANDS_HZ = {
    "radio": 2.4e9,
    "sub-mm": 2.4e11,
    "infrared": 2.4e13,
    "optical/UV": 4.8e15,
    "X-ray": 2.4e18,
    "gamma-ray": 2.4e22,
}


def energy_eV_to_frequency_hz(energy_eV):
    """Convert photon energy in eV to frequency in Hz."""
    return np.asarray(energy_eV) / h_eV_s


def add_observing_band_labels(axis):
    """Add the fixed broad-band guide labels used by the tutorial SED plots."""
    return [
        axis.text(
            frequency_hz,
            1.02,
            label,
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9,
            clip_on=False,
        )
        for label, frequency_hz in OBSERVING_BANDS_HZ.items()
    ]

# ── tutorial interface customisation ─────────────────────────────────────────

# Leave VISIBLE_PARAMETERS as None to show every parameter. To show only a
# smaller teaching-focused selection, replace None with a set of names from
# PARAMETER_LABELS, for example: {"lg_10_jet_power_eddington", "lg_10_z_dissipation"}.
VISIBLE_PARAMETERS = {
    "lg_10_mass_bh", "lg_10_jet_power_eddington", "lg_10_r_initial",
    "lg_10_z_end_of_acceleration", "link_particle_to_bulk_acceleration",
    "lg_10_z_dissipation", "lg_10_z_max_calculation",
    "lg_10_sigma_final", "gamma_final", "lg_10_electron_temperature_jet_base",
    "lg_10_fraction_nonthermal_electrons", "lg_10_factor_break_electrons",
    "lg_10_factor_max_energy_electrons", "index_injected_electrons", "theta_obs",
    "lg_10_distance", 
    "bb_enable", "bb_add_to_total", "lg_10_bb_temperature", "lg_10_bb_energy_density", 
    "lg_10_bb_luminosity",
}

# Edit values (the text on the right) to change labels in the tutorial widget.
PARAMETER_LABELS = {
    "lg_10_mass_bh": "log10 mass_bh [solar masses]",
    "lg_10_jet_power_eddington": "log10 jet power eddington [L_Edd]",
    "theta_obs": "Viewing angle [degrees]",
    "lg_10_distance": "log10 distance [kpc]",
    "redshift": "Redshift",
    "lg_10_z_jet_launching": "log10 jet launching height [rg]",
    "lg_10_r_initial": "log10 jet base initial radius [rg]",
    "lg_10_z_end_of_acceleration": "log10 dist (z) end of bulk acceleration [rg]",
    "link_particle_to_bulk_acceleration": "z_dissipation = z_acceleration",
    "lg_10_z_dissipation": "log10 dist (z) particle dissipation [rg]",
    "lg_10_z_max_calculation": "log10 dist (z) max calculation [rg]",
    "lg_10_sigma_final": "log10 final magnetisation (sigma_final)",
    "gamma_final": "Bulk Lorentz factor of jet",
    "plasma_beta_jet_base": "Plasma beta at jet base",
    "lg_10_electron_temperature_jet_base": "log10 electron temp. at jet base [keV]",
    "gamma_acceleration_exponent": "Acceleration exponent",
    "lg_10_opening_angle_constant": "log10 opening-angle constant",
    "lg_10_fraction_nonthermal_electrons": "log10 non-thermal electron fraction",
    "lg_10_factor_break_electrons": "log10 electron break-energy factor",
    "lg_10_factor_max_energy_electrons": "log10 maximum-electron-energy factor",
    "index_injected_electrons": "Injected electron index",
    "lg_10_fraction_nonthermal_protons": "log10 non-thermal proton fraction",
    "lg_10_factor_break_protons": "log10 proton break-energy factor",
    "lg_10_factor_max_energy_protons": "log10 maximum-proton-energy factor",
    "index_injected_protons": "Injected proton index",
    "bb_enable": "Enable black-body target",
    "bb_add_to_total": "Add black-body flux to total",
    "lg_10_bb_temperature": "log10 black-body temperature [keV]",
    "lg_10_bb_energy_density": "log10 black-body energy density [erg cm^-3]",
    "lg_10_bb_luminosity": "log10 black-body luminosity [erg s^-1]",
    "include_counterjet": "Include counterjet",
    "lg_10_compton_threshold": "log10 Compton threshold",
    "dlgz": "Grid spacing",
}

# Units displayed alongside the linear-value readout of each log10 slider.
# Parameters absent from this mapping are dimensionless.
PARAMETER_VALUE_UNITS = {
    "lg_10_mass_bh": " M_sun",
    "lg_10_jet_power_eddington": " L_Edd",
    "lg_10_distance": " kpc",
    "lg_10_z_jet_launching": " r_g",
    "lg_10_r_initial": " r_g",
    "lg_10_z_end_of_acceleration": " r_g",
    "lg_10_z_dissipation": " r_g",
    "lg_10_z_max_calculation": " r_g",
    "lg_10_electron_temperature_jet_base": " keV",
    "lg_10_bb_temperature": " keV",
    "lg_10_bb_energy_density": " erg cm^-3",
    "lg_10_bb_luminosity": " erg s^-1",
}

# Parameter names in this set are given a bold label in the widget panel.
EMPHASISED_PARAMETERS = frozenset() #{"lg_10_jet_power_eddington"}

PARAMETER_LIMITS = {
    "lg_10_jet_power_eddington": (-9.0, -1.0, 0.05),
    "theta_obs": (1, 89, 2),
    "index_injected_electrons": (1.2, 5, 0.2),
}

# The optional 3ML plugin still calls logarithmic parameters "lg_...". Map
# tutorial-facing "lg_10_..." names to those internal plugin names.
_UI_PLUGIN_PARAMETER_NAMES = {
    "lg_10_mass_bh", "lg_10_jet_power_eddington",
    "lg_10_z_jet_launching", "lg_10_r_initial",
    "lg_10_z_end_of_acceleration", "lg_10_z_dissipation", "lg_10_z_max_calculation",
    "lg_10_sigma_final", "plasma_beta_jet_base",
    "lg_10_electron_temperature_jet_base", "gamma_acceleration_exponent",
    "lg_10_opening_angle_constant",
    "lg_10_fraction_nonthermal_electrons", "lg_10_fraction_nonthermal_protons",
    "lg_10_factor_break_electrons", "lg_10_factor_break_protons",
    "lg_10_factor_max_energy_electrons", "lg_10_factor_max_energy_protons",
    "index_injected_electrons", "index_injected_protons",
    "theta_obs", "lg_10_distance", "redshift",
}
_PLUGIN_PARAMETER_NAMES = {
    name: name.replace("lg_10_", "lg_") for name in _UI_PLUGIN_PARAMETER_NAMES
}


class OneZoneFittingExplorer:
    """Interactive single-emission-region SED fitter for a tutorial notebook."""

    def __init__(
        self,
        distance,
        redshift,
        theta_obs,
        Emin_eV=1e-9,
        Emax_eV=1e13,
        n_energy=220,
        defaults=None,
        observation_energy_eV=None,
        observation_sed=None,
        observation_sed_error=None,
        observation_label="observations",
    ):
        self.distance = distance
        self.redshift = redshift
        self.theta_obs = theta_obs
        self.E_eV = np.logspace(np.log10(Emin_eV), np.log10(Emax_eV), n_energy)
        self._defaults = {
            "magnetic_field": 0.2,
            "radius": 1e16,
            "height": 1e16,
            "geometry": "sphere",
            "bulk_momentum": 2.0,
            "electron_number_density": 1e3,
            "electron_temperature": 900.0,
            "fraction_nonthermal_electrons": 0.1,
            "factor_break_electrons": 0.1,
            "factor_max_energy_electrons": 1e-7,
            "index_injected_electrons": 1.8,
        }
        if defaults is not None:
            self._defaults.update(defaults)
        self._starting_defaults = self._defaults.copy()
        self._resetting = False

        self._observation_energy_eV = None
        self._observation_sed = None
        self._observation_sed_error = None
        self._observation_label = observation_label
        if observation_energy_eV is not None and observation_sed is not None:
            self._observation_energy_eV = np.asarray(observation_energy_eV, dtype=float)
            self._observation_sed = np.asarray(observation_sed, dtype=float)
            self._observation_sed_error = (
                None if observation_sed_error is None
                else np.asarray(observation_sed_error, dtype=float)
            )

        self._build_figure()
        self._build_widgets()
        self.calculate()

    def _build_figure(self):
        plt.ioff()
        self.fig, self.ax = plt.subplots(figsize=(10, 5))
        self.ax.set(
            xscale="log", yscale="log",
            xlim=(energy_eV_to_frequency_hz(self.E_eV[0]),
                  energy_eV_to_frequency_hz(self.E_eV[-1])),
            ylim=(1e-18, 1e-10),
            xlabel="Observed frequency [Hz]",
            ylabel=r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]",
        )
        self.ax.grid(alpha=0.3)
        self._band_labels = add_observing_band_labels(self.ax)
        if self._observation_energy_eV is not None:
            self.ax.errorbar(
                energy_eV_to_frequency_hz(self._observation_energy_eV), self._observation_sed,
                yerr=self._observation_sed_error,
                fmt="o", ms=4, color="black", label=self._observation_label,
            )
        self._model_line, = self.ax.loglog([np.nan], [np.nan], c="k", lw=2,
                                            label="one-zone total")
        self._synchrotron_line, = self.ax.loglog([np.nan], [np.nan], c="tab:blue",
                                                   ls="--", lw=1.5, label="synchrotron")
        self._ssc_line, = self.ax.loglog([np.nan], [np.nan], c="tab:red",
                                          ls="--", lw=1.5, label="SSC")
        self.ax.legend()
        plt.ion()

    def _build_widgets(self):
        d = self._defaults

        def log_slider(description, value, lower, upper):
            return widgets.FloatLogSlider(
                value=value, base=10, min=lower, max=upper, step=0.05,
                description=description, continuous_update=False,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="540px"),
            )

        def linear_slider(description, value, lower, upper, step):
            return widgets.FloatSlider(
                value=value, min=lower, max=upper, step=step,
                description=description, continuous_update=False,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="540px"),
            )

        self._widgets = {
            "geometry": widgets.Dropdown(
                options=["sphere", "cylinder"], value=d["geometry"],
                description="Geometry", style={"description_width": "250px"},
                layout=widgets.Layout(width="540px"),
            ),
            "magnetic_field": log_slider("Magnetic field [G, log scale]", d["magnetic_field"], -3, 2),
            "radius": log_slider("Radius [cm, log scale]", d["radius"], 12, 20),
            "height": log_slider("Height [cm, log scale]", d["height"], 12, 20),
            "bulk_momentum": linear_slider("Bulk momentum", d["bulk_momentum"], 1, 20, 0.1),
            "electron_number_density": log_slider("Electron density [cm$^{-3}$, log scale]", d["electron_number_density"], 0, 6),
            "electron_temperature": log_slider("Electron temperature [keV, log scale]", d["electron_temperature"], 1, 4),
            "fraction_nonthermal_electrons": log_slider("Non-thermal electron fraction [log scale]", d["fraction_nonthermal_electrons"], -2, 0),
            "factor_break_electrons": log_slider("Electron break factor [log scale]", d["factor_break_electrons"], -2, 0),
            "factor_max_energy_electrons": log_slider("Maximum-energy factor [log scale]", d["factor_max_energy_electrons"], -9, 0),
            "index_injected_electrons": linear_slider("Injected electron index", d["index_injected_electrons"], 1.5, 3.5, 0.05),
        }
        for widget in self._widgets.values():
            widget.observe(self._on_widget_change, names="value")

        emission_region = VBox([
            self._widgets["geometry"], self._widgets["magnetic_field"],
            self._widgets["radius"], self._widgets["height"],
            self._widgets["bulk_momentum"],
        ])
        electrons = VBox([
            self._widgets["electron_number_density"],
            self._widgets["electron_temperature"],
            self._widgets["fraction_nonthermal_electrons"],
            self._widgets["factor_break_electrons"],
            self._widgets["factor_max_energy_electrons"],
            self._widgets["index_injected_electrons"],
        ])
        self._tab = Tab(children=[emission_region, electrons])
        self._tab.set_title(0, "Emission region")
        self._tab.set_title(1, "Electrons")

        self._reset_button = widgets.Button(description="Reset to starting model")
        self._status = widgets.Label(value="Model updates when you release a control.")
        self._reset_button.on_click(lambda _: self.reset())
        self._controls = HBox([self._reset_button, self._status])

    def _parameters(self):
        return {name: widget.value for name, widget in self._widgets.items()}

    def _on_widget_change(self, change=None):
        if not self._resetting:
            self.calculate()

    def calculate(self):
        """Calculate and redraw the one-zone model."""
        zone = RadiationZone(
            distance=self.distance, redshift=self.redshift, theta_obs=self.theta_obs,
            include_counterjet=True, compton_threshold=1e-5, **self._parameters(),
        )
        zone.compute_particles()
        zone.compute_radiation((self.E_eV * eV2erg).tolist())
        self._zone = zone

        energy_erg = np.asarray(zone.get_observed_photon_energy_grid_total(), dtype=float)
        self._energy_eV = energy_erg * erg2eV
        self._sed = energy_erg * np.asarray(zone.get_observed_photon_flux_total(), dtype=float)

        synchrotron_energy_erg = np.asarray(
            zone.get_observed_photon_energy_grid_electron_cyclosyn(), dtype=float
        )
        self._synchrotron_energy_eV = synchrotron_energy_erg * erg2eV
        self._synchrotron_sed = synchrotron_energy_erg * np.asarray(
            zone.get_observed_photon_flux_electron_cyclosyn(), dtype=float
        )

        ssc_energy_erg = np.asarray(
            zone.get_observed_photon_energy_grid_electron_compton(), dtype=float
        )
        self._ssc_energy_eV = ssc_energy_erg * erg2eV
        self._ssc_sed = ssc_energy_erg * np.asarray(
            zone.get_observed_photon_flux_electron_compton(), dtype=float
        )

        self._model_line.set_data(energy_eV_to_frequency_hz(self._energy_eV), self._sed)
        self._synchrotron_line.set_data(
            energy_eV_to_frequency_hz(self._synchrotron_energy_eV), self._synchrotron_sed
        )
        self._ssc_line.set_data(energy_eV_to_frequency_hz(self._ssc_energy_eV), self._ssc_sed)
        self.fig.canvas.draw_idle()
        self._status.value = "Model updated."

    def reset(self):
        """Restore the tutorial starting values."""
        self._resetting = True
        for name, widget in self._widgets.items():
            widget.value = self._starting_defaults[name]
        self._resetting = False
        self.calculate()

    def get_total_sed(self):
        """Return the current model as energy [eV] and nu F_nu [cgs]."""
        return self._energy_eV, self._sed

    def get_component_seds(self):
        """Return the current synchrotron and SSC spectra in tutorial units."""
        return {
            "synchrotron": (self._synchrotron_energy_eV, self._synchrotron_sed),
            "ssc": (self._ssc_energy_eV, self._ssc_sed),
        }

    def get_parameters(self):
        """Return the current values selected in the widgets."""
        return self._parameters().copy()

    def set_parameters(self, parameters):
        """Apply saved parameters to the widgets and recalculate the model."""
        self._resetting = True
        try:
            for name, value in parameters.items():
                if name in self._widgets:
                    self._widgets[name].value = value
        finally:
            self._resetting = False
        self.calculate()

    def display(self):
        """Render the plot first, with tabbed controls directly below it."""
        display(VBox([self.fig.canvas, self._controls, self._tab]))


class BHJetExplorer:
    """
    Interactive BHJet SED explorer (pure classes + optional 3ML mode).

    Parameters
    ----------
    verbosity_level : int
    Emin_eV, Emax_eV : float
    n_energy : int
    """

    def __init__(
        self,
        verbosity_level=0,
        Emin_eV=1e-6,
        Emax_eV=1e12,
        n_energy=300,
        defaults=None,
        observation_energy_eV=None,
        observation_sed=None,
        observation_sed_error=None,
        observation_label="observations",
        manual_update=False,
    ):
        self._verbosity = verbosity_level
        self._manual_update = manual_update
        self._loading_parameters = False
        self._observation_energy_eV = None
        self._observation_sed = None
        self._observation_sed_error = None
        self._observation_label = observation_label

        if observation_energy_eV is not None and observation_sed is not None:
            self._observation_energy_eV = np.asarray(observation_energy_eV, dtype=float)
            self._observation_sed = np.asarray(observation_sed, dtype=float)
            self._observation_sed_error = (
                None if observation_sed_error is None
                else np.asarray(observation_sed_error, dtype=float)
            )

        self.E_eV  = np.logspace(np.log10(Emin_eV), np.log10(Emax_eV), n_energy)
        self.E_erg = self.E_eV * eV2erg
        self.E_keV = self.E_eV * 1e-3

        self._defaults = dict(
            theta_obs                          = 17.0,
            lg_10_distance                     = 4.0,
            redshift                           = 1e-3,
            lg_10_mass_bh                      = 9.0,
            lg_10_jet_power_eddington          = -5.0,
            lg_10_z_jet_launching              = 0.3,
            lg_10_r_initial                    = 0.7,
            lg_10_z_end_of_acceleration        = 2.5,
            lg_10_z_dissipation                = 2.5,
            link_particle_to_bulk_acceleration = True,
            lg_10_z_max_calculation            = 5.5,
            lg_10_sigma_final                  = 0.0,
            gamma_final                        = 5.0,
            plasma_beta_jet_base               = 0.02,
            calc_pair_content_from_plasma_beta = False,
            lg_10_electron_temperature_jet_base = 3.5,
            gamma_acceleration_exponent        = 0.5,
            lg_10_opening_angle_constant       = -0.82,
            lg_10_fraction_nonthermal_electrons = -1.0,
            lg_10_fraction_nonthermal_protons  = -1.0,
            lg_10_factor_break_electrons       = 0.0,
            lg_10_factor_break_protons         = 0.0,
            lg_10_factor_max_energy_electrons  = 0.0,
            lg_10_factor_max_energy_protons    = 0.0,
            index_injected_electrons           = 2.0,
            index_injected_protons             = 2.0,
            bb_enable                          = False,
            bb_add_to_total                    = True,
            lg_10_bb_temperature               = -3.5,
            lg_10_bb_energy_density            = -9.0,
            lg_10_bb_luminosity                = 44.0,
            include_counterjet                 = True,
            lg_10_compton_threshold            = -5.0,
            dlgz                               = 0.1,
        )
        if defaults is not None:
            self._defaults.update(defaults)
        self._starting_defaults = self._defaults.copy()

        self._mode             = "pure"
        self._plugin           = None
        self._bb_plugin_target = None
        self._bb_added         = False
        self.bljet             = None
        self.bhjet             = None
        self._bb_target        = None

        self._build_models()
        self._build_figure()
        self._build_widgets()

        for w in self._flat_widgets.values():
            if hasattr(w, "continuous_update"):
                w.continuous_update = False
            w.observe(self._on_widget_change, names="value")

    # ── model construction ────────────────────────────────────────────────────

    def _build_models(self):
        d = self._defaults
        self._build_pure_models(d)
        if _THREEML_AVAILABLE:
            self._build_3ml_plugin(d)

    def _build_pure_models(self, d):
        self.bljet = BLJet(
            mass_bh                            = 10**d["lg_10_mass_bh"],
            jet_power_eddington                = 10**d["lg_10_jet_power_eddington"],
            z_jet_launching                    = 10**d["lg_10_z_jet_launching"],
            r_initial                          = 10**d["lg_10_r_initial"],
            z_end_of_acceleration              = 10**d["lg_10_z_end_of_acceleration"],
            z_dissipation                      = 10**d["lg_10_z_dissipation"],
            z_max_calculation                  = 10**d["lg_10_z_max_calculation"],
            sigma_final                        = 10**d["lg_10_sigma_final"],
            gamma_final                        = d["gamma_final"],
            plasma_beta_jet_base               = d["plasma_beta_jet_base"],
            calc_pair_content_from_plasma_beta = d["calc_pair_content_from_plasma_beta"],
            electron_temperature_jet_base      = 10**d["lg_10_electron_temperature_jet_base"],
            gamma_acceleration_exponent        = d["gamma_acceleration_exponent"],
            opening_angle_constant             = 10**d["lg_10_opening_angle_constant"],
            fraction_nonthermal_electrons      = 10**d["lg_10_fraction_nonthermal_electrons"],
            fraction_nonthermal_protons        = 10**d["lg_10_fraction_nonthermal_protons"],
            factor_break_electrons             = 10**d["lg_10_factor_break_electrons"],
            factor_break_protons               = 10**d["lg_10_factor_break_protons"],
            factor_max_energy_electrons        = 10**d["lg_10_factor_max_energy_electrons"],
            factor_max_energy_protons          = 10**d["lg_10_factor_max_energy_protons"],
            index_injected_electrons           = d["index_injected_electrons"],
            index_injected_protons             = d["index_injected_protons"],
            dlgz                               = d["dlgz"],
            verbosity_level                    = self._verbosity,
        )
        self.bhjet = BHJet(
            theta_obs          = d["theta_obs"],
            distance           = 10**d["lg_10_distance"],
            redshift           = d["redshift"],
            include_counterjet = d["include_counterjet"],
            compton_threshold  = 10**d["lg_10_compton_threshold"],
            verbosity_level    = self._verbosity,
        )
        self.bhjet.init_jet_dynamics(self.bljet)
        self.bhjet.clear_target_photon_fields()
        self._bb_target = None
        self._run_pure()

    def _build_3ml_plugin(self, d):
        self._plugin = BHJetPlugin()
        self._plugin._setup(
            Emin_eV            = self.E_eV[0],
            Emax_eV            = self.E_eV[-1],
            dlgz               = d["dlgz"],
            compton_threshold  = 10**d["lg_10_compton_threshold"],
            verbosity_level    = self._verbosity,
            include_counterjet = d["include_counterjet"],
            cache_rtol         = 0,
        )
        for ui_name, plugin_name in _PLUGIN_PARAMETER_NAMES.items():
            if plugin_name in self._plugin.parameters:
                self._plugin.parameters[plugin_name].value = d[ui_name]
        if "lg_gamma_final" in self._plugin.parameters:
            self._plugin.parameters["lg_gamma_final"].value = np.log10(d["gamma_final"])
        self._bb_plugin_target = None
        self._bb_added         = False
        self._plugin(self.E_keV)

    def _run_pure(self):
        self.bhjet.compute_full_jet(photon_energy_grid=self.E_erg.tolist())

    def _run_3ml(self):
        self._plugin(self.E_keV)

    # ── active model references ───────────────────────────────────────────────

    @property
    def _active_bljet(self):
        if self._mode == "3ml" and self._plugin is not None:
            return self._plugin.bljet
        return self.bljet

    @property
    def _active_bhjet(self):
        if self._mode == "3ml" and self._plugin is not None:
            return self._plugin.bhjet
        return self.bhjet

    @property
    def _active_bb(self):
        if self._mode == "3ml":
            return self._bb_plugin_target.bb if self._bb_plugin_target is not None else None
        return self._bb_target

    # ── SED figure ────────────────────────────────────────────────────────────

    def _build_figure(self):
        plt.ioff()
        self.fig = plt.figure(figsize=(12, 5))
        gs = gridspec.GridSpec(1, 2, width_ratios=[8.5, 1.5], wspace=0.05,
                               figure=self.fig)
        self.ax  = self.fig.add_subplot(gs[0])
        self.cax = self.fig.add_subplot(gs[1])

        self.ax.set_xlabel("Observed frequency [Hz]", fontsize=12)
        self.ax.set_ylabel(r"$\nu F_\nu$  [erg cm$^{-2}$ s$^{-1}$]", fontsize=12)
        self.ax.set_xscale("log")
        self.ax.set_yscale("log")
        self.ax.set_xlim(energy_eV_to_frequency_hz(self.E_eV[0]),
                         energy_eV_to_frequency_hz(self.E_eV[-1]))
        self.ax.set_ylim(1e-20, 1e-8)
        self.ax.set_aspect("equal")
        self.ax.grid(alpha=0.3)

        # A widget-backed Matplotlib canvas cannot reliably render an added
        # secondary axis. These fixed labels therefore form a visual top-axis
        # guide, while keeping the SED itself as one interactive canvas.
        self._band_labels = add_observing_band_labels(self.ax)

        self.cax.set_ylabel(r"$\log_{10}\,z/r_g$")
        self.cax.set_xlabel(r"$\log_{10}\,r/r_g$")
        self.cax.yaxis.tick_right()
        self.cax.yaxis.set_label_position("right")

        dummy = np.full_like(self.E_eV, np.nan)
        self._lines = {}
        frequency_hz = energy_eV_to_frequency_hz(self.E_eV)
        self._lines["total"],   = self.ax.loglog(frequency_hz, dummy, c="k",          lw=2,    label="total")
        self._lines["presyn"],  = self.ax.loglog(frequency_hz, dummy, c="tab:blue",   ls="--", lw=1.5, label="pre-acc. syn.")
        self._lines["precom"],  = self.ax.loglog(frequency_hz, dummy, c="tab:green",  ls="--", lw=1.5, label="pre-acc. IC")
        self._lines["postsyn"], = self.ax.loglog(frequency_hz, dummy, c="tab:orange", ls=":",  lw=1.5, label="post-acc. syn.")
        self._lines["postcom"], = self.ax.loglog(frequency_hz, dummy, c="tab:red",    ls=":",  lw=1.5, label="post-acc. IC")
        self._lines["bb"],      = self.ax.loglog(frequency_hz, dummy, c="tab:purple", ls="-.", lw=1.5, label="BB target")
        self._lines["bb"].set_visible(False)
        if self._observation_energy_eV is not None:
            self._observation_artist = self.ax.errorbar(
                energy_eV_to_frequency_hz(self._observation_energy_eV),
                self._observation_sed,
                yerr=self._observation_sed_error,
                fmt="o",
                ms=4,
                color="black",
                label=self._observation_label,
            )
        else:
            self._observation_artist = None
        self.ax.legend(fontsize=10)

        self._jet_shape_line, = self.cax.plot([], [], c="k",      lw=1.5, marker="_")
        self._zacc_line       = self.cax.axhline(np.nan, ls=":",  c="grey",       lw=1, label=r"end acc.")
        self._zdiss_line      = self.cax.axhline(np.nan, ls=":", c="tab:orange", lw=1, label=r"dissip.")
        self._zmax_line       = self.cax.axhline(np.nan, ls="-",  c="k",          lw=1)
        self._zone_marker,    = self.cax.plot([], [], marker="o", ms=6,
                                               c="tab:red", zorder=5)
        self.cax.legend(fontsize=10, handlelength=0.5, handletextpad=0.2)
        self._draw_spectra()
        self._draw_jet_shape()
        plt.ion()

    # ── SED draw helpers ──────────────────────────────────────────────────────

    def _draw_spectra(self):
        bj = self._active_bhjet
        bl = self._active_bljet
        bb = self._active_bb

        E         = np.array(bj.get_observed_photon_energy_grid())
        frequency_hz = energy_eV_to_frequency_hz(E * erg2eV)
        z_diss_cm = bl.z_dissipation * bl.r_g
        z_max_cm  = bl.z_max_calculation * bl.r_g

        self._lines["total"].set_data(
            frequency_hz, E * np.array(bj.get_observed_photon_flux_total()))
        self._lines["presyn"].set_data(
            frequency_hz, E * np.array(bj.get_observed_photon_integrated_flux_electron_cyclosyn(0.0, z_diss_cm)))
        self._lines["precom"].set_data(
            frequency_hz, E * np.array(bj.get_observed_photon_integrated_flux_electron_compton(0.0, z_diss_cm)))
        self._lines["postsyn"].set_data(
            frequency_hz, E * np.array(bj.get_observed_photon_integrated_flux_electron_cyclosyn(z_diss_cm, z_max_cm)))
        self._lines["postcom"].set_data(
            frequency_hz, E * np.array(bj.get_observed_photon_integrated_flux_electron_compton(z_diss_cm, z_max_cm)))

        if bb is not None:
            if self._mode == "3ml":
                bb.update_observed_flux()
            self._lines["bb"].set_data(
                energy_eV_to_frequency_hz(np.array(bb.get_observed_energy()) * erg2eV),
                np.array(bb.get_observed_energy_flux()))
            self._lines["bb"].set_visible(True)
        else:
            self._lines["bb"].set_visible(False)

    def _draw_jet_shape(self):
        bl   = self._active_bljet
        logz = np.log10(np.asarray(bl.get_z_center_grid(), dtype=float) / bl.r_g)
        logr = np.log10(np.asarray(bl.get_radius_grid(), dtype=float)   / bl.r_g)
        self._jet_shape_line.set_data(logr, logz)

        self.cax.set_xlim(max(logr) + 0.5, min(logr) - 0.5)
        self.cax.set_ylim(logz[0], logz[-1])
        self.cax.set_xticks(np.arange(int(min(logr)), int(max(logr)) + 1))
        self.cax.set_yticks(np.arange(int(logz[0]),   int(logz[-1]) + 1))

        self._zdiss_line.set_ydata([np.log10(bl.z_dissipation)] * 2)
        self._zacc_line.set_ydata( [np.log10(bl.z_end_of_acceleration)] * 2)
        self._zmax_line.set_ydata( [np.log10(bl.z_max_calculation)] * 2)

    # ── zone inspector ────────────────────────────────────────────────────────

    def _build_zone_figure(self):
        plt.ioff()
        self._fig_zone = plt.figure(figsize=(13, 10))
        gs = gridspec.GridSpec(2, 2, hspace=0.38, wspace=0.32,
                               figure=self._fig_zone)
        self._ax_target = self._fig_zone.add_subplot(gs[0, 0])
        self._ax_cool   = self._fig_zone.add_subplot(gs[0, 1])
        self._ax_part   = self._fig_zone.add_subplot(gs[1, 0])
        self._ax_emit   = self._fig_zone.add_subplot(gs[1, 1])

        nan1 = np.array([np.nan])

        ax = self._ax_target
        self._zt_total, = ax.loglog(nan1, nan1, c="k",          lw=2,    label="total")
        self._zt_syn,   = ax.loglog(nan1, nan1, c="tab:blue",   ls="--", lw=1.5, label="cyclosyn.")
        self._zt_bb,    = ax.loglog(nan1, nan1, c="tab:purple", ls="--", lw=1.5, label="BB target")
        self._zt_bb.set_visible(False)
        ax.set_xlabel("Comoving energy [eV]")
        ax.set_ylabel(r"$E^2\,dN/dE$ [erg cm$^{-3}$]")
        ax.set_title("Target photon field (comoving)")
        ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=9)

        ax = self._ax_cool
        self._zc_momspan = ax.axvspan(1, 2, color="grey", alpha=0.2, label="zone mom. grid")
        self._zc_acc,   = ax.loglog(nan1, nan1, c="tab:blue",   lw=1.5, label="acceleration")
        self._zc_adi,   = ax.loglog(nan1, nan1, c="tab:orange", lw=1.5, label="adiabatic")
        self._zc_syn,   = ax.loglog(nan1, nan1, c="tab:green",  lw=1.5, label="cyclosyn.")
        self._zc_comKN, = ax.loglog(nan1, nan1, c="tab:red",    lw=1.5, label="Compton (KN)")
        self._zc_comTh, = ax.loglog(nan1, nan1, c="tab:purple", lw=1.5, ls="--", label="Compton (Thomson)")
        self._zc_vth  = ax.axvline(np.nan, ls=":", c="k",         lw=1.2, label=r"$p_\mathrm{thermal}$")
        self._zc_vbr  = ax.axvline(np.nan, ls=":", c="tab:red",   lw=1.2, label=r"$p_\mathrm{break}$")
        self._zc_vmax = ax.axvline(np.nan, ls=":", c="tab:green", lw=1.2, label=r"$p_\mathrm{max}$")
        ax.set_xlabel("p [eV/c]"); ax.set_ylabel("t [s]")
        ax.set_title("Cooling timescales (comoving)")
        ax.set_aspect("equal"); ax.grid(alpha=0.3); 
        ax.legend(fontsize=9, framealpha=1,
          loc="upper left", bbox_to_anchor=(1.01, 1), borderaxespad=0)

        ax = self._ax_part
        self._zp_elec, = ax.loglog(nan1, nan1, c="tab:blue", lw=1.5)
        self._zp_vth  = ax.axvline(np.nan, ls=":", c="k",         lw=1.2, label=r"$p_\mathrm{thermal}$")
        self._zp_vbr  = ax.axvline(np.nan, ls=":", c="tab:red",   lw=1.2, label=r"$p_\mathrm{break}$")
        self._zp_vmax = ax.axvline(np.nan, ls=":", c="tab:green", lw=1.2, label=r"$p_\mathrm{max}$")
        ax.set_xlabel("p [eV/c]"); ax.set_ylabel(r"$p\,dN/dp$ [cm$^{-3}$]")
        ax.set_title("Electron distribution (comoving)")
        ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=9)

        ax = self._ax_emit
        self._ze_total, = ax.loglog(nan1, nan1, c="k",          lw=2,    label="total")
        self._ze_syn,   = ax.loglog(nan1, nan1, c="tab:blue",   ls="--", lw=1.5, label="cyclosyn.")
        self._ze_com,   = ax.loglog(nan1, nan1, c="tab:red",    ls="--", lw=1.5, label="Compton")
        self._ze_bb,    = ax.loglog(nan1, nan1, c="tab:purple", ls="-.", lw=1.5, label="BB (obs.)")
        self._ze_bb.set_visible(False)
        self._ze_vic_syn, = ax.plot([np.nan], [1e-20], ls=":", c="tab:blue",   lw=1)
        self._ze_vic_bb,  = ax.plot([np.nan], [1e-20], ls=":", c="tab:purple", lw=1)
        self._ze_vic_bb.set_visible(False)
        ax.set_xlabel("Observed frequency [Hz]")
        ax.set_ylabel(r"$\nu F_\nu$ [erg cm$^{-2}$ s$^{-1}$]")
        ax.set_title("Observed emission (this zone)")
        ax.set_aspect("equal"); ax.grid(alpha=0.3); ax.legend(fontsize=9)

        # self._zone_suptitle = self._fig_zone.suptitle("", fontsize=11, y=1.01)
        plt.ion()

    def _draw_zone(self, izone):
        bl  = self._active_bljet
        bj  = self._active_bhjet
        bb  = self._active_bb
        rad = bj.radiation_zones[izone]

        logz = np.log10(np.asarray(bl.get_z_center_grid(), dtype=float) / bl.r_g)
        logr = np.log10(np.asarray(bl.get_radius_grid(), dtype=float)   / bl.r_g)
        self._zone_marker.set_data([logr[izone]], [logz[izone]])

        z_rg = bl.get_z_center_grid()[izone] / bl.r_g
        # self._zone_suptitle.set_text(
        #     f"Zone {izone}  |  z = {z_rg:.2e} rg  |  "
        #     f"B = {rad.magnetic_field:.2e} G  |  "
        #     f"T_e = {rad.electron_temperature:.2e} keV  |  "
        #     f"δ = {rad.doppler_factor_bulk:.2f}  |  "
        #     f"Compton: {'yes' if rad.compton_calculation_necessary() else 'no'}  |  "
        #     f"mode: {self._mode}"
        # )

        p_th  = rad.electron_temperature * 1e3
        p_br  = rad.get_electron_break_momentum() * erg2eV * c_cm_s
        p_max = rad.get_electron_max_momentum()   * erg2eV * c_cm_s

        # panel 1: target field
        E_tgt = np.array(rad.get_total_photon_target_energy())
        u_tgt = np.array(rad.get_total_photon_target_energy_density())
        self._zt_total.set_data(E_tgt * erg2eV, u_tgt)

        E_s, u_s = np.array([np.nan]), np.array([np.nan])
        try:
            E_s = np.array(rad.get_photon_target_energy("cyclosynchrotron"))
            u_s = np.array(rad.get_photon_target_energy_density("cyclosynchrotron"))
            self._zt_syn.set_data(E_s * erg2eV, u_s)
        except Exception:
            self._zt_syn.set_data([np.nan], [np.nan])

        if bb is not None:
            try:
                E_b = np.array(rad.get_photon_target_energy("bb_target"))
                u_b = np.array(rad.get_photon_target_energy_density("bb_target"))
                self._zt_bb.set_data(E_b * erg2eV, u_b)
                self._zt_bb.set_visible(True)
            except Exception:
                self._zt_bb.set_visible(False)
        else:
            self._zt_bb.set_visible(False)

        # panel 2: cooling timescales
        mom   = _MOM_GRID
        p_eVc = mom * erg2eV * c_cm_s
        self._zc_acc.set_data(  p_eVc, rad.get_timescale_electron_acceleration(mom))
        self._zc_adi.set_data(  p_eVc, rad.get_timescale_electron_adiabatic(mom))
        self._zc_syn.set_data(  p_eVc, rad.get_timescale_electron_cyclosyn(mom))
        self._zc_comKN.set_data(p_eVc, rad.get_timescale_electron_compton(mom))
        self._zc_comTh.set_data(p_eVc, rad.get_timescale_electron_compton_thomson(mom))
        for vline, val in [(self._zc_vth, p_th), (self._zc_vbr, p_br), (self._zc_vmax, p_max)]:
            vline.set_xdata([val, val])

        mom_zone = np.array(rad.get_electron_momentum_grid())
        p_lo = mom_zone[0]  * erg2eV * c_cm_s
        p_hi = mom_zone[-1] * erg2eV * c_cm_s
        self._zc_momspan.set_xy([
            [p_lo, 0], [p_lo, 1], [p_hi, 1], [p_hi, 0], [p_lo, 0],
        ])

        # panel 3: electron distribution
        dn_dp = np.array(rad.get_electron_momentum_number_density())
        self._zp_elec.set_data(mom_zone * erg2eV * c_cm_s, dn_dp)
        for vline, val in [(self._zp_vth, p_th), (self._zp_vbr, p_br), (self._zp_vmax, p_max)]:
            vline.set_xdata([val, val])

        # panel 4: observed emission
        E_tot    = np.array(rad.get_observed_photon_energy_grid_total())
        tot_flux = E_tot * np.array(rad.get_observed_photon_flux_total())
        self._ze_total.set_data(energy_eV_to_frequency_hz(E_tot * erg2eV), tot_flux)

        E_syn = np.array(rad.get_observed_photon_energy_grid_electron_cyclosyn())
        self._ze_syn.set_data(energy_eV_to_frequency_hz(E_syn * erg2eV),
            E_syn * np.array(rad.get_observed_photon_flux_electron_cyclosyn()))

        E_com = np.array(rad.get_observed_photon_energy_grid_electron_compton())
        self._ze_com.set_data(energy_eV_to_frequency_hz(E_com * erg2eV),
            E_com * np.array(rad.get_observed_photon_flux_electron_compton()))

        gamma_grid = np.array(rad.get_electron_gamma_grid())
        main_lf    = gamma_grid[np.argmax(mom_zone * dn_dp)] if len(dn_dp) > 0 else 1.0

        valid_s = np.isfinite(E_s) & np.isfinite(u_s) & (u_s > 0)
        if valid_s.any():
            E_syn_peak = np.exp(np.average(np.log(E_s[valid_s] * erg2eV), weights=u_s[valid_s]))
            self._ze_vic_syn.set_xdata([
                energy_eV_to_frequency_hz(main_lf**2 * 4 * E_syn_peak * rad.doppler_factor_bulk)
            ] * 2)
        else:
            self._ze_vic_syn.set_xdata([np.nan, np.nan])

        if bb is not None:
            if self._mode == "3ml":
                bb.update_observed_flux()
            bb_T_eV = bb.temperature * 1e3
            self._ze_vic_bb.set_xdata([
                energy_eV_to_frequency_hz(main_lf**2 * 4 * bb_T_eV * rad.doppler_factor_bulk)
            ] * 2)
            self._ze_vic_bb.set_visible(True)
            self._ze_bb.set_data(energy_eV_to_frequency_hz(np.array(bb.get_observed_energy()) * erg2eV),
                                 np.array(bb.get_observed_energy_flux()))
            self._ze_bb.set_visible(True)
        else:
            self._ze_vic_bb.set_visible(False)
            self._ze_bb.set_visible(False)

        for ax, ref_line in [
            (self._ax_target, self._zt_total),
            (self._ax_cool,   self._zc_acc),
            (self._ax_part,   self._zp_elec),
            (self._ax_emit,   self._ze_total),
        ]:
            ax.relim(); ax.autoscale_view()
            yv    = np.asarray(ref_line.get_ydata()).ravel()
            valid = np.isfinite(yv) & (yv > 0)
            if valid.any():
                lp = np.log10(yv[valid].max())
                ax.set_ylim(10**(lp - 20), 10**(lp + 2))

        self._fig_zone.canvas.draw_idle()
        self.fig.canvas.draw_idle()

    # ── widget construction ───────────────────────────────────────────────────
    def _build_widgets(self):
        d = self._defaults
        _all_widgets = []   # (parameter name, widget) pairs for _flat_widgets

        def is_visible(desc):
            return VISIBLE_PARAMETERS is None or desc in VISIBLE_PARAMETERS

        def limits(desc, lo, hi, step):
            return PARAMETER_LIMITS.get(desc, (lo, hi, step))

        def display_label(desc):
            return PARAMETER_LABELS.get(desc, desc)

        def linear_value_label(desc, log_value):
            """Format a log10 slider's physical value and, where relevant, its unit."""
            return f"= {10**log_value:.2e}{PARAMETER_VALUE_UNITS.get(desc, '')}"

        def bold_label(desc):
            return widgets.HTML(
                value=f"<b>{escape(display_label(desc))}</b>",
                layout=widgets.Layout(width="250px"),
            )

        def widget_row(desc, slider, trailing=None):
            """Add a rich bold label when this parameter is emphasised."""
            if desc not in EMPHASISED_PARAMETERS:
                return HBox([slider, trailing]) if trailing is not None else slider

            slider.description = ""
            slider.style.description_width = "0px"
            slider.layout = widgets.Layout(width="230px")
            children = [bold_label(desc), slider]
            if trailing is not None:
                children.append(trailing)
            return HBox(children)

        def fls(desc, val, lo, hi, step=0.01, fmt=".2f"):
            """Log slider — returns HBox(slider, live 10^x label)."""
            if not is_visible(desc):
                return None
            lo, hi, step = limits(desc, lo, hi, step)
            s = widgets.FloatSlider(value=val, min=lo, max=hi, step=step,
                description=display_label(desc), readout_format=fmt,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="480px"))
            lbl = widgets.Label(value=linear_value_label(desc, val),
                                layout=widgets.Layout(width="155px"))
            s.observe(lambda c, l=lbl, d=desc: l.__setattr__(
                "value", linear_value_label(d, c["new"])
            ),
                    names="value")
            _all_widgets.append((desc, s))
            return widget_row(desc, s, lbl)

        def fs(desc, val, lo, hi, step=0.01, fmt=".3f"):
            """Plain slider — returned as-is."""
            if not is_visible(desc):
                return None
            lo, hi, step = limits(desc, lo, hi, step)
            s = widgets.FloatSlider(value=val, min=lo, max=hi, step=step,
                description=display_label(desc), readout_format=fmt,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="480px"))
            _all_widgets.append((desc, s))
            return widget_row(desc, s)
        
        def fs_log(desc, val, lo, hi, step=0.01, fmt=".1e"):
            """Plain slider — returned as-is."""
            if not is_visible(desc):
                return None
            lo, hi, step = limits(desc, lo, hi, step)
            s = widgets.FloatLogSlider(value=val, min=lo, max=hi, step=step,
                description=display_label(desc), readout_format=fmt,
                style={"description_width": "250px"},
                layout=widgets.Layout(width="480px"))
            _all_widgets.append((desc, s))
            return widget_row(desc, s)

        def cb(desc, val):
            if not is_visible(desc):
                return None
            s = widgets.Checkbox(value=val, description=display_label(desc))
            _all_widgets.append((desc, s))
            return widget_row(desc, s)

        param_groups = {
            "Source Parameters": [
                fls("lg_10_mass_bh",                     d["lg_10_mass_bh"],             0,  10),
                fls("lg_10_jet_power_eddington",         d["lg_10_jet_power_eddington"], -10,   2),
                fs( "theta_obs",                         d["theta_obs"],                 0,  89, step=1),
                fls("lg_10_distance",                    d["lg_10_distance"],            0,  10),
                fs( "redshift",                          d["redshift"],                  0,   1, step=0.0001, fmt=".5f"),
            ],
            "Jet Properties": [
                fls("lg_10_z_jet_launching",             d["lg_10_z_jet_launching"],     0,   5),
                fls("lg_10_r_initial",                   d["lg_10_r_initial"],           0,   3),
                fls("lg_10_z_end_of_acceleration",       d["lg_10_z_end_of_acceleration"], 0, 10),
                cb( "link_particle_to_bulk_acceleration", d["link_particle_to_bulk_acceleration"]),
                fls("lg_10_z_dissipation",               d["lg_10_z_dissipation"],       0,  10),
                fls("lg_10_sigma_final",                 d["lg_10_sigma_final"],        -4,   1),
                fs( "gamma_final",                       d["gamma_final"],             1, 50, step=1, fmt=".0f"),
                cb( "calc_pair_content_from_plasma_beta",d["calc_pair_content_from_plasma_beta"]),
                fs( "plasma_beta_jet_base",              d["plasma_beta_jet_base"],      0.01, 10, step=0.01, fmt=".2f"),
                fls("lg_10_electron_temperature_jet_base", d["lg_10_electron_temperature_jet_base"], 0, 4),
                fls("lg_10_z_max_calculation",           d["lg_10_z_max_calculation"],   4,  8),
                fs( "gamma_acceleration_exponent",       d["gamma_acceleration_exponent"], 1, 10, step=0.05),
                fls("lg_10_opening_angle_constant",      d["lg_10_opening_angle_constant"], -5, 5),
            ],
            "Electrons": [
                fls("lg_10_fraction_nonthermal_electrons", d["lg_10_fraction_nonthermal_electrons"], -2, np.log10(0.99)),
                fls("lg_10_factor_break_electrons",      d["lg_10_factor_break_electrons"], -2, 0),
                fls("lg_10_factor_max_energy_electrons", d["lg_10_factor_max_energy_electrons"], -9, -1),
                fs( "index_injected_electrons",          d["index_injected_electrons"],           0,  5, step=0.05),
            ],
            "Protons (Dummy)": [
                fls("lg_10_fraction_nonthermal_protons", d["lg_10_fraction_nonthermal_protons"], -10, 0),
                fls("lg_10_factor_break_protons",        d["lg_10_factor_break_protons"], -10, 10),
                fls("lg_10_factor_max_energy_protons",   d["lg_10_factor_max_energy_protons"], -10, 10),
                fs( "index_injected_protons",            d["index_injected_protons"],             0,  5, step=0.05),
            ],
            "BB Target": [
                cb( "bb_enable",                         d["bb_enable"]),
                cb( "bb_add_to_total",                   d["bb_add_to_total"]),
                fls("lg_10_bb_temperature",              d["lg_10_bb_temperature"], 0, 3),
                fls("lg_10_bb_energy_density",           d["lg_10_bb_energy_density"], -9, -1),
                fls("lg_10_bb_luminosity",               d["lg_10_bb_luminosity"], 38, 46),
            ],
            "Meta-parameters": [
                cb( "include_counterjet",                d["include_counterjet"]),
                fls("lg_10_compton_threshold",           d["lg_10_compton_threshold"], -10, 0),
                fs_log( "dlgz",                          d["dlgz"],                  -3, 0, step=0.01, fmt=".1e"),
            ],
        }

        tab = Tab()
        children, titles = [], []
        for title, wlist in param_groups.items():
            wlist = [widget for widget in wlist if widget is not None]
            if not wlist:
                continue
            children.append(VBox(wlist))
            titles.append(title)
        tab.children = children
        for i, t in enumerate(titles):
            tab.set_title(i, t)

        self._tab          = tab
        self._flat_widgets = dict(_all_widgets)

        if self._manual_update:
            self._calculate_button = widgets.Button(
                description="Calculate model", button_style="primary",
                tooltip="Calculate the spectrum using the current controls.",
            )
            self._reset_button = widgets.Button(
                description="Reset to starting model",
                tooltip="Restore the tutorial starting values and calculate again.",
            )
            self._status = widgets.Label(value="Adjust parameters, then calculate the model.")
            self._calculate_button.on_click(lambda _: self.calculate())
            self._reset_button.on_click(lambda _: self.reset())
            self._manual_controls = HBox([
                self._calculate_button, self._reset_button, self._status,
            ])
        else:
            self._manual_controls = None

        if _THREEML_AVAILABLE:
            self._mode_toggle = widgets.ToggleButtons(
                options     = [("Pure classes", "pure"), ("3ML plugin", "3ml")],
                value       = "pure",
                description = "Mode:",
                style       = {"description_width": "60px", "button_width": "120px"},
            )
            self._mode_toggle.observe(self._on_mode_change, names="value")
        else:
            self._mode_toggle = widgets.Label("(3ML not available — pure mode only)")

    # ── meta-parameter check ──────────────────────────────────────────────────

    def _meta_params_changed(self, p):
        return (
            p["include_counterjet"]  != self._defaults["include_counterjet"]
            or p["lg_10_compton_threshold"] != self._defaults["lg_10_compton_threshold"]
            or p["dlgz"]             != self._defaults["dlgz"]
        )

    # ── update callbacks ──────────────────────────────────────────────────────

    def _on_mode_change(self, change):
        self._mode = change["new"]
        self.calculate()

    def _parameters_from_widgets(self):
        """Read the controls and keep the linked acceleration distances in sync."""
        p = self._defaults.copy()
        p.update({k: w.value for k, w in self._flat_widgets.items()})

        particle_slider = self._flat_widgets.get("lg_10_z_dissipation")
        if p["link_particle_to_bulk_acceleration"]:
            p["lg_10_z_dissipation"] = p["lg_10_z_end_of_acceleration"]
            if particle_slider is not None:
                particle_slider.disabled = True
                if particle_slider.value != p["lg_10_z_dissipation"]:
                    particle_slider.value = p["lg_10_z_dissipation"]
        elif particle_slider is not None:
            particle_slider.disabled = False
        return p

    def _on_widget_change(self, change=None):
        if self._loading_parameters:
            return
        if self._manual_update:
            self._parameters_from_widgets()
            return
        self.calculate()

    def calculate(self):
        """Calculate and redraw the model using the current widget values."""
        p = self._parameters_from_widgets()

        if self._meta_params_changed(p):
            for k in p:
                if k in self._defaults:
                    self._defaults[k] = p[k]
            self._build_models()
        else:
            if self._mode == "pure":
                self._update_pure(p)
                self._run_pure()
            else:
                self._update_3ml(p)
                self._run_3ml()

        self._draw_spectra()
        self._draw_jet_shape()
        self.fig.canvas.draw_idle()

        if hasattr(self, "_fig_zone") and plt.fignum_exists(self._fig_zone.number):
            self._draw_zone(self._zone_slider.value)

        if self._manual_update:
            self._status.value = "Model calculated."

    def _update_pure(self, p):
        bl = self.bljet
        bj = self.bhjet

        bl.mass_bh                            = 10**p["lg_10_mass_bh"]
        bl.jet_power_eddington                = 10**p["lg_10_jet_power_eddington"]
        bl.z_jet_launching                    = 10**p["lg_10_z_jet_launching"]
        bl.r_initial                          = 10**p["lg_10_r_initial"]
        bl.z_end_of_acceleration              = 10**p["lg_10_z_end_of_acceleration"]
        bl.z_dissipation                      = 10**p["lg_10_z_dissipation"]
        bl.z_max_calculation                  = 10**p["lg_10_z_max_calculation"]
        bl.sigma_final                        = 10**p["lg_10_sigma_final"]
        bl.gamma_final                        = p["gamma_final"]
        bl.plasma_beta_jet_base               = p["plasma_beta_jet_base"]
        bl.calc_pair_content_from_plasma_beta = p["calc_pair_content_from_plasma_beta"]
        bl.electron_temperature_jet_base      = 10**p["lg_10_electron_temperature_jet_base"]
        bl.gamma_acceleration_exponent        = p["gamma_acceleration_exponent"]
        bl.opening_angle_constant             = 10**p["lg_10_opening_angle_constant"]
        bl.fraction_nonthermal_electrons      = 10**p["lg_10_fraction_nonthermal_electrons"]
        bl.fraction_nonthermal_protons        = 10**p["lg_10_fraction_nonthermal_protons"]
        bl.factor_break_electrons             = 10**p["lg_10_factor_break_electrons"]
        bl.factor_break_protons               = 10**p["lg_10_factor_break_protons"]
        bl.factor_max_energy_electrons        = 10**p["lg_10_factor_max_energy_electrons"]
        bl.factor_max_energy_protons          = 10**p["lg_10_factor_max_energy_protons"]
        bl.index_injected_electrons           = p["index_injected_electrons"]
        bl.index_injected_protons             = p["index_injected_protons"]

        bj.theta_obs = p["theta_obs"]
        bj.distance  = 10**p["lg_10_distance"]
        bj.redshift  = p["redshift"]

        if p["bb_enable"]:
            if self._bb_target is None:
                self._bb_target = TargetBlackBody(
                    name              = "bb_target",
                    distance          = bj.distance,
                    redshift          = bj.redshift,
                    luminosity        = 10**p["lg_10_bb_luminosity"],
                    temperature       = 10**p["lg_10_bb_temperature"],
                    energy_density    = 10**p["lg_10_bb_energy_density"],
                    add_to_total_flux = p["bb_add_to_total"],
                )
                bj.add_target_photon_field(self._bb_target)
            else:
                self._bb_target.temperature       = 10**p["lg_10_bb_temperature"]
                self._bb_target.luminosity        = 10**p["lg_10_bb_luminosity"]
                self._bb_target.energy_density    = 10**p["lg_10_bb_energy_density"]
                self._bb_target.distance          = bj.distance
                self._bb_target.redshift          = bj.redshift
                self._bb_target.add_to_total_flux = p["bb_add_to_total"]
        else:
            if self._bb_target is not None:
                bj.remove_target_photon_field(self._bb_target)
                self._bb_target = None

    def _update_3ml(self, p):
        for ui_name, plugin_name in _PLUGIN_PARAMETER_NAMES.items():
            if plugin_name in self._plugin.parameters:
                self._plugin.parameters[plugin_name].value = p[ui_name]
        if "lg_gamma_final" in self._plugin.parameters:
            self._plugin.parameters["lg_gamma_final"].value = np.log10(p["gamma_final"])

        if p["bb_enable"]:
            if self._bb_plugin_target is None:
                self._bb_plugin_target = BHJetTargetBB()
                self._bb_plugin_target._setup(target_name="bb_target")
                self._plugin.add_target(self._bb_plugin_target, name="bb_target")
                self._bb_added = True
            bb_pv = self._bb_plugin_target.parameters
            bb_pv["lg_temperature"].value    = p["lg_10_bb_temperature"]
            bb_pv["lg_luminosity"].value     = p["lg_10_bb_luminosity"]
            bb_pv["lg_energy_density"].value = p["lg_10_bb_energy_density"]
            bb_pv["lg_distance"].value       = p["lg_10_distance"]
            bb_pv["redshift"].value          = p["redshift"]
            self._bb_plugin_target.bb.add_to_total_flux = p["bb_add_to_total"]
        else:
            if self._bb_added:
                self._plugin.remove_target("bb_target")
                self._bb_plugin_target = None
                self._bb_added         = False

    # ── public API ────────────────────────────────────────────────────────────

    def set_param(self, name, value):
        """Set a widget value programmatically."""
        self._flat_widgets[name].value = value

    def get_parameters(self):
        """Return the current parameters, including the values set in the widgets."""
        return self._parameters_from_widgets().copy()

    def set_parameters(self, parameters):
        """Apply saved parameters to the widgets and recalculate the model."""
        rebuild_model = False
        self._loading_parameters = True
        try:
            for name, value in parameters.items():
                if name in self._flat_widgets:
                    self._flat_widgets[name].value = value
                elif name in self._defaults:
                    self._defaults[name] = value
                    rebuild_model = True
        finally:
            self._loading_parameters = False

        if rebuild_model:
            self._build_models()
        self.calculate()

    def reset(self):
        """Restore the supplied starting model and redraw it."""
        for name, widget in self._flat_widgets.items():
            if name in self._starting_defaults:
                widget.value = self._starting_defaults[name]
        self.calculate()

    def get_total_sed(self):
        """Return the current total model as energy [eV] and nu F_nu [cgs]."""
        bhjet = self._active_bhjet
        energy_erg = np.asarray(bhjet.get_observed_photon_energy_grid(), dtype=float)
        sed = energy_erg * np.asarray(bhjet.get_observed_photon_flux_total(), dtype=float)
        return energy_erg * erg2eV, sed

    def display(self):
        """Render the SED UI: figure canvas + mode toggle + tab panel."""
        children = [self.fig.canvas]
        if self._manual_controls is not None:
            children.append(self._manual_controls)
        children.extend([self._mode_toggle, self._tab])
        display(VBox(children))

    def display_zones(self):
        """
        Render the zone inspector: 4-panel figure + zone index slider.
        Works in both pure and 3ML mode; suptitle shows current mode.
        """
        if not hasattr(self, "_fig_zone"):
            self._build_zone_figure()

        n_zones = self._active_bljet.n_zones
        self._zone_slider = widgets.IntSlider(
            value=n_zones // 2, min=0, max=n_zones - 1, step=1,
            description="Zone index",
            style={"description_width": "100px"},
            layout=widgets.Layout(width="500px"),
            continuous_update=True,
        )
        self._zone_label = widgets.Label(value="")

        def _on_zone(change):
            izone = change["new"]
            rad   = self._active_bhjet.radiation_zones[izone]
            z_rg  = self._active_bljet.get_z_center_grid()[izone] / self._active_bljet.r_g
            self._zone_label.value = (
                f"z = {z_rg:.2e} rg  |  "
                f"B = {rad.magnetic_field:.2e} G  |  "
                f"T_e = {rad.electron_temperature:.2e} keV  |  "
                f"δ = {rad.doppler_factor_bulk:.2f}  |  "
                f"Compton: {'yes' if rad.compton_calculation_necessary() else 'no'}"
            )
            self._draw_zone(izone)

        self._zone_slider.observe(_on_zone, names="value")
        self._draw_zone(self._zone_slider.value)

        display(VBox([
            self._fig_zone.canvas,
            HBox([self._zone_slider, self._zone_label]),
        ]))
