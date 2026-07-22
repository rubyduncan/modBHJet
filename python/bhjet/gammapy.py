from functools import lru_cache
import astropy.units as u
import numpy as np
import time
from gammapy.modeling import Parameter
from gammapy.modeling.models import SPECTRAL_MODEL_REGISTRY, SpectralModel

from . import BHJet, BLJet


SOLVER_ENERGY_MIN = 1e-8 * u.eV
SOLVER_ENERGY_MAX = 1e12 * u.eV
SOLVER_ENERGY_DEX_STEP = 0.1

SOLVER_ENERGY_GRID = 10 ** np.arange(
    np.log10(SOLVER_ENERGY_MIN.to_value(u.eV)),
    np.log10(SOLVER_ENERGY_MAX.to_value(u.eV)) + SOLVER_ENERGY_DEX_STEP,
    SOLVER_ENERGY_DEX_STEP,
) * u.eV
SOLVER_ENERGY_GRID_ERG = SOLVER_ENERGY_GRID.to_value(u.erg)


class BHJetSpectralModel(SpectralModel):
    '''Spectral model for BHJet -  accepts photon energies as an Astropy object, 
        converts internally to erg, and returns differential photon flux, dN/dE - which is what gammapy expects 
    '''

    tag = ["BHJetSpectralModel"]

    #jet dynamics - bljet parameters 
    mass_bh = Parameter("mass_bh", 1e9, frozen=True)
    jet_power_eddington = Parameter("jet_power_eddington", 1e-5)
    z_jet_launching = Parameter("z_jet_launching", 2.0, frozen=True)
    r_initial = Parameter("r_initial", 3.0)
    z_end_of_acceleration = Parameter("z_end_of_acceleration", 1e5)
    z_dissipation = Parameter("z_dissipation", 1e2)
    z_max_calculation = Parameter("z_max_calculation", 1e6)
    sigma_final = Parameter("sigma_final", 1.0)
    gamma_final = Parameter("gamma_final", 15.0)
    electron_temperature_jet_base = Parameter("electron_temperature_jet_base",1e4)
    gamma_acceleration_exponent = Parameter("gamma_acceleration_exponent",0.5,frozen=True)
    opening_angle_constant = Parameter("opening_angle_constant", 0.15,frozen=True)
    fraction_nonthermal_electrons = Parameter(
        "fraction_nonthermal_electrons",
        0.1,
    )
    fraction_nonthermal_protons = Parameter(
        "fraction_nonthermal_protons",
        0.1,
    )
    factor_break_electrons = Parameter(
        "factor_break_electrons",
        1.0,
    )
    factor_break_protons = Parameter(
        "factor_break_protons",
        1.0,
    )
    factor_max_energy_electrons = Parameter(
        "factor_max_energy_electrons",
        1.0,
    )
    factor_max_energy_protons = Parameter(
        "factor_max_energy_protons",
        1.0,
    )
    index_injected_electrons = Parameter(
        "index_injected_electrons",
        1.5,
    )
    index_injected_protons = Parameter(
        "index_injected_protons",
        2.0,
    )
    plasma_beta_jet_base = Parameter(
        "plasma_beta_jet_base",
        0.2,
    )
    dlgz = Parameter("dlgz", 0.1, frozen=True)
    cutoff_type = Parameter("cutoff_type", 0, frozen=True)

    #bhjet parameters 
    theta_obs = Parameter("theta_obs", 17.0, frozen=True)
    distance = Parameter("distance", 16e3, frozen=True)
    redshift = Parameter("redshift", 0.00428, frozen=True)
    compton_threshold = Parameter("compton_threshold",1e-5,frozen=True)


    @staticmethod
    def _compute_jet(**parameters):
        
        # t0 = time.perf_counter()

        dynamics = build_bljet(
            mass_bh=parameters["mass_bh"],
            jet_power_eddington=parameters["jet_power_eddington"],
            z_jet_launching=parameters["z_jet_launching"],
            r_initial=parameters["r_initial"],
            z_end_of_acceleration=parameters["z_end_of_acceleration"],
            z_dissipation=parameters["z_dissipation"],
            z_max_calculation=parameters["z_max_calculation"],
            sigma_final=parameters["sigma_final"],
            gamma_final=parameters["gamma_final"],
            electron_temperature_jet_base=parameters["electron_temperature_jet_base"],
            gamma_acceleration_exponent=parameters["gamma_acceleration_exponent"],
            opening_angle_constant=parameters["opening_angle_constant"],
            fraction_nonthermal_electrons=parameters["fraction_nonthermal_electrons"],
            fraction_nonthermal_protons=parameters["fraction_nonthermal_protons"],
            factor_break_electrons=parameters["factor_break_electrons"],
            factor_break_protons=parameters["factor_break_protons"],
            factor_max_energy_electrons=parameters["factor_max_energy_electrons"],
            factor_max_energy_protons=parameters["factor_max_energy_protons"],
            index_injected_electrons=parameters["index_injected_electrons"],
            index_injected_protons=parameters["index_injected_protons"],
            plasma_beta_jet_base=parameters["plasma_beta_jet_base"],
            dlgz=parameters["dlgz"],
            cutoff_type=parameters["cutoff_type"],

        )

        # t1 = time.perf_counter()

        jet = BHJet(
            theta_obs=float(parameters["theta_obs"].value),
            distance=float(parameters["distance"].value),
            redshift=float(parameters["redshift"].value),
            include_counterjet=True,
            compton_threshold=float(parameters["compton_threshold"].value),
            profile_time=False,
            verbosity_level=0,
        )
        jet.init_jet_dynamics(dynamics)

        # t2 = time.perf_counter()

        # t3 = time.perf_counter()

        photon_energy_grid = [float(value) for value in SOLVER_ENERGY_GRID_ERG]
        jet.compute_full_jet(photon_energy_grid=photon_energy_grid)

        # t4 = time.perf_counter()
        
        # print(f"BLJet construction: {t1-t0:.4f} s")
        # print(f"BHJet construction: {t2-t1:.4f} s")
        # print(f"energy_shape:   {t3-t2:.4f} s")
        # print(f"compute_full_jet:  {t4-t3:.4f} s")
        
        return jet, dynamics


    def _solution_cache_key(self, **parameters):
        ''' Make a key for the fixed solver grid and current model parameters. '''

        energy_grid_key = (
            SOLVER_ENERGY_MIN.to_value(u.eV),
            SOLVER_ENERGY_MAX.to_value(u.eV),
            SOLVER_ENERGY_DEX_STEP,
        )
        parameter_key = tuple(
            (name, float(parameter.value))
            for name, parameter in sorted(parameters.items())
        )

        return energy_grid_key, parameter_key


    def _get_jet(self, **parameters):
        ''' Reuse the most recent BHJet solution when the inputs are unchanged. '''

        cache_key = self._solution_cache_key(**parameters)

        if cache_key == getattr(self, "_last_solution_key", None):
            return self._last_solution

        solution = self._compute_jet(**parameters)
        self._last_solution_key = cache_key
        self._last_solution = solution

        return solution


    @staticmethod
    def _interpolate_flux(energy, flux):
        ''' Interpolate a fixed-grid BHJet flux to requested energies. '''

        energy_shape = energy.shape
        requested_energy_erg = np.asarray(
            energy.to_value(u.erg),
            dtype=float,
        )
        requested_energy_erg = np.atleast_1d(requested_energy_erg).reshape(-1)
        flux = np.asarray(flux, dtype=float)

        interpolated_flux = np.zeros_like(requested_energy_erg)
        valid = requested_energy_erg > 0
        safe_flux = np.where(
            np.isfinite(flux) & (flux > 0),
            flux,
            1e-100,
        )

        interpolated_flux[valid] = np.exp(
            np.interp(
                np.log(requested_energy_erg[valid]),
                np.log(SOLVER_ENERGY_GRID_ERG),
                np.log(safe_flux),
                left=-1e100,
                right=-1e100,
            )
        )

        if energy_shape == ():
            return interpolated_flux[0]

        return interpolated_flux.reshape(energy_shape)


    def __getstate__(self):
        ''' Do not copy the transient C++ BHJet cache with a model. '''

        state = self.__dict__.copy()
        state.pop("_last_solution_key", None)
        state.pop("_last_solution", None)

        return state


    def evaluate(self,
        energy,
        mass_bh,
        jet_power_eddington,
        z_jet_launching,
        r_initial,
        z_end_of_acceleration,
        z_dissipation,
        z_max_calculation,
        sigma_final,
        gamma_final,
        electron_temperature_jet_base,
        gamma_acceleration_exponent,
        opening_angle_constant,
        fraction_nonthermal_electrons,
        fraction_nonthermal_protons,
        factor_break_electrons,
        factor_break_protons,
        factor_max_energy_electrons,
        factor_max_energy_protons,
        index_injected_electrons,
        index_injected_protons,
        plasma_beta_jet_base,
        dlgz,
        cutoff_type,
        theta_obs,
        distance,
        redshift,
        compton_threshold,
    ):
        
        parameters = locals()
        parameters.pop("self")
        parameters.pop("energy")
        jet, dynamics = self._get_jet(**parameters)

    #the photon total flux from bhjet is returned as E dN/dE in cm-2 s-1
        e_dnde = self._interpolate_flux(
            energy,
            jet.get_observed_photon_flux_total(),
        )
        dnde = e_dnde / energy.to_value(u.erg) #then this needs to be converted to number flux 
        number_flux = dnde * u.Unit("cm-2 s-1 erg-1")

        return number_flux


    def evaluate_components(self, energy):

        ''' Calculate the individual BHJet components as E dN/dE. '''

        parameters = {parameter.name: parameter for parameter in self.parameters} 
        jet, dynamics = self._get_jet(**parameters)

        z_dissipation_cm = dynamics.z_dissipation * dynamics.r_g
        z_max_calculation_cm = dynamics.z_max_calculation * dynamics.r_g
        flux_unit = u.Unit("cm-2 s-1")

        def get_flux(getter, *args):
            return self._interpolate_flux(energy, getter(*args)) * flux_unit

        outputs = {
            "pre_syn": get_flux(
                jet.get_observed_photon_integrated_flux_electron_cyclosyn,
                0.0,
                z_dissipation_cm,
            ),
            "pre_com": get_flux(
                jet.get_observed_photon_integrated_flux_electron_compton,
                0.0,
                z_dissipation_cm,
            ),
            "post_syn": get_flux(
                jet.get_observed_photon_integrated_flux_electron_cyclosyn,
                z_dissipation_cm,
                z_max_calculation_cm,
            ),
            "post_com": get_flux(
                jet.get_observed_photon_integrated_flux_electron_compton,
                z_dissipation_cm,
                z_max_calculation_cm,
            ),
            "total": get_flux(jet.get_observed_photon_flux_total),
        }

        return outputs


    def plot_components(
        self,
        energy_bounds,
        ax=None,
        n_points=200,
        energy_power=2,
        components=None,
    ):
        ''' Plot the individual BHJet components. '''

        energy_bounds = u.Quantity(energy_bounds)

        if energy_bounds.size != 2:
            raise ValueError("energy_bounds must contain exactly two energies.")

        energy_min, energy_max = energy_bounds
        energy = np.geomspace(
            energy_min.to_value(energy_bounds.unit),
            energy_max.to_value(energy_bounds.unit),
            n_points,
        ) * energy_bounds.unit

        outputs = self.evaluate_components(energy)

        styles = {
            "pre_syn": {"color": "tab:blue", "label": "pre syn"},
            "pre_com": {"color": "tab:green", "label": "pre com"},
            "post_syn": {"color": "tab:orange", "label": "post syn"},
            "post_com": {"color": "tab:red", "label": "post com"},
            "total": {"color": "k", "ls": "--", "label": "total"},
        }

        if components is None:
            components = styles.keys()

        if ax is None:
            import matplotlib.pyplot as plt

            _, ax = plt.subplots()

        for name in components:
            if name not in outputs:
                raise ValueError(f"Unknown BHJet component: {name}")

            flux = energy ** (energy_power - 1) * outputs[name]

            if energy_power == 0:
                flux = flux.to("cm-2 s-1 TeV-1")
            elif energy_power == 1:
                flux = flux.to("cm-2 s-1")
            elif energy_power == 2:
                flux = flux.to("erg cm-2 s-1")

            ax.loglog(
                energy.to_value(energy.unit),
                flux.to_value(flux.unit),
                **styles[name],
            )

        return ax
    

def parameter_value(parameter):
    return float(parameter.value)


def build_bljet(
    #this is to keep the parameters for the different jet dynamics classes separate
    mass_bh,
    jet_power_eddington,
    z_jet_launching,
    r_initial,
    z_end_of_acceleration,
    z_dissipation,
    z_max_calculation,
    sigma_final,
    gamma_final,
    electron_temperature_jet_base,
    gamma_acceleration_exponent,
    opening_angle_constant,
    fraction_nonthermal_electrons,
    fraction_nonthermal_protons,
    factor_break_electrons,
    factor_break_protons,
    factor_max_energy_electrons,
    factor_max_energy_protons,
    index_injected_electrons,
    index_injected_protons,
    plasma_beta_jet_base,
    dlgz,
    cutoff_type,
):
    return BLJet(
        mass_bh= parameter_value (mass_bh),
        jet_power_eddington=parameter_value(jet_power_eddington),
        z_jet_launching=parameter_value(z_jet_launching),
        r_initial=parameter_value(r_initial),
        z_end_of_acceleration=parameter_value(z_end_of_acceleration),
        z_dissipation=parameter_value(z_dissipation),
        z_max_calculation=parameter_value(z_max_calculation),
        sigma_final=parameter_value(sigma_final),
        gamma_final=parameter_value(gamma_final),
        electron_temperature_jet_base=parameter_value(
            electron_temperature_jet_base
        ),
        gamma_acceleration_exponent=parameter_value(
            gamma_acceleration_exponent
        ),
        opening_angle_constant=parameter_value(
            opening_angle_constant
        ),
        fraction_nonthermal_electrons=parameter_value(
            fraction_nonthermal_electrons
        ),
        fraction_nonthermal_protons=parameter_value(
            fraction_nonthermal_protons
        ),
        factor_break_electrons=parameter_value(
            factor_break_electrons
        ),
        factor_break_protons=parameter_value(
            factor_break_protons
        ),
        factor_max_energy_electrons=parameter_value(
            factor_max_energy_electrons
        ),
        factor_max_energy_protons=parameter_value(
            factor_max_energy_protons
        ),
        index_injected_electrons=parameter_value(
            index_injected_electrons
        ),
        index_injected_protons=parameter_value(
            index_injected_protons
        ),
        plasma_beta_jet_base=parameter_value(
            plasma_beta_jet_base
        ),
        dlgz=parameter_value(dlgz),
        cutoff_type=int(cutoff_type.value),
        calc_pair_content_from_plasma_beta=False,
    )


if BHJetSpectralModel not in SPECTRAL_MODEL_REGISTRY:
    SPECTRAL_MODEL_REGISTRY.append(BHJetSpectralModel)
