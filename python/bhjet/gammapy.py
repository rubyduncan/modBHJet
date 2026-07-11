import astropy.units as u
import numpy as np
import time

from gammapy.modeling import Parameter
from gammapy.modeling.models import SpectralModel

from . import BHJet, BLJet

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
    def evaluate(energy,
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
        
        # t0 = time.perf_counter()

        dynamics = build_bljet(
            mass_bh=mass_bh,
            jet_power_eddington=jet_power_eddington,
            z_jet_launching=z_jet_launching,
            r_initial=r_initial,
            z_end_of_acceleration=z_end_of_acceleration,
            z_dissipation=z_dissipation,
            z_max_calculation=z_max_calculation,
            sigma_final=sigma_final,
            gamma_final=gamma_final,
            electron_temperature_jet_base=electron_temperature_jet_base,
            gamma_acceleration_exponent=gamma_acceleration_exponent,
            opening_angle_constant=opening_angle_constant,
            fraction_nonthermal_electrons=fraction_nonthermal_electrons,
            fraction_nonthermal_protons=fraction_nonthermal_protons,
            factor_break_electrons=factor_break_electrons,
            factor_break_protons=factor_break_protons,
            factor_max_energy_electrons=factor_max_energy_electrons,
            factor_max_energy_protons=factor_max_energy_protons,
            index_injected_electrons=index_injected_electrons,
            index_injected_protons=index_injected_protons,
            plasma_beta_jet_base=plasma_beta_jet_base,
            dlgz=dlgz,
            cutoff_type=cutoff_type,

        )

        # t1 = time.perf_counter()

        jet = BHJet(
            theta_obs=float(theta_obs.value),
            distance=float(distance.value),
            redshift=float(redshift.value),
            include_counterjet=True,
            compton_threshold=float(compton_threshold.value),
            profile_time=False,
            verbosity_level=0,
        )

        jet.init_jet_dynamics(dynamics)

        # t2 = time.perf_counter()

    # bhjet is expecting an energy grid in ergs
        energy_shape = energy.shape
        energy_erg = np.atleast_1d(np.asarray(energy.to_value(u.erg), dtype=float))

    # #will compute photon energies also using the erg grid 
        jet.compute_full_jet(photon_energy_grid=energy_erg.tolist())

        # t3 = time.perf_counter()

    #the photon total flux from bhjet is returned as E dN/dE in cm-2 s-1
        e_dnde = np.asarray(jet.get_observed_photon_flux_total(),dtype=float)
        dnde = e_dnde / energy_erg #then this needs to be converted to number flux 
        number_flux = dnde * u.Unit("cm-2 s-1 erg-1")

        # t4 = time.perf_counter()

        if energy_shape == ():
            return number_flux[0]
        
        # print(f"BLJet construction: {t1-t0:.4f} s")
        # print(f"BHJet construction: {t2-t1:.4f} s")
        # print(f"compute_full_jet:   {t3-t2:.4f} s")
        # print(f"Output conversion:  {t4-t3:.4f} s")
        
        return number_flux
    

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