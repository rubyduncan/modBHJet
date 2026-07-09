# <litgen_stub>
####################    <generated_from:JetDynamics.hpp>    ####################





class JetDynamics:
    # Member variables
    n_zones: int = defaults.DEFAULT_N_ZONES
    verbosity_level: int = defaults.DEFAULT_VERBOSITY_LEVEL

    def __init__(
        self,
        n_zones_: int = defaults.DEFAULT_N_ZONES,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass




    def compute_jet_dynamics(self) -> None:
        pass
    def info(self) -> str:
        pass

    # readout functions
    def get_z_min_grid(self) -> List[float]:
        pass
    def get_z_height_grid(self) -> List[float]:
        pass
    def get_z_center_grid(self) -> List[float]:
        pass
    def get_radius_grid(self) -> List[float]:
        pass
    def get_gamma_grid(self) -> List[float]:
        pass
    def get_beta_grid(self) -> List[float]:
        pass
    def get_beta_gamma_grid(self) -> List[float]:
        pass
    def get_magnetic_field_grid(self) -> List[float]:
        pass
    def get_electron_density_grid(self) -> List[float]:
        pass
    def get_proton_density_grid(self) -> List[float]:
        pass
    def get_electron_temperature_grid(self) -> List[float]:
        pass
    def get_proton_temperature_grid(self) -> List[float]:
        pass
    def get_fraction_nonthermal_electrons_grid(self) -> List[float]:
        pass
    def get_fraction_nonthermal_protons_grid(self) -> List[float]:
        pass
    def get_factor_break_electrons_grid(self) -> List[float]:
        pass
    def get_factor_break_protons_grid(self) -> List[float]:
        pass
    def get_factor_max_energy_electrons_grid(self) -> List[float]:
        pass
    def get_factor_max_energy_protons_grid(self) -> List[float]:
        pass
    def get_index_injected_electrons_grid(self) -> List[float]:
        pass
    def get_index_injected_protons_grid(self) -> List[float]:
        pass
    def get_cutoff_type_grid(self) -> List[int]:
        pass


####################    </generated_from:JetDynamics.hpp>    ####################


####################    <generated_from:TargetPhotonField.hpp>    ####################


class TargetPhotonField:

    # Member variables
    name: str
    target_type: str
    distance: float
    redshift: float
    add_to_total_flux: bool
    verbosity_level: int

    def __init__(
        self,
        name_: str,
        distance_: float = defaults.DISTANCE,
        redshift_: float = defaults.REDSHIFT,
        add_to_total_flux_: bool = defaults.DEFAULT_ADD_TO_TOTAL_FLUX,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass

    def get_target_energy_grid_and_density(
        self,
        z: float,
        bulk_momentum: float,
        theta_obs: float,
        min_energy: float,
        max_energy: float
        ) -> Tuple[List[float], List[float]]:
        pass


    # readout functions for total flux
    def update_observed_flux(self) -> None:
        pass
    def get_observed_energy(self) -> List[float]:
        pass
    def get_observed_energy_flux(self) -> List[float]:
        pass
    def get_observed_number_flux(self) -> List[float]:
        pass



####################    </generated_from:TargetPhotonField.hpp>    ####################


####################    <generated_from:TargetBlackBody.hpp>    ####################


class TargetBlackBody:
    # Member variables
    luminosity: float = defaults.DEFAULT_LUMINOSTIY
    temperature: float = defaults.DEFAULT_TEMPERATURE
    energy_density: float = defaults.DEFAULT_ENERGY_DENSITY

    def __init__(
        self,
        name_: str,
        distance_: float = defaults.DISTANCE,
        redshift_: float = defaults.REDSHIFT,
        luminosity_: float = defaults.DEFAULT_LUMINOSTIY,
        temperature_: float = defaults.DEFAULT_TEMPERATURE,
        energy_density_: float = defaults.DEFAULT_ENERGY_DENSITY,
        add_to_total_flux_: bool = defaults.DEFAULT_ADD_TO_TOTAL_FLUX,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass

    def get_target_energy_grid_and_density(
        self,
        z: float,
        bulk_momentum: float,
        theta_obs: float,
        min_energy: float,
        max_energy: float
        ) -> Tuple[List[float], List[float]]:
        pass

    def update_observed_flux(self) -> None:
        pass



####################    </generated_from:TargetBlackBody.hpp>    ####################


####################    <generated_from:RadiationZone.hpp>    ####################



# struct TargetBlackBody
# {
#     double temperature, energy_density;
#     std::string name;
#     TargetBlackBody(double t, double u, std::string n)
#         : temperature(t), energy_density(u), name(n) {}
# };
# struct TargetDisk
# {
#     double mass_bh, inner_radius, outer_radius, luminosity, inclination;
# };

class RadiationZone:

    # total integral (for Thomson cooling)
    radiation_energy_density: float


    def add_target_photon_field(
        self,
        energies: List[float],
        energy_densities: List[float],
        name: str
        ) -> None:
        pass
    def remove_target_photon_field(self, name: str) -> None:
        pass
    # radiation targets per zone, for total luminosity track one for each in BHJet class
    # None add_target_black_body(double temperature, double energy_density, std::string name);
    # None add_target_disk(double Mbh, double inner_radius, double outer_radius, double luminosity, double inclination);
    # None add_target_field(std::vector<double> target_energy, std::vector<double> target_array);

    def compute_particles(self, reset_photon_targets: bool = True) -> None:
        pass

    @overload
    def compute_radiation(self) -> None:
        pass
    @overload
    def compute_radiation(self, obs_energy_grid: List[float]) -> None:
        pass


    # ----------------------------
    # Member variables
    # ----------------------------
    magnetic_field: float
    radius: float
    height: float
    bulk_momentum: float
    theta_obs: float
    distance: float
    redshift: float
    electron_number_density: float
    proton_number_density: float
    electron_temperature: float
    proton_temperature: float
    fraction_nonthermal_electrons: float
    fraction_nonthermal_protons: float
    factor_break_electrons: float
    factor_break_protons: float
    factor_max_energy_electrons: float
    factor_max_energy_protons: float
    index_injected_electrons: float
    index_injected_protons: float
    cutoff_type: int
    geometry: str
    include_counterjet: bool = defaults.DEFAULT_INCLUDE_COUNTERJET
    force_compton_calculation: bool = defaults.DEFAULT_FORCE_COMPTON_CALCULATION
    compton_switch: bool = defaults.DEFAULT_COMPTON_SWITCH
    compton_threshold: float = defaults.DEFAULT_COMPTON_THRESHOLD
    profile_time: bool = defaults.DEFAULT_PROFILE_TIME
    verbosity_level: int = defaults.DEFAULT_VERBOSITY_LEVEL

    # momentum grid
    n_bins_e: int
    # momentum grid
    n_bins_p: int




    syn_res: int = 10
    com_res: int = 6

    doppler_factor_bulk: float
    beta_bulk: float
    gamma_bulk: float

    # std::vector<TargetBlackBody> target_vector_blackbody;
    # std::vector<TargetDisk> target_vector_disk;
    # std::vector<TargetPhotonField> target_vector_photons;

    computation_times: List[float]

    def __init__(
        self,
        magnetic_field_: float = defaults.DEFAULT_MAGNETIC_FIELD,
        radius_: float = defaults.DEFAULT_RADIUS,
        height_: float = defaults.DEFAULT_HEIGHT,
        geometry_: str = defaults.DEFAULT_GEOMETRY,
        bulk_momentum_: float = defaults.DEFAULT_BULK_MOMENTUM,
        theta_obs_: float = defaults.DEFAULT_THETA_OBS,
        distance_: float = defaults.DEFAULT_DISTANCE,
        redshift_: float = defaults.DEFAULT_REDSHIFT,
        electron_number_density_: float = defaults.DEFAULT_ELECTRON_NUMBER_DENSITY,
        proton_number_density_: float = defaults.DEFAULT_PROTON_NUMBER_DENSITY,
        electron_temperature_: float = defaults.DEFAULT_ELECTRON_TEMPERATURE,
        proton_temperature_: float = defaults.DEFAULT_PROTON_TEMPERATURE,
        fraction_nonthermal_electrons_: float = defaults.DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
        fraction_nonthermal_protons_: float = defaults.DEFAULT_FRACTION_NONTHERMAL_PROTONS,
        factor_break_electrons_: float = defaults.DEFAULT_FACTOR_BREAK_ELECTRONS,
        factor_break_protons_: float = defaults.DEFAULT_FACTOR_BREAK_PROTONS,
        factor_max_energy_electrons_: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
        factor_max_energy_protons_: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
        index_injected_electrons_: float = defaults.DEFAULT_INDEX_INJECTED_ELECTRONS,
        index_injected_protons_: float = defaults.DEFAULT_INDEX_INJECTED_PROTONS,
        include_counterjet_: bool = defaults.DEFAULT_INCLUDE_COUNTERJET,
        force_compton_calculation_: bool = defaults.DEFAULT_FORCE_COMPTON_CALCULATION,
        compton_switch_: bool = defaults.DEFAULT_COMPTON_SWITCH,
        compton_threshold_: float = defaults.DEFAULT_COMPTON_THRESHOLD,
        profile_time_: bool = defaults.DEFAULT_PROFILE_TIME,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL,
        cutoff_type_: int = defaults.DEFAULT_CUTOFF_TYPE
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass

    def get_timescale_electron_cyclosyn(self, momentum: List[float]) -> List[float]:
        pass
    def get_timescale_electron_adiabatic(self, momentum: List[float]) -> List[float]:
        pass
    def get_timescale_electron_compton_thomson(self, momentum: List[float]) -> List[float]:
        pass
    def get_timescale_electron_compton(self, momentum: List[float]) -> List[float]:
        pass
    def get_timescale_electron_acceleration(self, momentum: List[float]) -> List[float]:
        pass

    def get_timescale_photon_cyclosyn_selfabsorption(
        self,
        momentum: List[float]
        ) -> List[float]:
        pass
    def get_timescale_photon_escape(self, momentum: List[float]) -> List[float]:
        """ read out the array, left for debugging
         std::vector<double> get_timescale_photon_cyclosyn_selfabsorption_array(){return cyclosyn_selfabsorption_rate;};
         std::vector<double> get_timescale_photon_cyclosyn_selfabsorption_energy(){return cyclosyn_energy;};
        """
        pass


    def get_electron_max_momentum(self) -> float:
        pass
    def get_electron_break_momentum(self) -> float:
        pass

    def get_electron_momentum_grid(self) -> List[float]:
        pass
    def get_electron_gamma_grid(self) -> List[float]:
        pass
    def get_electron_momentum_number_density(self) -> List[float]:
        pass
    def get_electron_gamma_number_density(self) -> List[float]:
        pass

    # std::vector<double> additional_target_field_energy, additional_target_field_energy_density;
    def get_total_photon_target_energy(self) -> List[float]:
        pass
    def get_total_photon_target_energy_density(self) -> List[float]:
        pass
    def get_photon_target_energy(self, name: str) -> List[float]:
        pass
    def get_photon_target_energy_density(self, name: str) -> List[float]:
        pass
    # std::vector<double> get_photon_target_energy_density_black_body(std::string name);
    # double get_target_black_body_temperature(std::string name);
    # double get_target_black_body_energy_density(std::string name);
    # None set_target_black_body_temperature(std::string name, double new_temperature);
    # None set_target_black_body_energy_density(std::string name, double new_energy_density);


    def get_observed_photon_energy_grid_electron_cyclosyn(self) -> List[float]:
        pass
    def get_observed_photon_luminosity_electron_cyclosyn(self) -> List[float]:
        pass
    def get_observed_photon_flux_electron_cyclosyn(self) -> List[float]:
        pass
    def get_observed_photon_energy_grid_electron_compton(self) -> List[float]:
        pass
    def get_observed_photon_luminosity_electron_compton(self) -> List[float]:
        pass
    def get_observed_photon_flux_electron_compton(self) -> List[float]:
        pass
    def get_observed_photon_energy_grid_total(self) -> List[float]:
        pass
    def get_observed_photon_luminosity_total(self) -> List[float]:
        pass
    def get_observed_photon_flux_total(self) -> List[float]:
        pass

    def get_computation_times(self) -> List[float]:
        pass

    def sum_jet_and_counterjet(
        self,
        size: int,
        input_en: List[float],
        input_lum: List[float],
        en: List[float],
        lum: List[float]
        ) -> None:
        pass
    def sum_jet_only(
        self,
        size: int,
        input_en: List[float],
        input_lum: List[float],
        en: List[float],
        lum: List[float]
        ) -> None:
        pass

    def compton_calculation_necessary(self) -> bool:
        pass
    def fic(self, gamma: float, alpha: float) -> float:
        """*
                 * @brief function to compute ker_ic[][].
                 *
                 * cf. Jones, "Inverse Compton Scattering of Cosmic-Ray Electrons", P.R. 1965, Eqn. 13-15
                 *
                 * @param gamma Energy of incoming particle
                 * @param alpha Energy of incoming photon / me c^2

        """
        pass
    # utility functions to construct Fic
    def fic_1(self, z: float) -> float:
        """< utility functions to construct Fic"""
        pass
    def fic_2(self, z: float) -> float:
        """< utility functions to construct Fic"""
        pass
    def polylog(self, z: float) -> float:
        """ utility functions to construct Fic, valid for -inf<z<1.0"""
        pass
    def polylog_smallz(self, z: float) -> float:
        """< utility functions to calc polylog"""
        pass
    def polylog_rk4(self, z: float) -> float:
        """< utility functions to construct polylog"""
        pass
    def polylog_largez(self, z: float) -> float:
        """< utility functions to construct polylog"""
        pass


####################    </generated_from:RadiationZone.hpp>    ####################


####################    <generated_from:BLJet.hpp>    ####################



class BLJet:
    # Member variables
    mass_bh: float = defaults.DEFAULT_MASS_BH
    jet_power_eddington: float = defaults.DEFAULT_JET_POWER_EDDINGTON
    z_jet_launching: float = defaults.DEFAULT_Z_JET_LAUNCHING
    r_initial: float = defaults.DEFAULT_R_INITIAL
    z_end_of_acceleration: float = defaults.DEFAULT_Z_END_OF_ACCELERATION
    z_dissipation: float = defaults.DEFAULT_Z_DISSIPATION
    z_max_calculation: float = defaults.DEFAULT_Z_MAX_CALCULATION
    sigma_final: float = defaults.DEFAULT_SIGMA_FINAL
    gamma_final: float = defaults.DEFAULT_GAMMA_FINAL
    electron_temperature_jet_base: float = defaults.DEFAULT_ELECTRON_TEMPERATURE_JET_BASE
    gamma_acceleration_exponent: float = defaults.DEFAULT_GAMMA_ACCELERATION_EXPONENT
    gamma_deceleration_exponent: float = defaults.DEFAULT_GAMMA_DECELERATION_EXPONENT
    opening_angle_constant: float = defaults.DEFAULT_OPENING_ANGLE_CONSTANT
    fraction_nonthermal_electrons: float = defaults.DEFAULT_FRACTION_NONTHERMAL_ELECTRONS
    fraction_nonthermal_protons: float = defaults.DEFAULT_FRACTION_NONTHERMAL_PROTONS
    factor_break_electrons: float = defaults.DEFAULT_FACTOR_BREAK_ELECTRONS
    factor_break_protons: float = defaults.DEFAULT_FACTOR_BREAK_PROTONS
    factor_max_energy_electrons: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS
    factor_max_energy_protons: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_PROTONS
    index_injected_electrons: float = defaults.DEFAULT_INDEX_INJECTED_ELECTRONS
    index_injected_protons: float = defaults.DEFAULT_INDEX_INJECTED_PROTONS
    cutoff_type: int = defaults.DEFAULT_CUTOFF_TYPE
    calc_pair_content_from_plasma_beta: bool = defaults.DEFAULT_CALC_PAIR_CONTENT_FROM_PLASMA_BETA
    plasma_beta_jet_base: float = defaults.DEFAULT_PLASMA_BETA_JET_BASE
    dlgz: float = defaults.DEFAULT_DLGZ

    # internal variables
    eddington_luminosity: float
    # internal variables
    r_g: float
    # internal variables
    zmin: float

    # hardcoded values
    # size_t n_bins_speed = 54;
    # double jet_opening_constant = 0.15;

    def __init__(
        self,
        mass_bh_: float = defaults.DEFAULT_MASS_BH,
        jet_power_eddington_: float = defaults.DEFAULT_JET_POWER_EDDINGTON,
        z_jet_launching_: float = defaults.DEFAULT_Z_JET_LAUNCHING,
        r_initial_: float = defaults.DEFAULT_R_INITIAL,
        z_end_of_acceleration_: float = defaults.DEFAULT_Z_END_OF_ACCELERATION,
        z_dissipation_: float = defaults.DEFAULT_Z_DISSIPATION,
        z_max_calculation_: float = defaults.DEFAULT_Z_MAX_CALCULATION,
        sigma_final_: float = defaults.DEFAULT_SIGMA_FINAL,
        gamma_final_: float = defaults.DEFAULT_GAMMA_FINAL,
        electron_temperature_jet_base_: float = defaults.DEFAULT_ELECTRON_TEMPERATURE_JET_BASE,
        gamma_acceleration_exponent_: float = defaults.DEFAULT_GAMMA_ACCELERATION_EXPONENT,
        gamma_deceleration_exponent_: float = defaults.DEFAULT_GAMMA_DECELERATION_EXPONENT,
        opening_angle_constant_: float = defaults.DEFAULT_OPENING_ANGLE_CONSTANT,
        fraction_nonthermal_electrons_: float = defaults.DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
        fraction_nonthermal_protons_: float = defaults.DEFAULT_FRACTION_NONTHERMAL_PROTONS,
        factor_break_electrons_: float = defaults.DEFAULT_FACTOR_BREAK_ELECTRONS,
        factor_break_protons_: float = defaults.DEFAULT_FACTOR_BREAK_PROTONS,
        factor_max_energy_electrons_: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
        factor_max_energy_protons_: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
        index_injected_electrons_: float = defaults.DEFAULT_INDEX_INJECTED_ELECTRONS,
        index_injected_protons_: float = defaults.DEFAULT_INDEX_INJECTED_PROTONS,
        calc_pair_content_from_plasma_beta_: bool = defaults.DEFAULT_CALC_PAIR_CONTENT_FROM_PLASMA_BETA,
        plasma_beta_jet_base_: float = defaults.DEFAULT_PLASMA_BETA_JET_BASE,
        dlgz_: float = defaults.DEFAULT_DLGZ,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL,
        cutoff_type_: int = defaults.DEFAULT_CUTOFF_TYPE
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass



    def calc_velocity_profile_magnetized_jet(self) -> None:
        pass
    def calc_nozzle_energetics_equipartition(self) -> None:
        pass

    def calc_grid_next_zone(self, i: int, cut: int, zcut: float) -> None:
        pass
    def calc_zone_properties(self, i: int) -> None:
        pass

    def compute_jet_dynamics(self) -> None:
        pass


####################    </generated_from:BLJet.hpp>    ####################


####################    <generated_from:IsoJet.hpp>    ####################




class IsoJet:
    # Member variables
    mass_bh: float = defaults.DEFAULT_MASS_BH
    jet_power_eddington: float = defaults.DEFAULT_JET_POWER_EDDINGTON
    z_jet_launching: float = defaults.DEFAULT_Z_JET_LAUNCHING
    r_initial: float = defaults.DEFAULT_R_INITIAL
    z_dissipation: float = defaults.DEFAULT_Z_DISSIPATION
    z_max_calculation: float = defaults.DEFAULT_Z_MAX_CALCULATION
    electron_temperature_jet_base: float = defaults.DEFAULT_ELECTRON_TEMPERATURE_JET_BASE
    plasma_beta_jet_base: float = defaults.DEFAULT_PLASMA_BETA_JET_BASE
    fraction_nonthermal_electrons: float = defaults.DEFAULT_FRACTION_NONTHERMAL_ELECTRONS
    fraction_nonthermal_protons: float = defaults.DEFAULT_FRACTION_NONTHERMAL_PROTONS
    factor_break_electrons: float = defaults.DEFAULT_FACTOR_BREAK_ELECTRONS
    factor_break_protons: float = defaults.DEFAULT_FACTOR_BREAK_PROTONS
    factor_max_energy_electrons: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS
    factor_max_energy_protons: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_PROTONS
    index_injected_electrons: float = defaults.DEFAULT_INDEX_INJECTED_ELECTRONS
    index_injected_protons: float = defaults.DEFAULT_INDEX_INJECTED_PROTONS
    cutoff_type: int = defaults.DEFAULT_CUTOFF_TYPE
    dlgz: float = defaults.DEFAULT_DLGZ

    # internal variables
    eddington_luminosity: float
    # internal variables
    r_g: float
    # internal variables
    zmin: float

    # hardcoded values
    # size_t n_bins_speed = 54;
    # double jet_opening_constant = 0.15;

    def __init__(
        self,
        mass_bh_: float = defaults.DEFAULT_MASS_BH,
        jet_power_eddington_: float = defaults.DEFAULT_JET_POWER_EDDINGTON,
        z_jet_launching_: float = defaults.DEFAULT_Z_JET_LAUNCHING,
        r_initial_: float = defaults.DEFAULT_R_INITIAL,
        z_dissipation_: float = defaults.DEFAULT_Z_DISSIPATION,
        z_max_calculation_: float = defaults.DEFAULT_Z_MAX_CALCULATION,
        electron_temperature_jet_base_: float = defaults.DEFAULT_ELECTRON_TEMPERATURE_JET_BASE,
        plasma_beta_jet_base_: float = defaults.DEFAULT_PLASMA_BETA_JET_BASE,
        fraction_nonthermal_electrons_: float = defaults.DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
        fraction_nonthermal_protons_: float = defaults.DEFAULT_FRACTION_NONTHERMAL_PROTONS,
        factor_break_electrons_: float = defaults.DEFAULT_FACTOR_BREAK_ELECTRONS,
        factor_break_protons_: float = defaults.DEFAULT_FACTOR_BREAK_PROTONS,
        factor_max_energy_electrons_: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
        factor_max_energy_protons_: float = defaults.DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
        index_injected_electrons_: float = defaults.DEFAULT_INDEX_INJECTED_ELECTRONS,
        index_injected_protons_: float = defaults.DEFAULT_INDEX_INJECTED_PROTONS,
        dlgz_: float = defaults.DEFAULT_DLGZ,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL,
        cutoff_type_: int = defaults.DEFAULT_CUTOFF_TYPE
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass



    def calc_velocity_profile_iso_jet(self) -> None:
        pass
    def calc_nozzle_energetics_equipartition(self) -> None:
        pass

    def calc_zone_properties(self, i: int) -> None:
        pass

    def compute_jet_dynamics(self) -> None:
        pass


####################    </generated_from:IsoJet.hpp>    ####################


####################    <generated_from:BHJet.hpp>    ####################


# struct TargetFieldBlackBody
# {
#     double luminosity, temperature, energy_density;
#     std::string name;
#     TargetFieldBlackBody(double l, double t, double u, std::string n)
#         : luminosity(l), temperature(t), energy_density(u), name(n) {}
# };

class BHJet:

    target_list: List[TargetPhotonField]
    def add_target_photon_field(self, target: TargetPhotonField) -> None:
        pass
    def remove_target_photon_field(self, target: TargetPhotonField) -> None:
        pass
    def clear_target_photon_fields(self) -> None:
        pass

    def init_jet_dynamics(self, jet_dynamics_: JetDynamics) -> None:
        pass
    def get_jet_dynamics(self) -> JetDynamics:
        pass
    def compute_full_jet(self, photon_energy_grid: List[float]) -> None:
        pass


    # ----------------------------
    # Member variables
    # ----------------------------
    theta_obs: float
    distance: float
    redshift: float
    compton_threshold: float
    include_counterjet: bool
    profile_time: bool
    verbosity_level: int

    def __init__(
        self,
        theta_obs_: float = defaults.DEFAULT_THETA_OBS,
        distance_: float = defaults.DEFAULT_DISTANCE,
        redshift_: float = defaults.DEFAULT_REDSHIFT,
        include_counterjet_: bool = defaults.DEFAULT_INCLUDE_COUNTERJET,
        compton_threshold_: float = defaults.DEFAULT_COMPTON_THRESHOLD,
        profile_time_: bool = defaults.DEFAULT_PROFILE_TIME,
        verbosity_level_: int = defaults.DEFAULT_VERBOSITY_LEVEL
        ) -> None:
        """ ----------------------------
         Constructor with defaults
         ----------------------------
        """
        pass


    radiation_zones: List[RadiationZone]


    def get_observed_target_photon_energy(self, name: str) -> List[float]:
        pass
    def get_observed_target_photon_flux(self, name: str) -> List[float]:
        pass

    def get_observed_photon_energy_grid(self) -> List[float]:
        pass
    def get_observed_photon_flux_total(self) -> List[float]:
        pass
    def get_observed_photon_integrated_flux_total(
        self,
        z_min: float,
        z_max: float
        ) -> List[float]:
        pass
    def get_observed_photon_integrated_flux_electron_cyclosyn(
        self,
        z_min: float,
        z_max: float
        ) -> List[float]:
        pass
    def get_observed_photon_integrated_flux_electron_compton(
        self,
        z_min: float,
        z_max: float
        ) -> List[float]:
        pass

    computation_times: List[float]
    def get_computation_times(self) -> List[float]:
        pass

####################    </generated_from:BHJet.hpp>    ####################

# </litgen_stub>