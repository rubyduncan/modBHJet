#include <memory>
#include <vector>
#include <string>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/pair.h>

#include "../include/JetDynamics.hpp"
#include "../include/TargetPhotonField.hpp"
#include "../include/TargetBlackBody.hpp"
#include "../include/BLJet.hpp"
#include "../include/IsoJet.hpp"
#include "../include/BHJet.hpp"

namespace nb = nanobind;

// <litgen_glue_code>

// </litgen_glue_code>

void py_init_module_modbhjet(nb::module_ &m)
{
    // <litgen_pydef>
    ////////////////////    <generated_from:JetDynamics.hpp>    ////////////////////
    auto pyClassJetDynamics =
        nb::class_<bhjet::JetDynamics>
            (m, "JetDynamics", "")
        .def_rw("n_zones", &bhjet::JetDynamics::n_zones, "")
        .def_rw("verbosity_level", &bhjet::JetDynamics::verbosity_level, "")
        .def(nb::init<size_t, size_t>(),
            nb::arg("n_zones") = bhjet::defaults::DEFAULT_N_ZONES, nb::arg("verbosity_level") = bhjet::defaults::DEFAULT_VERBOSITY_LEVEL,
            " ----------------------------\n Constructor with defaults\n ----------------------------")
        .def("compute_jet_dynamics",
            &bhjet::JetDynamics::compute_jet_dynamics)
        .def("info",
            &bhjet::JetDynamics::info)
        .def("get_z_min_grid",
            &bhjet::JetDynamics::get_z_min_grid)
        .def("get_z_height_grid",
            &bhjet::JetDynamics::get_z_height_grid)
        .def("get_z_center_grid",
            &bhjet::JetDynamics::get_z_center_grid)
        .def("get_radius_grid",
            &bhjet::JetDynamics::get_radius_grid)
        .def("get_gamma_grid",
            &bhjet::JetDynamics::get_gamma_grid)
        .def("get_beta_grid",
            &bhjet::JetDynamics::get_beta_grid)
        .def("get_beta_gamma_grid",
            &bhjet::JetDynamics::get_beta_gamma_grid)
        .def("get_magnetic_field_grid",
            &bhjet::JetDynamics::get_magnetic_field_grid)
        .def("get_electron_density_grid",
            &bhjet::JetDynamics::get_electron_density_grid)
        .def("get_proton_density_grid",
            &bhjet::JetDynamics::get_proton_density_grid)
        .def("get_electron_temperature_grid",
            &bhjet::JetDynamics::get_electron_temperature_grid)
        .def("get_proton_temperature_grid",
            &bhjet::JetDynamics::get_proton_temperature_grid)
        .def("get_fraction_nonthermal_electrons_grid",
            &bhjet::JetDynamics::get_fraction_nonthermal_electrons_grid)
        .def("get_fraction_nonthermal_protons_grid",
            &bhjet::JetDynamics::get_fraction_nonthermal_protons_grid)
        .def("get_factor_break_electrons_grid",
            &bhjet::JetDynamics::get_factor_break_electrons_grid)
        .def("get_factor_break_protons_grid",
            &bhjet::JetDynamics::get_factor_break_protons_grid)
        .def("get_factor_max_energy_electrons_grid",
            &bhjet::JetDynamics::get_factor_max_energy_electrons_grid)
        .def("get_factor_max_energy_protons_grid",
            &bhjet::JetDynamics::get_factor_max_energy_protons_grid)
        .def("get_index_injected_electrons_grid",
            &bhjet::JetDynamics::get_index_injected_electrons_grid)
        .def("get_index_injected_protons_grid",
            &bhjet::JetDynamics::get_index_injected_protons_grid)
        .def("get_cutoff_type_grid",
            &bhjet::JetDynamics::get_cutoff_type_grid)
        ;
    ////////////////////    </generated_from:JetDynamics.hpp>    ////////////////////


    ////////////////////    <generated_from:TargetPhotonField.hpp>    ////////////////////
    auto pyClassTargetPhotonField =
        nb::class_<bhjet::TargetPhotonField>
            (m, "TargetPhotonField", "")
        .def_rw("name", &bhjet::TargetPhotonField::name, "")
        .def_rw("target_type", &bhjet::TargetPhotonField::target_type, "")
        .def_rw("distance", &bhjet::TargetPhotonField::distance, "")
        .def_rw("redshift", &bhjet::TargetPhotonField::redshift, "")
        .def_rw("add_to_total_flux", &bhjet::TargetPhotonField::add_to_total_flux, "")
        .def_rw("verbosity_level", &bhjet::TargetPhotonField::verbosity_level, "")
        .def("get_target_energy_grid_and_density",
            &bhjet::TargetPhotonField::get_target_energy_grid_and_density, nb::arg("z"), nb::arg("bulk_momentum"), nb::arg("theta_obs"), nb::arg("min_energy"), nb::arg("max_energy"))
        .def("update_observed_flux",
            &bhjet::TargetPhotonField::update_observed_flux)
        .def("get_observed_energy",
            &bhjet::TargetPhotonField::get_observed_energy)
        .def("get_observed_energy_flux",
            &bhjet::TargetPhotonField::get_observed_energy_flux)
        .def("get_observed_number_flux",
            &bhjet::TargetPhotonField::get_observed_number_flux)
        ;
    ////////////////////    </generated_from:TargetPhotonField.hpp>    ////////////////////


    ////////////////////    <generated_from:TargetBlackBody.hpp>    ////////////////////
    auto pyClassTargetBlackBody =
        nb::class_<bhjet::TargetBlackBody, bhjet::TargetPhotonField>
            (m, "TargetBlackBody", "")
        .def_rw("luminosity", &bhjet::TargetBlackBody::luminosity, "")
        .def_rw("temperature", &bhjet::TargetBlackBody::temperature, "")
        .def_rw("energy_density", &bhjet::TargetBlackBody::energy_density, "")
        .def(nb::init<std::string, double, double, double, double, double, bool, size_t>(),
            nb::arg("name"), nb::arg("distance") = bhjet::defaults::DISTANCE, nb::arg("redshift") = bhjet::defaults::REDSHIFT, nb::arg("luminosity") = bhjet::defaults::DEFAULT_LUMINOSTIY, nb::arg("temperature") = bhjet::defaults::DEFAULT_TEMPERATURE, nb::arg("energy_density") = bhjet::defaults::DEFAULT_ENERGY_DENSITY, nb::arg("add_to_total_flux") = bhjet::defaults::DEFAULT_ADD_TO_TOTAL_FLUX, nb::arg("verbosity_level") = bhjet::defaults::DEFAULT_VERBOSITY_LEVEL,
            " ----------------------------\n Constructor with defaults\n ----------------------------")
        .def("get_target_energy_grid_and_density",
            &bhjet::TargetBlackBody::get_target_energy_grid_and_density, nb::arg("z"), nb::arg("bulk_momentum"), nb::arg("theta_obs"), nb::arg("min_energy"), nb::arg("max_energy"))
        .def("update_observed_flux",
            &bhjet::TargetBlackBody::update_observed_flux)
        ;
    ////////////////////    </generated_from:TargetBlackBody.hpp>    ////////////////////


    ////////////////////    <generated_from:RadiationZone.hpp>    ////////////////////
    auto pyClassRadiationZone =
        nb::class_<bhjet::RadiationZone>
            (m, "RadiationZone", "")
        .def_rw("radiation_energy_density", &bhjet::RadiationZone::radiation_energy_density, "total integral (for Thomson cooling)")
        .def("add_target_photon_field",
            &bhjet::RadiationZone::add_target_photon_field, nb::arg("energies"), nb::arg("energy_densities"), nb::arg("name"))
        .def("remove_target_photon_field",
            &bhjet::RadiationZone::remove_target_photon_field, nb::arg("name"))
        .def("compute_particles",
            &bhjet::RadiationZone::compute_particles, nb::arg("reset_photon_targets") = true)
        .def("compute_radiation",
            [](bhjet::RadiationZone & self) { return self.compute_radiation(); })
        .def("compute_radiation",
            nb::overload_cast<std::vector<double>>(&bhjet::RadiationZone::compute_radiation), nb::arg("obs_energy_grid"))
        .def_rw("magnetic_field", &bhjet::RadiationZone::magnetic_field, "")
        .def_rw("radius", &bhjet::RadiationZone::radius, "")
        .def_rw("height", &bhjet::RadiationZone::height, "")
        .def_rw("bulk_momentum", &bhjet::RadiationZone::bulk_momentum, "")
        .def_rw("theta_obs", &bhjet::RadiationZone::theta_obs, "")
        .def_rw("distance", &bhjet::RadiationZone::distance, "")
        .def_rw("redshift", &bhjet::RadiationZone::redshift, "")
        .def_rw("electron_number_density", &bhjet::RadiationZone::electron_number_density, "")
        .def_rw("proton_number_density", &bhjet::RadiationZone::proton_number_density, "")
        .def_rw("electron_temperature", &bhjet::RadiationZone::electron_temperature, "")
        .def_rw("proton_temperature", &bhjet::RadiationZone::proton_temperature, "")
        .def_rw("fraction_nonthermal_electrons", &bhjet::RadiationZone::fraction_nonthermal_electrons, "")
        .def_rw("fraction_nonthermal_protons", &bhjet::RadiationZone::fraction_nonthermal_protons, "")
        .def_rw("factor_break_electrons", &bhjet::RadiationZone::factor_break_electrons, "")
        .def_rw("factor_break_protons", &bhjet::RadiationZone::factor_break_protons, "")
        .def_rw("factor_max_energy_electrons", &bhjet::RadiationZone::factor_max_energy_electrons, "")
        .def_rw("factor_max_energy_protons", &bhjet::RadiationZone::factor_max_energy_protons, "")
        .def_rw("index_injected_electrons", &bhjet::RadiationZone::index_injected_electrons, "")
        .def_rw("index_injected_protons", &bhjet::RadiationZone::index_injected_protons, "")
        .def_rw("cutoff_type", &bhjet::RadiationZone::cutoff_type, "")
        .def_rw("geometry", &bhjet::RadiationZone::geometry, "")
        .def_rw("include_counterjet", &bhjet::RadiationZone::include_counterjet, "")
        .def_rw("force_compton_calculation", &bhjet::RadiationZone::force_compton_calculation, "")
        .def_rw("compton_switch", &bhjet::RadiationZone::compton_switch, "")
        .def_rw("compton_threshold", &bhjet::RadiationZone::compton_threshold, "")
        .def_rw("profile_time", &bhjet::RadiationZone::profile_time, "")
        .def_rw("verbosity_level", &bhjet::RadiationZone::verbosity_level, "")
        .def_rw("n_bins_e", &bhjet::RadiationZone::n_bins_e, "momentum grid")
        .def_rw("n_bins_p", &bhjet::RadiationZone::n_bins_p, "momentum grid")
        .def_rw("syn_res", &bhjet::RadiationZone::syn_res, "")
        .def_rw("com_res", &bhjet::RadiationZone::com_res, "")
        .def_rw("doppler_factor_bulk", &bhjet::RadiationZone::doppler_factor_bulk, "")
        .def_rw("beta_bulk", &bhjet::RadiationZone::beta_bulk, "")
        .def_rw("gamma_bulk", &bhjet::RadiationZone::gamma_bulk, "")
        .def_rw("computation_times", &bhjet::RadiationZone::computation_times, "")
        .def(nb::init<double, double, double, std::string, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, bool, bool, bool, double, bool, size_t, int>(),
            nb::arg("magnetic_field") = bhjet::defaults::DEFAULT_MAGNETIC_FIELD, nb::arg("radius") = bhjet::defaults::DEFAULT_RADIUS, nb::arg("height") = bhjet::defaults::DEFAULT_HEIGHT, nb::arg("geometry") = bhjet::defaults::DEFAULT_GEOMETRY, nb::arg("bulk_momentum") = bhjet::defaults::DEFAULT_BULK_MOMENTUM, nb::arg("theta_obs") = bhjet::defaults::DEFAULT_THETA_OBS, nb::arg("distance") = bhjet::defaults::DEFAULT_DISTANCE, nb::arg("redshift") = bhjet::defaults::DEFAULT_REDSHIFT, nb::arg("electron_number_density") = bhjet::defaults::DEFAULT_ELECTRON_NUMBER_DENSITY, nb::arg("proton_number_density") = bhjet::defaults::DEFAULT_PROTON_NUMBER_DENSITY, nb::arg("electron_temperature") = bhjet::defaults::DEFAULT_ELECTRON_TEMPERATURE, nb::arg("proton_temperature") = bhjet::defaults::DEFAULT_PROTON_TEMPERATURE, nb::arg("fraction_nonthermal_electrons") = bhjet::defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS, nb::arg("fraction_nonthermal_protons") = bhjet::defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS, nb::arg("factor_break_electrons") = bhjet::defaults::DEFAULT_FACTOR_BREAK_ELECTRONS, nb::arg("factor_break_protons") = bhjet::defaults::DEFAULT_FACTOR_BREAK_PROTONS, nb::arg("factor_max_energy_electrons") = bhjet::defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS, nb::arg("factor_max_energy_protons") = bhjet::defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS, nb::arg("index_injected_electrons") = bhjet::defaults::DEFAULT_INDEX_INJECTED_ELECTRONS, nb::arg("index_injected_protons") = bhjet::defaults::DEFAULT_INDEX_INJECTED_PROTONS, nb::arg("include_counterjet") = bhjet::defaults::DEFAULT_INCLUDE_COUNTERJET, nb::arg("force_compton_calculation") = bhjet::defaults::DEFAULT_FORCE_COMPTON_CALCULATION, nb::arg("compton_switch") = bhjet::defaults::DEFAULT_COMPTON_SWITCH, nb::arg("compton_threshold") = bhjet::defaults::DEFAULT_COMPTON_THRESHOLD, nb::arg("profile_time") = bhjet::defaults::DEFAULT_PROFILE_TIME, nb::arg("verbosity_level") = bhjet::defaults::DEFAULT_VERBOSITY_LEVEL, nb::arg("cutoff_type") = bhjet::defaults::DEFAULT_CUTOFF_TYPE,
            " ----------------------------\n Constructor with defaults\n ----------------------------")
        .def("get_timescale_electron_cyclosyn",
            &bhjet::RadiationZone::get_timescale_electron_cyclosyn, nb::arg("momentum"))
        .def("get_timescale_electron_adiabatic",
            &bhjet::RadiationZone::get_timescale_electron_adiabatic, nb::arg("momentum"))
        .def("get_timescale_electron_compton_thomson",
            &bhjet::RadiationZone::get_timescale_electron_compton_thomson, nb::arg("momentum"))
        .def("get_timescale_electron_compton",
            &bhjet::RadiationZone::get_timescale_electron_compton, nb::arg("momentum"))
        .def("get_timescale_electron_acceleration",
            &bhjet::RadiationZone::get_timescale_electron_acceleration, nb::arg("momentum"))
        .def("get_timescale_photon_cyclosyn_selfabsorption",
            &bhjet::RadiationZone::get_timescale_photon_cyclosyn_selfabsorption, nb::arg("momentum"))
        .def("get_timescale_photon_escape",
            &bhjet::RadiationZone::get_timescale_photon_escape,
            nb::arg("momentum"),
            " read out the array, left for debugging\n std::vector<double> get_timescale_photon_cyclosyn_selfabsorption_array(){return cyclosyn_selfabsorption_rate;};\n std::vector<double> get_timescale_photon_cyclosyn_selfabsorption_energy(){return cyclosyn_energy;};")
        .def("get_electron_max_momentum",
            &bhjet::RadiationZone::get_electron_max_momentum)
        .def("get_electron_break_momentum",
            &bhjet::RadiationZone::get_electron_break_momentum)
        .def("get_electron_momentum_grid",
            &bhjet::RadiationZone::get_electron_momentum_grid)
        .def("get_electron_gamma_grid",
            &bhjet::RadiationZone::get_electron_gamma_grid)
        .def("get_electron_momentum_number_density",
            &bhjet::RadiationZone::get_electron_momentum_number_density)
        .def("get_electron_gamma_number_density",
            &bhjet::RadiationZone::get_electron_gamma_number_density)
        .def("get_total_photon_target_energy",
            &bhjet::RadiationZone::get_total_photon_target_energy)
        .def("get_total_photon_target_energy_density",
            &bhjet::RadiationZone::get_total_photon_target_energy_density)
        .def("get_photon_target_energy",
            &bhjet::RadiationZone::get_photon_target_energy, nb::arg("name"))
        .def("get_photon_target_energy_density",
            &bhjet::RadiationZone::get_photon_target_energy_density, nb::arg("name"))
        .def("get_observed_photon_energy_grid_electron_cyclosyn",
            &bhjet::RadiationZone::get_observed_photon_energy_grid_electron_cyclosyn)
        .def("get_observed_photon_luminosity_electron_cyclosyn",
            &bhjet::RadiationZone::get_observed_photon_luminosity_electron_cyclosyn)
        .def("get_observed_photon_flux_electron_cyclosyn",
            &bhjet::RadiationZone::get_observed_photon_flux_electron_cyclosyn)
        .def("get_observed_photon_energy_grid_electron_compton",
            &bhjet::RadiationZone::get_observed_photon_energy_grid_electron_compton)
        .def("get_observed_photon_luminosity_electron_compton",
            &bhjet::RadiationZone::get_observed_photon_luminosity_electron_compton)
        .def("get_observed_photon_flux_electron_compton",
            &bhjet::RadiationZone::get_observed_photon_flux_electron_compton)
        .def("get_observed_photon_energy_grid_total",
            &bhjet::RadiationZone::get_observed_photon_energy_grid_total)
        .def("get_observed_photon_luminosity_total",
            &bhjet::RadiationZone::get_observed_photon_luminosity_total)
        .def("get_observed_photon_flux_total",
            &bhjet::RadiationZone::get_observed_photon_flux_total)
        .def("get_computation_times",
            &bhjet::RadiationZone::get_computation_times)
        .def("sum_jet_and_counterjet",
            &bhjet::RadiationZone::sum_jet_and_counterjet, nb::arg("size"), nb::arg("input_en"), nb::arg("input_lum"), nb::arg("en"), nb::arg("lum"))
        .def("sum_jet_only",
            &bhjet::RadiationZone::sum_jet_only, nb::arg("size"), nb::arg("input_en"), nb::arg("input_lum"), nb::arg("en"), nb::arg("lum"))
        .def("compton_calculation_necessary",
            &bhjet::RadiationZone::compton_calculation_necessary)
        .def("fic",
            &bhjet::RadiationZone::Fic,
            nb::arg("gamma"), nb::arg("alpha"),
            "*\n         * @brief function to compute ker_ic[][].\n         *\n         * cf. Jones, \"Inverse Compton Scattering of Cosmic-Ray Electrons\", P.R. 1965, Eqn. 13-15\n         *\n         * @param gamma Energy of incoming particle\n         * @param alpha Energy of incoming photon / me c^2\n")
        .def("fic_1",
            &bhjet::RadiationZone::Fic_1,
            nb::arg("z"),
            "< utility functions to construct Fic")
        .def("fic_2",
            &bhjet::RadiationZone::Fic_2,
            nb::arg("z"),
            "< utility functions to construct Fic")
        .def("polylog",
            &bhjet::RadiationZone::polylog,
            nb::arg("z"),
            "utility functions to construct Fic, valid for -inf<z<1.0")
        .def("polylog_smallz",
            &bhjet::RadiationZone::polylog_smallz,
            nb::arg("z"),
            "< utility functions to calc polylog")
        .def("polylog_rk4",
            &bhjet::RadiationZone::polylog_rk4,
            nb::arg("z"),
            "< utility functions to construct polylog")
        .def("polylog_largez",
            &bhjet::RadiationZone::polylog_largez,
            nb::arg("z"),
            "< utility functions to construct polylog")
        ;
    ////////////////////    </generated_from:RadiationZone.hpp>    ////////////////////


    ////////////////////    <generated_from:BLJet.hpp>    ////////////////////
    auto pyClassBLJet =
        nb::class_<bhjet::BLJet, bhjet::JetDynamics>
            (m, "BLJet", "")
        .def_rw("mass_bh", &bhjet::BLJet::mass_bh, "")
        .def_rw("jet_power_eddington", &bhjet::BLJet::jet_power_eddington, "")
        .def_rw("z_jet_launching", &bhjet::BLJet::z_jet_launching, "")
        .def_rw("r_initial", &bhjet::BLJet::r_initial, "")
        .def_rw("z_end_of_acceleration", &bhjet::BLJet::z_end_of_acceleration, "")
        .def_rw("z_dissipation", &bhjet::BLJet::z_dissipation, "")
        .def_rw("z_max_calculation", &bhjet::BLJet::z_max_calculation, "")
        .def_rw("sigma_final", &bhjet::BLJet::sigma_final, "")
        .def_rw("gamma_final", &bhjet::BLJet::gamma_final, "")
        .def_rw("electron_temperature_jet_base", &bhjet::BLJet::electron_temperature_jet_base, "")
        .def_rw("gamma_acceleration_exponent", &bhjet::BLJet::gamma_acceleration_exponent, "")
        .def_rw("gamma_deceleration_exponent", &bhjet::BLJet::gamma_deceleration_exponent, "")
        .def_rw("opening_angle_constant", &bhjet::BLJet::opening_angle_constant, "")
        .def_rw("fraction_nonthermal_electrons", &bhjet::BLJet::fraction_nonthermal_electrons, "")
        .def_rw("fraction_nonthermal_protons", &bhjet::BLJet::fraction_nonthermal_protons, "")
        .def_rw("factor_break_electrons", &bhjet::BLJet::factor_break_electrons, "")
        .def_rw("factor_break_protons", &bhjet::BLJet::factor_break_protons, "")
        .def_rw("factor_max_energy_electrons", &bhjet::BLJet::factor_max_energy_electrons, "")
        .def_rw("factor_max_energy_protons", &bhjet::BLJet::factor_max_energy_protons, "")
        .def_rw("index_injected_electrons", &bhjet::BLJet::index_injected_electrons, "")
        .def_rw("index_injected_protons", &bhjet::BLJet::index_injected_protons, "")
        .def_rw("cutoff_type", &bhjet::BLJet::cutoff_type, "")
        .def_rw("calc_pair_content_from_plasma_beta", &bhjet::BLJet::calc_pair_content_from_plasma_beta, "")
        .def_rw("plasma_beta_jet_base", &bhjet::BLJet::plasma_beta_jet_base, "")
        .def_rw("dlgz", &bhjet::BLJet::dlgz, "")
        .def_rw("eddington_luminosity", &bhjet::BLJet::eddington_luminosity, "internal variables")
        .def_rw("r_g", &bhjet::BLJet::r_g, "internal variables")
        .def_rw("zmin", &bhjet::BLJet::zmin, "internal variables")
        .def(nb::init<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, bool, double, double, size_t, int>(),
            nb::arg("mass_bh") = bhjet::defaults::DEFAULT_MASS_BH, nb::arg("jet_power_eddington") = bhjet::defaults::DEFAULT_JET_POWER_EDDINGTON, nb::arg("z_jet_launching") = bhjet::defaults::DEFAULT_Z_JET_LAUNCHING, nb::arg("r_initial") = bhjet::defaults::DEFAULT_R_INITIAL, nb::arg("z_end_of_acceleration") = bhjet::defaults::DEFAULT_Z_END_OF_ACCELERATION, nb::arg("z_dissipation") = bhjet::defaults::DEFAULT_Z_DISSIPATION, nb::arg("z_max_calculation") = bhjet::defaults::DEFAULT_Z_MAX_CALCULATION, nb::arg("sigma_final") = bhjet::defaults::DEFAULT_SIGMA_FINAL, nb::arg("gamma_final") = bhjet::defaults::DEFAULT_GAMMA_FINAL, nb::arg("electron_temperature_jet_base") = bhjet::defaults::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE, nb::arg("gamma_acceleration_exponent") = bhjet::defaults::DEFAULT_GAMMA_ACCELERATION_EXPONENT, nb::arg("gamma_deceleration_exponent") = bhjet::defaults::DEFAULT_GAMMA_DECELERATION_EXPONENT, nb::arg("opening_angle_constant") = bhjet::defaults::DEFAULT_OPENING_ANGLE_CONSTANT, nb::arg("fraction_nonthermal_electrons") = bhjet::defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS, nb::arg("fraction_nonthermal_protons") = bhjet::defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS, nb::arg("factor_break_electrons") = bhjet::defaults::DEFAULT_FACTOR_BREAK_ELECTRONS, nb::arg("factor_break_protons") = bhjet::defaults::DEFAULT_FACTOR_BREAK_PROTONS, nb::arg("factor_max_energy_electrons") = bhjet::defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS, nb::arg("factor_max_energy_protons") = bhjet::defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS, nb::arg("index_injected_electrons") = bhjet::defaults::DEFAULT_INDEX_INJECTED_ELECTRONS, nb::arg("index_injected_protons") = bhjet::defaults::DEFAULT_INDEX_INJECTED_PROTONS, nb::arg("calc_pair_content_from_plasma_beta") = bhjet::defaults::DEFAULT_CALC_PAIR_CONTENT_FROM_PLASMA_BETA, nb::arg("plasma_beta_jet_base") = bhjet::defaults::DEFAULT_PLASMA_BETA_JET_BASE, nb::arg("dlgz") = bhjet::defaults::DEFAULT_DLGZ, nb::arg("verbosity_level") = bhjet::defaults::DEFAULT_VERBOSITY_LEVEL, nb::arg("cutoff_type") = bhjet::defaults::DEFAULT_CUTOFF_TYPE,
            " ----------------------------\n Constructor with defaults\n ----------------------------")
        .def("calc_velocity_profile_magnetized_jet",
            &bhjet::BLJet::calc_velocity_profile_magnetized_jet)
        .def("calc_nozzle_energetics_equipartition",
            &bhjet::BLJet::calc_nozzle_energetics_equipartition)
        .def("calc_grid_next_zone",
            &bhjet::BLJet::calc_grid_next_zone, nb::arg("i"), nb::arg("cut"), nb::arg("zcut"))
        .def("calc_zone_properties",
            &bhjet::BLJet::calc_zone_properties, nb::arg("i"))
        .def("compute_jet_dynamics",
            &bhjet::BLJet::compute_jet_dynamics)
        ;
    ////////////////////    </generated_from:BLJet.hpp>    ////////////////////


    ////////////////////    <generated_from:IsoJet.hpp>    ////////////////////
    auto pyClassIsoJet =
        nb::class_<bhjet::IsoJet, bhjet::JetDynamics>
            (m, "IsoJet", "")
        .def_rw("mass_bh", &bhjet::IsoJet::mass_bh, "")
        .def_rw("jet_power_eddington", &bhjet::IsoJet::jet_power_eddington, "")
        .def_rw("z_jet_launching", &bhjet::IsoJet::z_jet_launching, "")
        .def_rw("r_initial", &bhjet::IsoJet::r_initial, "")
        .def_rw("z_dissipation", &bhjet::IsoJet::z_dissipation, "")
        .def_rw("z_max_calculation", &bhjet::IsoJet::z_max_calculation, "")
        .def_rw("electron_temperature_jet_base", &bhjet::IsoJet::electron_temperature_jet_base, "")
        .def_rw("plasma_beta_jet_base", &bhjet::IsoJet::plasma_beta_jet_base, "")
        .def_rw("fraction_nonthermal_electrons", &bhjet::IsoJet::fraction_nonthermal_electrons, "")
        .def_rw("fraction_nonthermal_protons", &bhjet::IsoJet::fraction_nonthermal_protons, "")
        .def_rw("factor_break_electrons", &bhjet::IsoJet::factor_break_electrons, "")
        .def_rw("factor_break_protons", &bhjet::IsoJet::factor_break_protons, "")
        .def_rw("factor_max_energy_electrons", &bhjet::IsoJet::factor_max_energy_electrons, "")
        .def_rw("factor_max_energy_protons", &bhjet::IsoJet::factor_max_energy_protons, "")
        .def_rw("index_injected_electrons", &bhjet::IsoJet::index_injected_electrons, "")
        .def_rw("index_injected_protons", &bhjet::IsoJet::index_injected_protons, "")
        .def_rw("cutoff_type", &bhjet::IsoJet::cutoff_type, "")
        .def_rw("dlgz", &bhjet::IsoJet::dlgz, "")
        .def_rw("eddington_luminosity", &bhjet::IsoJet::eddington_luminosity, "internal variables")
        .def_rw("r_g", &bhjet::IsoJet::r_g, "internal variables")
        .def_rw("zmin", &bhjet::IsoJet::zmin, "internal variables")
        .def(nb::init<double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double, size_t, int>(),
            nb::arg("mass_bh") = bhjet::defaults::DEFAULT_MASS_BH, nb::arg("jet_power_eddington") = bhjet::defaults::DEFAULT_JET_POWER_EDDINGTON, nb::arg("z_jet_launching") = bhjet::defaults::DEFAULT_Z_JET_LAUNCHING, nb::arg("r_initial") = bhjet::defaults::DEFAULT_R_INITIAL, nb::arg("z_dissipation") = bhjet::defaults::DEFAULT_Z_DISSIPATION, nb::arg("z_max_calculation") = bhjet::defaults::DEFAULT_Z_MAX_CALCULATION, nb::arg("electron_temperature_jet_base") = bhjet::defaults::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE, nb::arg("plasma_beta_jet_base") = bhjet::defaults::DEFAULT_PLASMA_BETA_JET_BASE, nb::arg("fraction_nonthermal_electrons") = bhjet::defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS, nb::arg("fraction_nonthermal_protons") = bhjet::defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS, nb::arg("factor_break_electrons") = bhjet::defaults::DEFAULT_FACTOR_BREAK_ELECTRONS, nb::arg("factor_break_protons") = bhjet::defaults::DEFAULT_FACTOR_BREAK_PROTONS, nb::arg("factor_max_energy_electrons") = bhjet::defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS, nb::arg("factor_max_energy_protons") = bhjet::defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS, nb::arg("index_injected_electrons") = bhjet::defaults::DEFAULT_INDEX_INJECTED_ELECTRONS, nb::arg("index_injected_protons") = bhjet::defaults::DEFAULT_INDEX_INJECTED_PROTONS, nb::arg("dlgz") = bhjet::defaults::DEFAULT_DLGZ, nb::arg("verbosity_level") = bhjet::defaults::DEFAULT_VERBOSITY_LEVEL, nb::arg("cutoff_type") = bhjet::defaults::DEFAULT_CUTOFF_TYPE,
            " ----------------------------\n Constructor with defaults\n ----------------------------")
        .def("calc_velocity_profile_iso_jet",
            &bhjet::IsoJet::calc_velocity_profile_iso_jet)
        .def("calc_nozzle_energetics_equipartition",
            &bhjet::IsoJet::calc_nozzle_energetics_equipartition)
        .def("calc_zone_properties",
            &bhjet::IsoJet::calc_zone_properties, nb::arg("i"))
        .def("compute_jet_dynamics",
            &bhjet::IsoJet::compute_jet_dynamics)
        ;
    ////////////////////    </generated_from:IsoJet.hpp>    ////////////////////


    ////////////////////    <generated_from:BHJet.hpp>    ////////////////////
    auto pyClassBHJet =
        nb::class_<bhjet::BHJet>
            (m, "BHJet", "")
        .def_rw("target_list", &bhjet::BHJet::target_list, "")
        .def("add_target_photon_field",
            &bhjet::BHJet::add_target_photon_field, nb::arg("target"))
        .def("remove_target_photon_field",
            &bhjet::BHJet::remove_target_photon_field, nb::arg("target"))
        .def("clear_target_photon_fields",
            &bhjet::BHJet::clear_target_photon_fields)
        .def("init_jet_dynamics",
            &bhjet::BHJet::init_jet_dynamics, nb::arg("jet_dynamics"))
        .def("get_jet_dynamics",
            &bhjet::BHJet::get_jet_dynamics)
        .def("compute_full_jet",
            &bhjet::BHJet::compute_full_jet, nb::arg("photon_energy_grid"))
        .def_rw("theta_obs", &bhjet::BHJet::theta_obs, "")
        .def_rw("distance", &bhjet::BHJet::distance, "")
        .def_rw("redshift", &bhjet::BHJet::redshift, "")
        .def_rw("compton_threshold", &bhjet::BHJet::compton_threshold, "")
        .def_rw("include_counterjet", &bhjet::BHJet::include_counterjet, "")
        .def_rw("profile_time", &bhjet::BHJet::profile_time, "")
        .def_rw("verbosity_level", &bhjet::BHJet::verbosity_level, "")
        .def(nb::init<double, double, double, bool, double, bool, size_t>(),
            nb::arg("theta_obs") = bhjet::defaults::DEFAULT_THETA_OBS, nb::arg("distance") = bhjet::defaults::DEFAULT_DISTANCE, nb::arg("redshift") = bhjet::defaults::DEFAULT_REDSHIFT, nb::arg("include_counterjet") = bhjet::defaults::DEFAULT_INCLUDE_COUNTERJET, nb::arg("compton_threshold") = bhjet::defaults::DEFAULT_COMPTON_THRESHOLD, nb::arg("profile_time") = bhjet::defaults::DEFAULT_PROFILE_TIME, nb::arg("verbosity_level") = bhjet::defaults::DEFAULT_VERBOSITY_LEVEL,
            " ----------------------------\n Constructor with defaults\n ----------------------------")
        .def_rw("radiation_zones", &bhjet::BHJet::radiation_zones, "")
        .def("get_observed_target_photon_energy",
            &bhjet::BHJet::get_observed_target_photon_energy, nb::arg("name"))
        .def("get_observed_target_photon_flux",
            &bhjet::BHJet::get_observed_target_photon_flux, nb::arg("name"))
        .def("get_observed_photon_energy_grid",
            &bhjet::BHJet::get_observed_photon_energy_grid)
        .def("get_observed_photon_flux_total",
            &bhjet::BHJet::get_observed_photon_flux_total)
        .def("get_observed_photon_integrated_flux_total",
            &bhjet::BHJet::get_observed_photon_integrated_flux_total, nb::arg("z_min"), nb::arg("z_max"))
        .def("get_observed_photon_integrated_flux_electron_cyclosyn",
            &bhjet::BHJet::get_observed_photon_integrated_flux_electron_cyclosyn, nb::arg("z_min"), nb::arg("z_max"))
        .def("get_observed_photon_integrated_flux_electron_compton",
            &bhjet::BHJet::get_observed_photon_integrated_flux_electron_compton, nb::arg("z_min"), nb::arg("z_max"))
        .def_rw("computation_times", &bhjet::BHJet::computation_times, "")
        .def("get_computation_times",
            &bhjet::BHJet::get_computation_times)
        ;
    ////////////////////    </generated_from:BHJet.hpp>    ////////////////////

    // </litgen_pydef>
}