#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "JetDynamics.hpp"
#include "BLJet.hpp"
#include "IsoJet.hpp"
#include "RadiationZone.hpp"
#include "BHJet.hpp"
#include "TargetPhotonField.hpp"
#include "TargetBlackBody.hpp"
#include "utils.hpp"
#include "default_values.hpp"

namespace py = pybind11;
using namespace bhjet;

// ----------------------------
// Macros for parameter lists
// ----------------------------

// This is needed to build a list of the arguments with no comma at the last one
#define SEP_COMMA ,

#define RADIATIONZONE_PARAMS                                                                                                                                                                                                                     \
    X(magnetic_field, double, defaults::DEFAULT_MAGNETIC_FIELD, "Turbulent magnetic field [G]", SEP_COMMA)                                                                                                                                  \
    X(radius, double, defaults::DEFAULT_RADIUS, "Radius of the zone [cm]", SEP_COMMA)                                                                                                                                                       \
    X(height, double, defaults::DEFAULT_HEIGHT, "Height of the zone (only used for cylindrical geometry) [cm]", SEP_COMMA)                                                                                                                  \
    X(geometry, std::string, defaults::DEFAULT_GEOMETRY, "Geometry of zone. Default: 'sphere', alternative 'cylinder'", SEP_COMMA)                                                                                                          \
    X(bulk_momentum, double, defaults::DEFAULT_BULK_MOMENTUM, "Bulk speed as beta*gamma, with beta=speed/c and gamma^2=1/(1-beta^2)", SEP_COMMA)                                                                                            \
    X(theta_obs, double, defaults::DEFAULT_THETA_OBS, "Observation angle [degree]", SEP_COMMA)                                                                                                                                              \
    X(distance, double, defaults::DEFAULT_DISTANCE, "Distance to zone [kpc]", SEP_COMMA)                                                                                                                                                    \
    X(redshift, double, defaults::DEFAULT_REDSHIFT, "Redshift of zone", SEP_COMMA)                                                                                                                                                          \
    X(electron_number_density, double, defaults::DEFAULT_ELECTRON_NUMBER_DENSITY, "Number density of electrons in zone [1/cm³]", SEP_COMMA)                                                                                                 \
    X(proton_number_density, double, defaults::DEFAULT_PROTON_NUMBER_DENSITY, "Number density of protons in zone [1/cm³]", SEP_COMMA)                                                                                                       \
    X(electron_temperature, double, defaults::DEFAULT_ELECTRON_TEMPERATURE, "Temperature of electrons in zone [keV]", SEP_COMMA)                                                                                                            \
    X(proton_temperature, double, defaults::DEFAULT_PROTON_TEMPERATURE, "Temperature of protons in zone [keV]", SEP_COMMA)                                                                                                                  \
    X(fraction_nonthermal_electrons, double, defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS, "Fraction of energy in non-thermal electron tail", SEP_COMMA)                                                                                 \
    X(fraction_nonthermal_protons, double, defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS, "Fraction of energy in non-thermal proton tail", SEP_COMMA)                                                                                       \
    X(factor_break_electrons, double, defaults::DEFAULT_FACTOR_BREAK_ELECTRONS, "Scaling factor for electron adiabtic timescale", SEP_COMMA)                                                                                                \
    X(factor_break_protons, double, defaults::DEFAULT_FACTOR_BREAK_PROTONS, "Scaling factor for proton adiabtic timescale", SEP_COMMA)                                                                                                      \
    X(factor_max_energy_electrons, double, defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS, "Scaling factor for electron acceleration rate, translating to an effective scaling of the maximum energy", SEP_COMMA)                            \
    X(factor_max_energy_protons, double, defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS, "Scaling factor for proton acceleration rate, translating to an effective scaling of the maximum energy", SEP_COMMA)                                  \
    X(index_injected_electrons, double, defaults::DEFAULT_INDEX_INJECTED_ELECTRONS, "Injected electron spectral index (dlogN/dlogE), ie. before cooling", SEP_COMMA)                                                                        \
    X(index_injected_protons, double, defaults::DEFAULT_INDEX_INJECTED_PROTONS, "Injected proton spectral index (dlogN/dlogE), ie. before cooling", SEP_COMMA)                                                                              \
    X(include_counterjet, bool, defaults::DEFAULT_INCLUDE_COUNTERJET, "True: Includes the emission of the counterjet; False: Includes only one jet", SEP_COMMA)                                                                             \
    X(force_compton_calculation, bool, defaults::DEFAULT_FORCE_COMPTON_CALCULATION, "True: Forces the Compton emission to be computed or not based on compton_switch; False: Uses internal criteria and ignores compton_switch", SEP_COMMA) \
    X(compton_switch, bool, defaults::DEFAULT_COMPTON_SWITCH, "True: Forces the Compton emission to be computed; False: Skips Compton emission computation", SEP_COMMA)                                                                     \
    X(compton_threshold, double, defaults::DEFAULT_COMPTON_THRESHOLD, "Internal rough threshold for Compton emission computation: L_com/L_syn > compton_threshold ? Do calculation. Otheriwse skip.", SEP_COMMA)                            \
    X(profile_time, bool, defaults::DEFAULT_PROFILE_TIME, "True: measures computation time of multiple emission processes; False: Does nothing", SEP_COMMA)                                                                                 \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, "Regulates print output of the code. 0: No output; 1: Only important warnings; 2: More output; 3: Debugging output", SEP_COMMA) \
    X(cutoff_type, int, defaults::DEFAULT_CUTOFF_TYPE, "Particle cutoff shape. 0: exp(-x), 1: exp(-x^2), 2: sech(x)^2.", ) // leave the last one empty

#define JETDYNAMICS_PARAMS                                \
    X(n_zones, size_t, defaults::DEFAULT_N_ZONES, SEP_COMMA) \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, )

#define BLJET_PARAMS                                                                                                                                                                                                  \
    X(mass_bh, double, defaults::DEFAULT_MASS_BH, "Black hole mass in units of solar mass", SEP_COMMA)                                                                                                                   \
    X(jet_power_eddington, double, defaults::DEFAULT_JET_POWER_EDDINGTON, "Jet power in units of the Eddington luminosity", SEP_COMMA)                                                                                   \
    X(z_jet_launching, double, defaults::DEFAULT_Z_JET_LAUNCHING, "Distance from black hole where the jet starts [rg]", SEP_COMMA)                                                                                       \
    X(r_initial, double, defaults::DEFAULT_R_INITIAL, "Jet radius when the jet starts [rg]", SEP_COMMA)                                                                                                                  \
    X(z_end_of_acceleration, double, defaults::DEFAULT_Z_END_OF_ACCELERATION, "Distance from black hole where the bulk jet accelerated to gamma_final [rg]", SEP_COMMA)                                                  \
    X(z_dissipation, double, defaults::DEFAULT_Z_DISSIPATION, "Distance from black hole where the jet starts to dissipate energy into non-thermal particles [rg]", SEP_COMMA)                                            \
    X(z_max_calculation, double, defaults::DEFAULT_Z_MAX_CALCULATION, "Distance from black hole where the calulation stops [rg]", SEP_COMMA)                                                                             \
    X(sigma_final, double, defaults::DEFAULT_SIGMA_FINAL, "Magnetisation at z_end_of_acceleration", SEP_COMMA)                                                                                                           \
    X(gamma_final, double, defaults::DEFAULT_GAMMA_FINAL, "Bulk Lorentz factor at z_end_of_acceleration", SEP_COMMA)                                                                                                     \
    X(electron_temperature_jet_base, double, defaults::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE, "Electron temperature value at z_jet_launching (the jet base) [keV]", SEP_COMMA)                                           \
    X(gamma_acceleration_exponent, double, defaults::DEFAULT_GAMMA_ACCELERATION_EXPONENT, "Exponent alpha of jet acceleration profile, gamma_bulk propto z^alpha", SEP_COMMA)                                            \
    X(gamma_deceleration_exponent, double, defaults::DEFAULT_GAMMA_DECELERATION_EXPONENT, "  ", SEP_COMMA)                                                                                                               \
    X(opening_angle_constant, double, defaults::DEFAULT_OPENING_ANGLE_CONSTANT, "Jet opening angle used to convert Lorentz factor to jet radius, r=opening_angle_constant / gamma", SEP_COMMA)                           \
    X(fraction_nonthermal_electrons, double, defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS, "Fraction of energy in non-thermal electron tail", SEP_COMMA)                                                      \
    X(fraction_nonthermal_protons, double, defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS, "Fraction of energy in non-thermal proton tail", SEP_COMMA)                                                            \
    X(factor_break_electrons, double, defaults::DEFAULT_FACTOR_BREAK_ELECTRONS, "Scaling factor for electron adiabtic timescale", SEP_COMMA)                                                                     \
    X(factor_break_protons, double, defaults::DEFAULT_FACTOR_BREAK_PROTONS, "Scaling factor for proton adiabtic timescale", SEP_COMMA)                                                                           \
    X(factor_max_energy_electrons, double, defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS, "Scaling factor for electron acceleration rate, translating to an effective scaling of the maximum energy", SEP_COMMA) \
    X(factor_max_energy_protons, double, defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS, "Scaling factor for proton acceleration rate, translating to an effective scaling of the maximum energy", SEP_COMMA)       \
    X(index_injected_electrons, double, defaults::DEFAULT_INDEX_INJECTED_ELECTRONS, "Injected electron spectral index (dlogN/dlogE), ie. before cooling", SEP_COMMA)                                             \
    X(index_injected_protons, double, defaults::DEFAULT_INDEX_INJECTED_PROTONS, "Injected proton spectral index (dlogN/dlogE), ie. before cooling", SEP_COMMA)                                                   \
    X(calc_pair_content_from_plasma_beta, bool, defaults::DEFAULT_CALC_PAIR_CONTENT_FROM_PLASMA_BETA, "True: use the plasma beta variable to estimate the pair content (n_p/n_e) at the jet base; False: Assume same number density for electrons and protons and ignore plasma_beta_jet_base", SEP_COMMA)     \
    X(plasma_beta_jet_base, double, defaults::DEFAULT_PLASMA_BETA_JET_BASE, "Plasma beta value at z_jet_launching (the jet base)", SEP_COMMA)                                                                            \
    X(dlgz, double, defaults::DEFAULT_DLGZ, "log. grid spacing log10(z[i+1]) - log10(z[i])", SEP_COMMA) \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, "Regulates print output of the code. 0: No output; 1: Only important warnings; 2: More output; 3: Debugging output", SEP_COMMA) \
    X(cutoff_type, int, defaults::DEFAULT_CUTOFF_TYPE, "Particle cutoff shape. 0: exp(-x), 1: exp(-x^2), 2: sech(x)^2.", )

#define ISOJET_PARAMS                                                                                                                                                                                                  \
    X(mass_bh, double, defaults::DEFAULT_MASS_BH, "Black hole mass in units of solar mass", SEP_COMMA)                                                                                                                   \
    X(jet_power_eddington, double, defaults::DEFAULT_JET_POWER_EDDINGTON, "Jet power in units of the Eddington luminosity", SEP_COMMA)                                                                                   \
    X(z_jet_launching, double, defaults::DEFAULT_Z_JET_LAUNCHING, "Distance from black hole where the jet starts [rg]", SEP_COMMA)                                                                                       \
    X(r_initial, double, defaults::DEFAULT_R_INITIAL, "Jet radius when the jet starts [rg]", SEP_COMMA)                                                                                                                  \
    X(z_dissipation, double, defaults::DEFAULT_Z_DISSIPATION, "Distance from black hole where the jet starts to dissipate energy into non-thermal particles [rg]", SEP_COMMA)                                            \
    X(z_max_calculation, double, defaults::DEFAULT_Z_MAX_CALCULATION, "Distance from black hole where the calulation stops [rg]", SEP_COMMA)                                                                             \
    X(electron_temperature_jet_base, double, defaults::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE, "Electron temperature value at z_jet_launching (the jet base) [keV]", SEP_COMMA)                                           \
    X(plasma_beta_jet_base, double, defaults::DEFAULT_PLASMA_BETA_JET_BASE, "Plasma beta value at z_jet_launching (the jet base)", SEP_COMMA)                                                                            \
    X(fraction_nonthermal_electrons, double, defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS, "Fraction of energy in non-thermal electron tail", SEP_COMMA)                                                      \
    X(fraction_nonthermal_protons, double, defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS, "Fraction of energy in non-thermal proton tail", SEP_COMMA)                                                            \
    X(factor_break_electrons, double, defaults::DEFAULT_FACTOR_BREAK_ELECTRONS, "Scaling factor for electron adiabtic timescale", SEP_COMMA)                                                                     \
    X(factor_break_protons, double, defaults::DEFAULT_FACTOR_BREAK_PROTONS, "Scaling factor for proton adiabtic timescale", SEP_COMMA)                                                                           \
    X(factor_max_energy_electrons, double, defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS, "Scaling factor for electron acceleration rate, translating to an effective scaling of the maximum energy", SEP_COMMA) \
    X(factor_max_energy_protons, double, defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS, "Scaling factor for proton acceleration rate, translating to an effective scaling of the maximum energy", SEP_COMMA)       \
    X(index_injected_electrons, double, defaults::DEFAULT_INDEX_INJECTED_ELECTRONS, "Injected electron spectral index (dlogN/dlogE), ie. before cooling", SEP_COMMA)                                             \
    X(index_injected_protons, double, defaults::DEFAULT_INDEX_INJECTED_PROTONS, "Injected proton spectral index (dlogN/dlogE), ie. before cooling", SEP_COMMA)                                                   \
    X(dlgz, double, defaults::DEFAULT_DLGZ, "log. grid spacing log10(z[i+1]) - log10(z[i])", SEP_COMMA) \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, "Regulates print output of the code. 0: No output; 1: Only important warnings; 2: More output; 3: Debugging output", SEP_COMMA) \
    X(cutoff_type, int, defaults::DEFAULT_CUTOFF_TYPE, "Particle cutoff shape. 0: exp(-x), 1: exp(-x^2), 2: sech(x)^2.", )

#define BHJET_PARAMS                                                                                                                                                 \
    X(theta_obs, double, defaults::DEFAULT_THETA_OBS, "Observation angle [degree]", SEP_COMMA)                                                                  \
    X(distance, double, defaults::DEFAULT_DISTANCE, "Distance to source [kpc]", SEP_COMMA)                                                                      \
    X(redshift, double, defaults::DEFAULT_REDSHIFT, "Redshift of source", SEP_COMMA)                                                                            \
    X(include_counterjet, bool, defaults::DEFAULT_INCLUDE_COUNTERJET, "True: Includes the emission of the counterjet; False: Includes only one jet", SEP_COMMA) \
    X(compton_threshold, double, defaults::DEFAULT_COMPTON_THRESHOLD, "Internal rough threshold for Compton emission computation: L_com/L_syn > compton_threshold ? Do calculation. Otheriwse skip.", SEP_COMMA) \
    X(profile_time, bool, defaults::DEFAULT_PROFILE_TIME, "True: measures computation time of each zone; False: Does nothing", SEP_COMMA)                               \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, "Regulates print output of the code. 0: No output; 1: Only important warnings; 2: More output; 3: Debugging output", )


#define TARGET_PHOTON_FIELD_PARAMS  \
    X_NODEF(name, std::string, SEP_COMMA) \
    X(distance, double, defaults::DISTANCE, SEP_COMMA) \
    X(redshift, double, defaults::REDSHIFT, SEP_COMMA) \
    X(add_to_total_flux, bool, defaults::DEFAULT_ADD_TO_TOTAL_FLUX, SEP_COMMA) \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, )

#define TARGET_BLACK_BODY_PARAMS                                \
    X_NODEF(name, std::string, "Name, has to be unique", SEP_COMMA) \
    X(distance, double, defaults::DISTANCE, "Distance to source [kpc]", SEP_COMMA ) \
    X(redshift, double, defaults::REDSHIFT, "Redshift of source", SEP_COMMA) \
    X(luminosity, double, defaults::DEFAULT_LUMINOSTIY, "Total luminosity of black body [erg/s]", SEP_COMMA) \
    X(temperature, double, defaults::DEFAULT_TEMPERATURE, "Temperature of black body [keV]", SEP_COMMA) \
    X(energy_density, double, defaults::DEFAULT_ENERGY_DENSITY, "Total energy density of black body [erg/cm³]", SEP_COMMA) \
    X(add_to_total_flux, bool, defaults::DEFAULT_ADD_TO_TOTAL_FLUX, "True: add the flux to the total model flux, False: only take into account as a target.", SEP_COMMA) \
    X(verbosity_level, size_t, defaults::DEFAULT_VERBOSITY_LEVEL, "Regulates print output of the code. 0: No output; 1: Only important warnings; 2: More output; 3: Debugging output", )

// function to convert the returned std::vector<type> from c++ function "function"
// into a numpy array
#define GET_ARGS_VEC(classtype, function, type) \
    [](classtype &self) {                                 \
        auto vec = self.function(); \
        return py::array_t<type>(vec.size(), vec.data()); }
#define GET_TIMESCALE(classtype, function, type, argname) \
    [](classtype &self, std::vector<double> argname) {                  \
        auto vec = self.function(argname);                 \
        return py::array_t<type>(vec.size(), vec.data()); }, py::arg("argname")


#define GET_NPARRAY_ARG(classtype, function, type, argname, argtype) \
    [](classtype &self, argtype argname) {                  \
        auto vec = self.function(argname);                 \
        return py::array_t<type>(vec.size(), vec.data()); }, py::arg("argname")

#define GET_NPARRAY_2ARG(classtype, function, type, argname, argtype, argname2, argtype2) \
    [](classtype &self, argtype argname, argtype2 argname2) {                  \
        auto vec = self.function(argname, argname2);                 \
        return py::array_t<type>(vec.size(), vec.data()); }, py::arg("argname"), py::arg("argname2")

PYBIND11_MODULE(bhjet, m)
{
    m.doc() = "Jet dynamics simulation module with Kariba backend";

    py::class_<JetDynamics, std::shared_ptr<JetDynamics>> jetdyn(m, "JetDynamics");
    jetdyn.def(py::init<
#define X(NAME, TYPE, DEFAULT, SEPARATOR) TYPE SEPARATOR
                   JETDYNAMICS_PARAMS
#undef X
                   >(),
// py::arg defaults
#define X(NAME, TYPE, DEFAULT, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
               JETDYNAMICS_PARAMS
#undef X
               )
// members
#define X(NAME, TYPE, DEFAULT, SEPARATOR) .def_readwrite(#NAME, &JetDynamics::NAME)
        JETDYNAMICS_PARAMS
#undef X
        ;
    jetdyn.def("compute_jet_dynamics", &JetDynamics::compute_jet_dynamics);
    jetdyn.def("get_z_min_grid", GET_ARGS_VEC(JetDynamics, get_z_min_grid, double), "Get array with z grid of zone start positions [r_g]");
    jetdyn.def("get_z_height_grid", GET_ARGS_VEC(JetDynamics, get_z_height_grid, double), "Get array of zone heights [r_g]");
    jetdyn.def("get_z_center_grid", GET_ARGS_VEC(JetDynamics, get_z_center_grid, double), "Get array of zone center positions [r_g]");
    jetdyn.def("get_radius_grid", GET_ARGS_VEC(JetDynamics, get_radius_grid, double), "Get array of zone widths [r_g]");
    jetdyn.def("get_gamma_grid", GET_ARGS_VEC(JetDynamics, get_gamma_grid, double), "Get array of zone bulk Lorentz factors");
    jetdyn.def("get_beta_grid", GET_ARGS_VEC(JetDynamics, get_beta_grid, double), "Get array of zone bulk velocities [speed of light]");
    jetdyn.def("get_beta_gamma_grid", GET_ARGS_VEC(JetDynamics, get_beta_gamma_grid, double), "Get array of zone bulk speeds (beta*gamma)");
    jetdyn.def("get_magnetic_field_grid", GET_ARGS_VEC(JetDynamics, get_magnetic_field_grid, double), "Get array of zone magnetic field strnegths [G]");
    jetdyn.def("get_electron_density_grid", GET_ARGS_VEC(JetDynamics, get_electron_density_grid, double), "Get array of zone electron densities [1/cm³]");
    jetdyn.def("get_electron_temperature_grid", GET_ARGS_VEC(JetDynamics, get_electron_temperature_grid, double), "Get array of zone electron temperatures [keV]");
    jetdyn.def("get_proton_temperature_grid", GET_ARGS_VEC(JetDynamics, get_proton_temperature_grid, double), "Get array of zone proton temperatures [keV]");
    jetdyn.def("get_fraction_nonthermal_electrons_grid", GET_ARGS_VEC(JetDynamics, get_fraction_nonthermal_electrons_grid, double), "Get array of non-thermal electron fractions");
    jetdyn.def("get_fraction_nonthermal_protons_grid", GET_ARGS_VEC(JetDynamics, get_fraction_nonthermal_protons_grid, double), "Get array of non-thermal proton fractions");
    jetdyn.def("get_factor_break_electrons_grid", GET_ARGS_VEC(JetDynamics, get_factor_break_electrons_grid, double), "Get array of scaling factors for electron adiabtic timescale");
    jetdyn.def("get_factor_break_protons_grid", GET_ARGS_VEC(JetDynamics, get_factor_break_protons_grid, double), "Get array of scaling factors for proton adiabtic timescale");
    jetdyn.def("get_factor_max_energy_electrons_grid", GET_ARGS_VEC(JetDynamics, get_factor_max_energy_electrons_grid, double), "Get array of scaling factors for electron acceleration rate, translating to an effective scaling of the maximum energy");
    jetdyn.def("get_factor_max_energy_protons_grid", GET_ARGS_VEC(JetDynamics, get_factor_max_energy_protons_grid, double), "Get array of  scaling factors for proton acceleration rate, translating to an effective scaling of the maximum energy");
    jetdyn.def("get_index_injected_electrons_grid", GET_ARGS_VEC(JetDynamics, get_index_injected_electrons_grid, double), "Get array of injected electron spectral indices");
    jetdyn.def("get_index_injected_protons_grid", GET_ARGS_VEC(JetDynamics, get_index_injected_protons_grid, double), "Get array of injected proton spectral indices");
    jetdyn.def("get_cutoff_type_grid", GET_ARGS_VEC(JetDynamics, get_cutoff_type_grid, int), "Get array of particle cutoff prescription types");
    jetdyn.def("info", &JetDynamics::info);

    py::class_<BLJet, JetDynamics, std::shared_ptr<BLJet>> bljet(m, "BLJet");
    bljet.def(py::init<
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) TYPE SEPARATOR
                  BLJET_PARAMS
#undef X
                  >(),
// py::arg defaults
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
              BLJET_PARAMS
#undef X
              )
// members
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) .def_readwrite(#NAME, &BLJet::NAME, DOC)
        BLJET_PARAMS
#undef X
        ;
    bljet.def_readonly("r_g", &BLJet::r_g, "Gravitational Radius [cm]");
    bljet.def_readonly("eddington_luminosity", &BLJet::eddington_luminosity, "Eddington Luminosity [erg/s]");
    bljet.def("compute_jet_dynamics", &BLJet::compute_jet_dynamics, "Computes the physical quantities along the jet");

    py::class_<IsoJet, JetDynamics, std::shared_ptr<IsoJet>> isojet(m, "IsoJet");
    isojet.def(py::init<
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) TYPE SEPARATOR
                  ISOJET_PARAMS
#undef X
                  >(),
// py::arg defaults
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
              ISOJET_PARAMS
#undef X
              )
// members
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) .def_readwrite(#NAME, &IsoJet::NAME, DOC)
        ISOJET_PARAMS
#undef X
        ;
    isojet.def_readonly("r_g", &IsoJet::r_g, "Gravitational Radius [cm]");
    isojet.def_readonly("eddington_luminosity", &IsoJet::eddington_luminosity, "Eddington Luminosity [erg/s]");
    isojet.def("compute_jet_dynamics", &IsoJet::compute_jet_dynamics, "Computes the physical quantities along the jet");


    py::class_<TargetPhotonField, std::shared_ptr<TargetPhotonField>> target(m, "TargetPhotonField");
    target.def("get_observed_energy", GET_ARGS_VEC(TargetPhotonField, get_observed_energy, double), "Returns the energy grid of the target photons [erg].");
    target.def("get_observed_energy_flux", GET_ARGS_VEC(TargetPhotonField, get_observed_energy_flux, double), "Returns the energy flux on the energy grid of the target photons [erg/(cm²s)].");


    py::class_<TargetBlackBody, TargetPhotonField, std::shared_ptr<TargetBlackBody>> target_bb(m, "TargetBlackBody");
    target_bb.def(py::init<
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR)       TYPE SEPARATOR
#define X_NODEF(NAME, TYPE, DOC, SEPARATOR)          TYPE SEPARATOR
        TARGET_BLACK_BODY_PARAMS
#undef X
#undef X_NODEF
                    >(),
    // py::arg defaults
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR)   py::arg(#NAME) = DEFAULT SEPARATOR
#define X_NODEF(NAME, TYPE, DOC, SEPARATOR)          py::arg(#NAME) SEPARATOR
        TARGET_BLACK_BODY_PARAMS
#undef X
#undef X_NODEF
    )
    // members
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR)   .def_readwrite(#NAME, &TargetBlackBody::NAME)
#define X_NODEF(NAME, TYPE, DOC, SEPARATOR)          .def_readwrite(#NAME, &TargetBlackBody::NAME)
        TARGET_BLACK_BODY_PARAMS
#undef X
#undef X_NODEF
    ;
    target_bb.def("update_observed_flux", &TargetBlackBody::update_observed_flux, "Recalculate the oberved flux arrays.");
    target_bb.def("get_target_energy_grid_and_density",
        [](TargetBlackBody &self, double z, double bulk_momentum, double theta_obs, double min_energy, double max_energy) {
            auto [a1, a2] = self.get_target_energy_grid_and_density(z, bulk_momentum, theta_obs, min_energy, max_energy);
            return py::make_tuple(
                py::array_t<double>(a1.size(), a1.data()),
                py::array_t<double>(a2.size(), a2.data())
            );
        },
        py::arg("z"), py::arg("bulk_momentum"), py::arg("theta_obs"), py::arg("min_energy"), py::arg("max_energy"),
        "Returns a tuple of the comoving energy grid [erg] with the corresponding energy density points [erg/cm³]"
    );


    py::class_<RadiationZone> radzone(m, "RadiationZone");
    // constructor
    radzone.def(py::init<
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) TYPE SEPARATOR
                    RADIATIONZONE_PARAMS
#undef X
                    >(),
// py::arg defaults
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
                RADIATIONZONE_PARAMS
#undef X
                ,
                "Sets up a zone with the given parameters. Next add target fields, compute the particle spectra and then the emission.")
// members
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) .def_readwrite(#NAME, &RadiationZone::NAME, DOC)
        RADIATIONZONE_PARAMS
#undef X
        ;
    radzone.def_readonly("beta_bulk", &RadiationZone::beta_bulk, "Bulk beta = speed/(speed of light) of the zone");
    radzone.def_readonly("gamma_bulk", &RadiationZone::gamma_bulk, "Bulk Lorentz factor of the zone");
    radzone.def_readonly("doppler_factor_bulk", &RadiationZone::doppler_factor_bulk, "Bulk Doppler factor of the zone");
    radzone.def("add_target_photon_field", &RadiationZone::add_target_photon_field, py::arg("energies"), py::arg("energy_densities"), py::arg("name"), "Takes as arguments the energies [erg] and energy densities [erg/cm³] of a target radiation field. Run before compute_particles() and compute_radiation().");
    radzone.def("remove_target_photon_field", &RadiationZone::remove_target_photon_field, py::arg("name"), "Removes the target photon field with name.");
    radzone.def("compute_particles", &RadiationZone::compute_particles, py::arg("reset_photon_targets")=true, "Computes the steady-state particle spectra (electrons). Run before compute_radiation().");
    radzone.def("compute_radiation", py::overload_cast<>(&RadiationZone::compute_radiation), "Computes the radiation from the particles (photons). Uses the default energy grid.");
    radzone.def("compute_radiation", py::overload_cast<const std::vector<double>>(&RadiationZone::compute_radiation), "Computes the radiation from the particles (photons). Takes as an argument the observed energy grid [erg].");
    // densities
    radzone.def("get_electron_momentum_number_density", GET_ARGS_VEC(RadiationZone, get_electron_momentum_number_density, double), "Get array with comoving electron momentum number density grid dN/dlnp [1/cm³]");
    radzone.def("get_electron_gamma_number_density", GET_ARGS_VEC(RadiationZone, get_electron_gamma_number_density, double), "Get array with comoving electron energy number density grid dN/dlnE [1/cm³]");
    radzone.def("get_electron_momentum_grid", GET_ARGS_VEC(RadiationZone, get_electron_momentum_grid, double), "Get array with comoving electron momentum grid [cm g / s]");
    radzone.def("get_electron_gamma_grid", GET_ARGS_VEC(RadiationZone, get_electron_gamma_grid, double), "Get array with comoving electron Lorentz factor grid");
    // time scales
    radzone.def("get_timescale_electron_cyclosyn", GET_TIMESCALE(RadiationZone, get_timescale_electron_cyclosyn, double, momentum), "Get array with comoving cyclosynchrotron cooling timescale [s]");
    radzone.def("get_timescale_electron_adiabatic", GET_TIMESCALE(RadiationZone, get_timescale_electron_adiabatic, double, momentum), "Get array with comoving adiabatic cooling timescale [s]");
    radzone.def("get_timescale_electron_compton_thomson", GET_TIMESCALE(RadiationZone, get_timescale_electron_compton_thomson, double, momentum), "Get array with comoving Compton cooling timescale in the Thomson approximation (used in the code) [s]");
    radzone.def("get_timescale_electron_compton", GET_TIMESCALE(RadiationZone, get_timescale_electron_compton, double, momentum), "Get array with comoving Compton cooling timescale (including Klein-Nishina effects) [s]");
    radzone.def("get_timescale_electron_acceleration", GET_TIMESCALE(RadiationZone, get_timescale_electron_acceleration, double, momentum), "Get array with comoving acceleration timescale [s]");
    radzone.def("get_timescale_photon_cyclosyn_selfabsorption", GET_TIMESCALE(RadiationZone, get_timescale_photon_cyclosyn_selfabsorption, double, momentum), "Get array with comoving photon cyclosyn-self-absorption timescale [s]");
    radzone.def("get_timescale_photon_escape", GET_TIMESCALE(RadiationZone, get_timescale_photon_escape, double, momentum), "Get array with comoving photon escape timescale [s]");
    // radzone.def("get_timescale_photon_cyclosyn_selfabsorption_array", GET_ARGS_VEC(RadiationZone, get_timescale_photon_cyclosyn_selfabsorption_array, double), " [s]");
    // radzone.def("get_timescale_photon_cyclosyn_selfabsorption_energy", GET_ARGS_VEC(RadiationZone, get_timescale_photon_cyclosyn_selfabsorption_energy, double), " [s]");
    // characteristic energies
    radzone.def("get_electron_max_momentum", &RadiationZone::get_electron_max_momentum, "Get maximum electron momentum from comparing acceleration with cooling times [cm g / s]");
    radzone.def("get_electron_break_momentum", &RadiationZone::get_electron_break_momentum, "Get break electron momentum from comparing adiabatic with cyclosyn./Compton cooling times [cm g / s]");

    // target fields
    radzone.def("get_total_photon_target_energy", GET_ARGS_VEC(RadiationZone, get_total_photon_target_energy, double), "Get array with target photon energy grid [erg]");
    radzone.def("get_total_photon_target_energy_density", GET_ARGS_VEC(RadiationZone, get_total_photon_target_energy_density, double), "Get array with total photon target field [erg/cm³]");
    radzone.def("get_photon_target_energy", GET_NPARRAY_ARG(RadiationZone, get_photon_target_energy, double, name, std::string), "Get array with target photon energy grid of a component with name [erg]");
    radzone.def("get_photon_target_energy_density", GET_NPARRAY_ARG(RadiationZone, get_photon_target_energy_density, double, name, std::string), "Get array with total photon target field of a component with name [erg/cm³]");

    // observed radiation spectra
    radzone.def("get_observed_photon_energy_grid_electron_cyclosyn", GET_ARGS_VEC(RadiationZone, get_observed_photon_energy_grid_electron_cyclosyn, double), "Get array with energy grid of electron cyclosynchrotron luminosity/flux [erg]");
    radzone.def("get_observed_photon_energy_grid_electron_compton", GET_ARGS_VEC(RadiationZone, get_observed_photon_energy_grid_electron_compton, double), "Get array with energy grid of electron Compton luminosity/flux [erg]");
    radzone.def("get_observed_photon_energy_grid_total", GET_ARGS_VEC(RadiationZone, get_observed_photon_energy_grid_total, double), "Get array with energy grid of total luminosity/flux [erg]");
    radzone.def("get_observed_photon_flux_electron_cyclosyn", GET_ARGS_VEC(RadiationZone, get_observed_photon_flux_electron_cyclosyn, double), "Get array with electron cyclosynchrotron flux [1/(cm²s)]");
    radzone.def("get_observed_photon_flux_electron_compton", GET_ARGS_VEC(RadiationZone, get_observed_photon_flux_electron_compton, double), "Get array with electron Compton flux [1/(cm²s)]");
    radzone.def("get_observed_photon_flux_total", GET_ARGS_VEC(RadiationZone, get_observed_photon_flux_total, double), "Get array with total flux [1/(cm²s)]");
    // These rates below are intermediate results
    // radzone.def("get_observed_photon_luminosity_electron_cyclosyn", GET_ARGS_VEC(RadiationZone, get_observed_photon_luminosity_electron_cyclosyn, double), "Get array with electron cyclosynchrotron photon rate [erg/(sHz)]");
    // radzone.def("get_observed_photon_luminosity_electron_compton", GET_ARGS_VEC(RadiationZone, get_observed_photon_luminosity_electron_compton, double), "Get array with electron Compton photon rate [erg/(sHz)]");
    // radzone.def("get_observed_photon_luminosity_total", GET_ARGS_VEC(RadiationZone, get_observed_photon_luminosity_total, double), "Get array with total photon rate [erg/(sHz)]");

    radzone.def("compton_calculation_necessary", &RadiationZone::compton_calculation_necessary, "True: Compton emission is calculated; False: Compton emission is not calculated.");
    radzone.def("get_computation_times", GET_ARGS_VEC(RadiationZone, get_computation_times, double), "Get array with computation times [ns]");

    py::class_<BHJet> bhjet(m, "BHJet");
    bhjet.def(py::init<
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) TYPE SEPARATOR
                  BHJET_PARAMS
#undef X
                  >(),
// py::arg defaults
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) py::arg(#NAME) = DEFAULT SEPARATOR
              BHJET_PARAMS
#undef X
              ,
              "Sets up a jet with multiple RadiationZones with the given radiation parameters as well as \
            properties derived from a JetDynamics object. Next init the jet dynamics, add target fields, \
            compute the full jet emission.")
// members
#define X(NAME, TYPE, DEFAULT, DOC, SEPARATOR) .def_readwrite(#NAME, &BHJet::NAME, DOC)
        BHJET_PARAMS
#undef X
        ;
    bhjet.def("init_jet_dynamics", &BHJet::init_jet_dynamics, "Assign a class of type JetDynamics \
              or a derivative (e.g. BLJet) to be used for the dynamics and properties of the jet");
    bhjet.def("get_jet_dynamics", &BHJet::get_jet_dynamics, "Returns the JetDynamics pointer");
    // bhjet.def("add_target_constant_black_body", &BHJet::add_target_constant_black_body,
    //           py::arg("luminosity"), py::arg("temperature"), py::arg("energy_density"), py::arg("name"),
    //           "Add a target black body radiation field, at rest in the black hole frame indexed with a given name. \
    //           Total luminosity [erg/s] is added to the total emission, temperature [keV] and energy_density [erg/cm³] \
    //           are boosted to each RadiationZone.");
    // bhjet.def("add_target_constant_bulge", &BHJet::add_target_constant_bulge,
    //           py::arg("luminosity"), py::arg("temperature"), py::arg("radius"), py::arg("name"),
    //           "Add a target black body radiation field, at rest in the black hole frame indexed with a given name. \
    //           Total luminosity [erg/s] is added to the total emission, temperature [keV] and radius (converted to a \
    //           homogeneous energy_density) [erg/cm³]  are boosted to each RadiationZone.");
    // bhjet.def("add_target_cmb", &BHJet::add_target_cmb,
    //           "Add a target black body radiation field with the CMB properties called CMB, at rest in the black hole \
    //           frame indexed with a given name.");
              
    // bhjet.def("remove_target_black_body", &BHJet::remove_target_black_body, "Removes a target photon field with a given name.");
    bhjet.def("add_target_photon_field", &BHJet::add_target_photon_field, "Add a target photon field.");
    bhjet.def("remove_target_photon_field", &BHJet::remove_target_photon_field, "Removes a target photon field.");
    bhjet.def("clear_target_photon_fields", &BHJet::clear_target_photon_fields, "Removes all target photon fields");

    bhjet.def("compute_full_jet", &BHJet::compute_full_jet,
              py::arg("photon_energy_grid"),
              "Iterates over each RadiationZone, determines its properties and performs the emission modelling. \
              Takes the photon_energy_grid to compute the emission spectra");
    // bhjet.def_readonly("radiation_zones", &BHJet::radiation_zones, "Array with the RadiationZone objects");
    bhjet.def_property_readonly(
        "radiation_zones",
        [](BHJet &self) -> std::vector<RadiationZone> &
        {
            return self.radiation_zones;
        },
        py::return_value_policy::reference_internal,
        "Array with the RadiationZone objects");

    bhjet.def("get_observed_photon_energy_grid", GET_ARGS_VEC(BHJet, get_observed_photon_energy_grid, double),
              "Get array with observed photon energies [erg]");
    bhjet.def("get_observed_photon_flux_total", GET_ARGS_VEC(BHJet, get_observed_photon_flux_total, double),
              "Get array with total observed photon flux [1/(cm²s)]");

    bhjet.def("get_observed_photon_integrated_flux_total",
              GET_NPARRAY_2ARG(BHJet, get_observed_photon_integrated_flux_total, double, z_min, double, z_max, double),
              "Get array with with total observed photon flux [1/(cm²s)] integrated between z_min [cm] and z_max [cm]");
    bhjet.def("get_observed_photon_integrated_flux_electron_cyclosyn",
              GET_NPARRAY_2ARG(BHJet, get_observed_photon_integrated_flux_electron_cyclosyn, double, z_min, double, z_max, double),
              "Get array with with electron cyclosyn. observed photon flux [1/(cm²s)] integrated between z_min [cm] and z_max [cm]");
    bhjet.def("get_observed_photon_integrated_flux_electron_compton",
              GET_NPARRAY_2ARG(BHJet, get_observed_photon_integrated_flux_electron_compton, double, z_min, double, z_max, double),
              "Get array with with electron compton observed photon flux [1/(cm²s)] integrated between z_min [cm] and z_max [cm]");


    // bhjet.def("get_target_black_body_temperature", &BHJet::get_target_black_body_temperature, "Get the temperature of the black body target with the given name [keV]");
    // bhjet.def("get_target_black_body_energy_density", &BHJet::get_target_black_body_energy_density, "Get the integrated energy density of the black body target with the given name [erg/cm³]");
    // bhjet.def("get_target_black_body_luminosity", &BHJet::get_target_black_body_luminosity, "Get the luminosity of the black body target with the given name [erg/s]");
    // bhjet.def("set_target_black_body_temperature", &BHJet::set_target_black_body_temperature, "Set the temperature of the black body target with the given name [keV]");
    // bhjet.def("set_target_black_body_energy_density", &BHJet::set_target_black_body_energy_density, "Set the integrated energy density of the black body target with the given name [erg/cm³]");
    // bhjet.def("set_target_black_body_luminosity", &BHJet::set_target_black_body_luminosity, "Set the luminosity of the black body target with the given name [erg/s]");
    // bhjet.def("get_observed_photon_energy_grid_black_body", GET_NPARRAY_ARG(BHJet, get_observed_photon_energy_grid_black_body, double, name, std::string),
    //           "Get array with observed photon energies for the black body with a given name [erg]");
    // bhjet.def("get_observed_photon_flux_black_body", GET_NPARRAY_ARG(BHJet, get_observed_photon_flux_black_body, double, name, std::string),
    //           "Get array with total observed photon flux for the black body with a given name [1/(cm²s)]");
    bhjet.def("get_computation_times", GET_ARGS_VEC(BHJet, get_computation_times, double), "Get array with computation times [ns]");
}
