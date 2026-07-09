#pragma once
#include <cstddef>

namespace bhjet {
    namespace defaults {

        static constexpr size_t DEFAULT_VERBOSITY_LEVEL = 1;

        // BHJet

        // TargetPhotonField
        static constexpr double DISTANCE = 1e3;
        static constexpr double REDSHIFT = 0.0;
        static constexpr bool DEFAULT_ADD_TO_TOTAL_FLUX = true;
        // Target Black Body
        static constexpr double DEFAULT_LUMINOSTIY = 1e39;
        static constexpr double DEFAULT_TEMPERATURE = 1e-3; // keV
        static constexpr double DEFAULT_ENERGY_DENSITY = 1e-9; // erg/cm³


        // Jet Dynamics
        static constexpr size_t DEFAULT_N_ZONES = 100;
        // BLJet
        static constexpr double DEFAULT_MASS_BH = 1e9;
        static constexpr double DEFAULT_JET_POWER_EDDINGTON = 1e-5;
        static constexpr double DEFAULT_Z_JET_LAUNCHING = 2;
        static constexpr double DEFAULT_R_INITIAL = 3;
        static constexpr double DEFAULT_Z_END_OF_ACCELERATION = 1e5;
        static constexpr double DEFAULT_Z_DISSIPATION = 1e2;
        static constexpr double DEFAULT_Z_MAX_CALCULATION = 1e6;
        static constexpr double DEFAULT_SIGMA_FINAL = 1;
        static constexpr double DEFAULT_GAMMA_FINAL = 15;
        static constexpr double DEFAULT_ELECTRON_TEMPERATURE_JET_BASE = 1e3;
        static constexpr double DEFAULT_GAMMA_ACCELERATION_EXPONENT = 0.5;
        static constexpr double DEFAULT_GAMMA_DECELERATION_EXPONENT = 0.0;
        static constexpr double DEFAULT_OPENING_ANGLE_CONSTANT = 0.15;
        static constexpr bool DEFAULT_CALC_PAIR_CONTENT_FROM_PLASMA_BETA = false;
        static constexpr double DEFAULT_PLASMA_BETA_JET_BASE = 1;
        static constexpr double DEFAULT_DLGZ = 0.1;


        // Radiation zone
        static constexpr double DEFAULT_MAGNETIC_FIELD = 1.0;
        static constexpr double DEFAULT_RADIUS = 1e15;
        static constexpr double DEFAULT_HEIGHT = 1e15;
        static constexpr const char *DEFAULT_GEOMETRY = "sphere";
        static constexpr double DEFAULT_BULK_MOMENTUM = 1.0;
        static constexpr double DEFAULT_THETA_OBS = 0.0;
        static constexpr double DEFAULT_DISTANCE = 1e3;
        static constexpr double DEFAULT_REDSHIFT = 0.0;
        static constexpr double DEFAULT_ELECTRON_NUMBER_DENSITY = 1;
        static constexpr double DEFAULT_PROTON_NUMBER_DENSITY = 1;
        static constexpr double DEFAULT_ELECTRON_TEMPERATURE = 1e3;
        static constexpr double DEFAULT_PROTON_TEMPERATURE = 1e3;
        static constexpr double DEFAULT_FRACTION_NONTHERMAL_ELECTRONS = 0.1;
        static constexpr double DEFAULT_FRACTION_NONTHERMAL_PROTONS = 0.1;
        static constexpr double DEFAULT_FACTOR_BREAK_ELECTRONS = 1.;
        static constexpr double DEFAULT_FACTOR_BREAK_PROTONS = 1.;
        static constexpr double DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS = 1.;
        static constexpr double DEFAULT_FACTOR_MAX_ENERGY_PROTONS = 1.;
        static constexpr double DEFAULT_INDEX_INJECTED_ELECTRONS = 2.;
        static constexpr double DEFAULT_INDEX_INJECTED_PROTONS = 2.;
        static constexpr int DEFAULT_CUTOFF_TYPE = 0;
        static constexpr bool DEFAULT_INCLUDE_COUNTERJET = true;
        static constexpr bool DEFAULT_FORCE_COMPTON_CALCULATION = false;
        static constexpr bool DEFAULT_COMPTON_SWITCH = true;
        static constexpr double DEFAULT_COMPTON_THRESHOLD = 1e-2;
        static constexpr bool DEFAULT_PROFILE_TIME = false;
    }
}
