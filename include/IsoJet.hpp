#pragma once
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include "kariba/Thermal.hpp"
#include "default_values.hpp"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

namespace bhjet
{


    class IsoJet : public JetDynamics
    {
    public:
        // Member variables
        double mass_bh = defaults::DEFAULT_MASS_BH;
        double jet_power_eddington = defaults::DEFAULT_JET_POWER_EDDINGTON;
        double z_jet_launching = defaults::DEFAULT_Z_JET_LAUNCHING;
        double r_initial = defaults::DEFAULT_R_INITIAL;
        double z_dissipation = defaults::DEFAULT_Z_DISSIPATION;
        double z_max_calculation = defaults::DEFAULT_Z_MAX_CALCULATION;
        double electron_temperature_jet_base = defaults::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE;
        double plasma_beta_jet_base = defaults::DEFAULT_PLASMA_BETA_JET_BASE;
        double fraction_nonthermal_electrons = defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS;
        double fraction_nonthermal_protons = defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS;
        double factor_break_electrons = defaults::DEFAULT_FACTOR_BREAK_ELECTRONS;
        double factor_break_protons = defaults::DEFAULT_FACTOR_BREAK_PROTONS;
        double factor_max_energy_electrons = defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS;
        double factor_max_energy_protons = defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS;
        double index_injected_electrons = defaults::DEFAULT_INDEX_INJECTED_ELECTRONS;
        double index_injected_protons = defaults::DEFAULT_INDEX_INJECTED_PROTONS;
        int cutoff_type = defaults::DEFAULT_CUTOFF_TYPE;
        double dlgz = defaults::DEFAULT_DLGZ;

        // internal variables
        double eddington_luminosity, r_g, zmin;

        // hardcoded values
        // size_t n_bins_speed = 54;
        // double jet_opening_constant = 0.15;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        IsoJet(
            double mass_bh_ = defaults::DEFAULT_MASS_BH,
            double jet_power_eddington_ = defaults::DEFAULT_JET_POWER_EDDINGTON,
            double z_jet_launching_ = defaults::DEFAULT_Z_JET_LAUNCHING,
            double r_initial_ = defaults::DEFAULT_R_INITIAL,
            double z_dissipation_ = defaults::DEFAULT_Z_DISSIPATION,
            double z_max_calculation_ = defaults::DEFAULT_Z_MAX_CALCULATION,
            double electron_temperature_jet_base_ = defaults::DEFAULT_ELECTRON_TEMPERATURE_JET_BASE,
            double plasma_beta_jet_base_ = defaults::DEFAULT_PLASMA_BETA_JET_BASE,
            double fraction_nonthermal_electrons_ = defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
            double fraction_nonthermal_protons_ = defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS,
            double factor_break_electrons_ = defaults::DEFAULT_FACTOR_BREAK_ELECTRONS,
            double factor_break_protons_ = defaults::DEFAULT_FACTOR_BREAK_PROTONS,
            double factor_max_energy_electrons_ = defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
            double factor_max_energy_protons_ = defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
            double index_injected_electrons_ = defaults::DEFAULT_INDEX_INJECTED_ELECTRONS,
            double index_injected_protons_ = defaults::DEFAULT_INDEX_INJECTED_PROTONS,
            double dlgz_ = defaults::DEFAULT_DLGZ,
            size_t verbosity_level_ = defaults::DEFAULT_VERBOSITY_LEVEL,
            int cutoff_type_ = defaults::DEFAULT_CUTOFF_TYPE)
            : JetDynamics(defaults::DEFAULT_N_ZONES, verbosity_level_),
              mass_bh(mass_bh_),
              jet_power_eddington(jet_power_eddington_),
              z_jet_launching(z_jet_launching_), r_initial(r_initial_),
              z_dissipation(z_dissipation_), z_max_calculation(z_max_calculation_),
              electron_temperature_jet_base(electron_temperature_jet_base_),
              plasma_beta_jet_base(plasma_beta_jet_base_),
              fraction_nonthermal_electrons(fraction_nonthermal_electrons_),
              fraction_nonthermal_protons(fraction_nonthermal_protons_),
              factor_break_electrons(factor_break_electrons_), factor_break_protons(factor_break_protons_),
              factor_max_energy_electrons(factor_max_energy_electrons_),
              factor_max_energy_protons(factor_max_energy_protons_),
              index_injected_electrons(index_injected_electrons_), index_injected_protons(index_injected_protons_),
              cutoff_type(cutoff_type_),
              dlgz(dlgz_)
        {
        }
        ~IsoJet() override = default;

        gsl_spline *spline_speed = nullptr;
        gsl_interp_accel *spline_speed_accel = nullptr;

        jet_dynpars jet_dyn;    // structure with jet dynamical parameters
        jet_enpars nozzle_ener; // structure with jet energetic parameters

        void calc_velocity_profile_iso_jet();
        void calc_nozzle_energetics_equipartition();

        void calc_zone_properties(size_t i);

        void compute_jet_dynamics() override;
    };

} // namespace bhjet
