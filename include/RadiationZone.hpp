#pragma once
#include <string>
#include <vector>
#include <memory>
#include "utils.hpp"
#include "kariba/Mixed.hpp"
#include "kariba/Thermal.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/Bknpower.hpp"
#include "kariba/BBody.hpp"
#include "kariba/ShSDisk.hpp"
#include "TargetPhotonField.hpp"
#include "kariba/constants.hpp"
#include "default_values.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{

    // struct TargetBlackBody
    // {
    //     double temperature, energy_density;
    //     std::string name;
    //     TargetBlackBody(double t, double u, std::string n)
    //         : temperature(t), energy_density(u), name(n) {}
    // };
    // struct TargetDisk
    // {
    //     double mass_bh, inner_radius, outer_radius, luminosity, inclination;
    // };

    class RadiationZone
    {
    public:
        ~RadiationZone();

        std::vector<std::string> target_component_names;
        std::vector<std::vector<double>> 
            target_component_radiation_energy_density, 
            target_component_radiation_energy_grid;
        // total integral (for Thomson cooling)
        double radiation_energy_density;
        // total vectors for convenience on same grid
        std::vector<double> total_target_radiation_energy_density, total_target_radiation_energy_grid;


        void add_target_photon_field(std::vector<double> energies, std::vector<double> energy_densities, std::string name);
        void remove_target_photon_field(std::string name);
        // radiation targets per zone, for total luminosity track one for each in BHJet class
        // void add_target_black_body(double temperature, double energy_density, std::string name);
        // void add_target_disk(double Mbh, double inner_radius, double outer_radius, double luminosity, double inclination);
        // void add_target_field(std::vector<double> target_energy, std::vector<double> target_array);
        
        void compute_particles(bool reset_photon_targets=true);
        
        void compute_radiation();
        void compute_radiation(std::vector<double> obs_energy_grid);


        // ----------------------------
        // Member variables
        // ----------------------------
        double magnetic_field,
            radius,
            height,
            bulk_momentum,
            theta_obs,
            distance,
            redshift,
            electron_number_density, proton_number_density,
            electron_temperature, proton_temperature,
            fraction_nonthermal_electrons, fraction_nonthermal_protons,
            factor_break_electrons, factor_break_protons,
            factor_max_energy_electrons, factor_max_energy_protons,
            index_injected_electrons, index_injected_protons;
        int cutoff_type;
        std::string geometry;
        bool include_counterjet = defaults::DEFAULT_INCLUDE_COUNTERJET;
        bool force_compton_calculation = defaults::DEFAULT_FORCE_COMPTON_CALCULATION;
        bool compton_switch = defaults::DEFAULT_COMPTON_SWITCH;
        double compton_threshold = defaults::DEFAULT_COMPTON_THRESHOLD;
        bool profile_time = defaults::DEFAULT_PROFILE_TIME;
        size_t verbosity_level = defaults::DEFAULT_VERBOSITY_LEVEL;

        // momentum grid
        size_t n_bins_e, n_bins_p;

        kariba::Thermal electrons_thermal;
        kariba::Bknpower electrons_bpl;
        kariba::Mixed electrons_mixed;
        kariba::Powerlaw electrons_pl;

        gsl_spline *spline_electrons;
        gsl_interp_accel *spline_electrons_accel;

        gsl_spline *spline_electrons_derivative;
        gsl_interp_accel *spline_electrons_derivative_accel;

        size_t syn_res = 10;
        size_t com_res = 6;

        double doppler_factor_bulk, beta_bulk, gamma_bulk;

        // std::vector<TargetBlackBody> target_vector_blackbody;
        // std::vector<TargetDisk> target_vector_disk;
        // std::vector<TargetPhotonField> target_vector_photons;

        std::vector<double> computation_times;

        // ----------------------------
        // Constructor with defaults
        // ----------------------------
        RadiationZone(
            double magnetic_field_ = defaults::DEFAULT_MAGNETIC_FIELD,
            double radius_ = defaults::DEFAULT_RADIUS,
            double height_ = defaults::DEFAULT_HEIGHT,
            std::string geometry_ = defaults::DEFAULT_GEOMETRY,
            double bulk_momentum_ = defaults::DEFAULT_BULK_MOMENTUM,
            double theta_obs_ = defaults::DEFAULT_THETA_OBS,
            double distance_ = defaults::DEFAULT_DISTANCE,
            double redshift_ = defaults::DEFAULT_REDSHIFT,
            double electron_number_density_ = defaults::DEFAULT_ELECTRON_NUMBER_DENSITY,
            double proton_number_density_ = defaults::DEFAULT_PROTON_NUMBER_DENSITY,
            double electron_temperature_ = defaults::DEFAULT_ELECTRON_TEMPERATURE,
            double proton_temperature_ = defaults::DEFAULT_PROTON_TEMPERATURE,
            double fraction_nonthermal_electrons_ = defaults::DEFAULT_FRACTION_NONTHERMAL_ELECTRONS,
            double fraction_nonthermal_protons_ = defaults::DEFAULT_FRACTION_NONTHERMAL_PROTONS,
            double factor_break_electrons_ = defaults::DEFAULT_FACTOR_BREAK_ELECTRONS,
            double factor_break_protons_ = defaults::DEFAULT_FACTOR_BREAK_PROTONS,
            double factor_max_energy_electrons_ = defaults::DEFAULT_FACTOR_MAX_ENERGY_ELECTRONS,
            double factor_max_energy_protons_ = defaults::DEFAULT_FACTOR_MAX_ENERGY_PROTONS,
            double index_injected_electrons_ = defaults::DEFAULT_INDEX_INJECTED_ELECTRONS,
            double index_injected_protons_ = defaults::DEFAULT_INDEX_INJECTED_PROTONS,
            bool include_counterjet_ = defaults::DEFAULT_INCLUDE_COUNTERJET,
            bool force_compton_calculation_ = defaults::DEFAULT_FORCE_COMPTON_CALCULATION,
            bool compton_switch_ = defaults::DEFAULT_COMPTON_SWITCH,
            double compton_threshold_ = defaults::DEFAULT_COMPTON_THRESHOLD,
            bool profile_time_ = defaults::DEFAULT_PROFILE_TIME,
            size_t verbosity_level_ = defaults::DEFAULT_VERBOSITY_LEVEL,
            int cutoff_type_ = defaults::DEFAULT_CUTOFF_TYPE)
            : magnetic_field(magnetic_field_), radius(radius_), height(height_), geometry(geometry_),
              bulk_momentum(bulk_momentum_), theta_obs(theta_obs_), distance(distance_),
              redshift(redshift_), electron_number_density(electron_number_density_),
              proton_number_density(proton_number_density_),
              electron_temperature(electron_temperature_), proton_temperature(proton_temperature_),
              fraction_nonthermal_electrons(fraction_nonthermal_electrons_),
              fraction_nonthermal_protons(fraction_nonthermal_protons_),
              factor_break_electrons(factor_break_electrons_), factor_break_protons(factor_break_protons_),
              factor_max_energy_electrons(factor_max_energy_electrons_),
              factor_max_energy_protons(factor_max_energy_protons_),
              index_injected_electrons(index_injected_electrons_), index_injected_protons(index_injected_protons_),
              cutoff_type(cutoff_type_),
              include_counterjet(include_counterjet_),
              force_compton_calculation(force_compton_calculation_), compton_switch(compton_switch_),
              compton_threshold(compton_threshold_),
              profile_time(profile_time_),
              verbosity_level(verbosity_level_),
              n_bins_e(100), n_bins_p(100), radiation_energy_density(0.),
              target_component_names(), target_component_radiation_energy_grid(), 
              target_component_radiation_energy_density(),
              electrons_thermal(kariba::Thermal(0)), electrons_mixed(kariba::Mixed(0)),
              electrons_bpl(kariba::Bknpower(0)), electrons_pl(kariba::Powerlaw(0)),
              spline_electrons(nullptr), spline_electrons_accel(nullptr),
              spline_electrons_derivative(nullptr), spline_electrons_derivative_accel(nullptr),
              computation_times(3, 0.) //, additional_target_field_energy(), additional_target_field_energy_density()
        {
            gamma_bulk = pow(1 + bulk_momentum * bulk_momentum, 0.5);
            beta_bulk = pow(1 - 1 / (gamma_bulk * gamma_bulk), 0.5);
            doppler_factor_bulk = 1. / (gamma_bulk * (1. - beta_bulk * std::cos(theta_obs * karcst::pi / 180.)));
        }

        std::vector<double> get_timescale_electron_cyclosyn(std::vector<double> momentum);
        std::vector<double> get_timescale_electron_adiabatic(std::vector<double> momentum);
        std::vector<double> get_timescale_electron_compton_thomson(std::vector<double> momentum);
        std::vector<double> get_timescale_electron_compton(std::vector<double> momentum);
        std::vector<double> get_timescale_electron_acceleration(std::vector<double> momentum);

        std::vector<double> cyclosyn_selfabsorption_rate, cyclosyn_energy;
        std::vector<double> get_timescale_photon_cyclosyn_selfabsorption(std::vector<double> momentum);
        // read out the array, left for debugging
        // std::vector<double> get_timescale_photon_cyclosyn_selfabsorption_array(){return cyclosyn_selfabsorption_rate;};
        // std::vector<double> get_timescale_photon_cyclosyn_selfabsorption_energy(){return cyclosyn_energy;};
        std::vector<double> get_timescale_photon_escape(std::vector<double> momentum);


        double get_electron_max_momentum();
        double get_electron_break_momentum();

        std::vector<double> get_electron_momentum_grid();
        std::vector<double> get_electron_gamma_grid();
        std::vector<double> get_electron_momentum_number_density();
        std::vector<double> get_electron_gamma_number_density();

        // std::vector<double> additional_target_field_energy, additional_target_field_energy_density;
        std::vector<double> get_total_photon_target_energy();
        std::vector<double> get_total_photon_target_energy_density();
        std::vector<double> get_photon_target_energy(std::string name);
        std::vector<double> get_photon_target_energy_density(std::string name);
        // std::vector<double> get_photon_target_energy_density_black_body(std::string name);
        // double get_target_black_body_temperature(std::string name);
        // double get_target_black_body_energy_density(std::string name);
        // void set_target_black_body_temperature(std::string name, double new_temperature);
        // void set_target_black_body_energy_density(std::string name, double new_energy_density);

        std::vector<double>
            photon_energy_grid_electron_cyclosyn, photon_observed_luminosity_electron_cyclosyn,
            photon_energy_grid_electron_compton, photon_observed_luminosity_electron_compton,
            photon_energy_grid_total, photon_observed_luminosity_total;

        std::vector<double> get_observed_photon_energy_grid_electron_cyclosyn();
        std::vector<double> get_observed_photon_luminosity_electron_cyclosyn();
        std::vector<double> get_observed_photon_flux_electron_cyclosyn();
        std::vector<double> get_observed_photon_energy_grid_electron_compton();
        std::vector<double> get_observed_photon_luminosity_electron_compton();
        std::vector<double> get_observed_photon_flux_electron_compton();
        std::vector<double> get_observed_photon_energy_grid_total();
        std::vector<double> get_observed_photon_luminosity_total();
        std::vector<double> get_observed_photon_flux_total();

        std::vector<double> get_computation_times();

        void sum_jet_and_counterjet(size_t size, const std::vector<double> &input_en,
                                    const std::vector<double> &input_lum, std::vector<double> &en,
                                    std::vector<double> &lum);
        void sum_jet_only(size_t size, const std::vector<double> &input_en,
                          const std::vector<double> &input_lum, std::vector<double> &en,
                          std::vector<double> &lum);

        bool compton_calculation_necessary();
        /**
         * @brief function to compute ker_ic[][].
         *
         * cf. Jones, "Inverse Compton Scattering of Cosmic-Ray Electrons", P.R. 1965, Eqn. 13-15
         *
         * @param gamma Energy of incoming particle
         * @param alpha Energy of incoming photon / me c^2
         */
        double Fic(double gamma, double alpha);
        // utility functions to construct Fic
        double Fic_1(double z);          //< utility functions to construct Fic
        double Fic_2(double z);          //< utility functions to construct Fic
        double polylog(double z);        // utility functions to construct Fic, valid for -inf<z<1.0
        double polylog_smallz(double z); //< utility functions to calc polylog
        double polylog_rk4(double z);    //< utility functions to construct polylog
        double polylog_largez(double z); //< utility functions to construct polylog
    };

} // namespace bhjet
