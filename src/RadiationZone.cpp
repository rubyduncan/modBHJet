#include "RadiationZone.hpp"
#include "utils.hpp"
#include <iostream>
#include <chrono>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include "kariba/Mixed.hpp"
#include "kariba/Thermal.hpp"
#include "kariba/Powerlaw.hpp"
#include "kariba/Bknpower.hpp"
#include "kariba/Cyclosyn.hpp"
#include "kariba/Compton.hpp"
#include "kariba/constants.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{

    RadiationZone::~RadiationZone()
    {
        if (verbosity_level > 2)
            std::cout << "cleaning up spline memory" << std::endl;
        if (spline_electrons)
            gsl_spline_free(spline_electrons);
        if (spline_electrons_accel)
            gsl_interp_accel_free(spline_electrons_accel);
        if (spline_electrons_derivative)
            gsl_spline_free(spline_electrons_derivative);
        if (spline_electrons_derivative_accel)
            gsl_interp_accel_free(spline_electrons_derivative_accel);
    }

    void RadiationZone::add_target_photon_field(
        std::vector<double> energies, std::vector<double> energy_densities, 
        std::string name)
    {
        if (verbosity_level > 2)
            std::cout << "adding target: " << name << std::endl;
        target_component_names.emplace_back(name);
        target_component_radiation_energy_density.emplace_back(energy_densities);
        target_component_radiation_energy_grid.emplace_back(energies);
    }

    void RadiationZone::remove_target_photon_field(std::string name)
    {
        if (verbosity_level > 1)
            std::cout << "RadiationZone: removing target photon field: " << name << std::endl;

        for (size_t i=0;  i<target_component_names.size(); i++){
            if(target_component_names[i] == name) {
                target_component_names.erase(target_component_names.begin() + i);
                target_component_radiation_energy_density.erase(target_component_radiation_energy_density.begin() + i);
                target_component_radiation_energy_grid.erase(target_component_radiation_energy_grid.begin() + i);
                if (verbosity_level > 1)
                    std::cout << "RadiationZone: removed target photon field "<< i << ": " << name << std::endl;
            }
        }
    }

    void RadiationZone::compute_particles(bool reset_photon_targets)
    {
        // reset computation times
        for (size_t i = 0; i < 3; i++)
            computation_times[i] = 0.0;
        auto tstart = std::chrono::steady_clock::now();

        // clean up existing gsl interpolation
        if (verbosity_level > 2)
            std::cout << "cleaning up spline memory" << std::endl;
        if (spline_electrons)
            gsl_spline_free(spline_electrons);
        if (spline_electrons_accel)
            gsl_interp_accel_free(spline_electrons_accel);
        if (spline_electrons_derivative)
            gsl_spline_free(spline_electrons_derivative);
        if (spline_electrons_derivative_accel)
            gsl_interp_accel_free(spline_electrons_derivative_accel);

        if (verbosity_level > 2)
            std::cout << "allocating spline memory" << std::endl;
        spline_electrons_accel = gsl_interp_accel_alloc();
        spline_electrons = gsl_spline_alloc(gsl_interp_steffen, n_bins_e);
        spline_electrons_derivative_accel = gsl_interp_accel_alloc();
        spline_electrons_derivative = gsl_spline_alloc(gsl_interp_steffen, n_bins_e);

        // calculate Urad from targets, could also extend this to KN cooling
        radiation_energy_density = 0;
        for (size_t i = 0; i < target_component_names.size(); i++)
        {

            double dlogE = std::log(
                target_component_radiation_energy_grid[i][1] / 
                target_component_radiation_energy_grid[i][0]);

            for (size_t e = 0; e < target_component_radiation_energy_grid[i].size(); e++)
            {
                radiation_energy_density += target_component_radiation_energy_density[i][e] * dlogE;
            }
            // use here a doppler factor ^2 instead of gamma^2 to account for anisotropy (roughly), see Dermer 1995
            // radiation_energy_density += std::pow(doppler_factor_bulk, 2) * target_vector_blackbody[i].energy_density;
        }

        if (fraction_nonthermal_electrons == 0.)
        {
            if (verbosity_level > 1)
                std::cout << "initializing a thermal distribution" << std::endl;
            // Thermal spectrum only
            electrons_thermal = kariba::Thermal(n_bins_e);
            electrons_thermal.set_temp_kev(electron_temperature);
            electrons_thermal.set_p();
            electrons_thermal.set_norm(electron_number_density);
            electrons_thermal.set_ndens();
            if (verbosity_level > 2)
                electrons_thermal.test();

            gsl_spline_init(spline_electrons, electrons_thermal.get_gamma().data(),
                            electrons_thermal.get_gdens().data(), n_bins_e);
            gsl_spline_init(spline_electrons_derivative,
                            electrons_thermal.get_gamma().data(),
                            electrons_thermal.get_pdensp2_diff_logp().data(), n_bins_e);
        }
        else if (fraction_nonthermal_electrons < 0.5)
        {
            if (verbosity_level > 1)
                std::cout << "initializing a mixed distribution (thermal + powerlaw)" << std::endl;
            // mixed thermal + non-thermal
            electrons_mixed = kariba::Mixed(n_bins_e);
            electrons_mixed.set_cutoff_type(cutoff_type);
            electrons_mixed.set_temp_kev(electron_temperature);
            electrons_mixed.set_pspec(index_injected_electrons);
            electrons_mixed.set_plfrac(fraction_nonthermal_electrons);
            electrons_mixed.set_p(radiation_energy_density, magnetic_field, factor_break_electrons, radius, factor_max_energy_electrons);
            electrons_mixed.set_norm(electron_number_density);
            electrons_mixed.set_ndens();
            if (verbosity_level > 2)
                electrons_mixed.test();
            electrons_mixed.cooling_steadystate(radiation_energy_density, electron_number_density, magnetic_field, radius, factor_break_electrons);

            gsl_spline_init(spline_electrons, electrons_mixed.get_gamma().data(),
                            electrons_mixed.get_gdens().data(), n_bins_e);
            gsl_spline_init(spline_electrons_derivative,
                            electrons_mixed.get_gamma().data(),
                            electrons_mixed.get_pdensp2_diff_logp().data(), n_bins_e);
        }
        else if (fraction_nonthermal_electrons < 1.)
        {
            // broken powerlaw approximation
            if (verbosity_level > 1)
                std::cout << "initializing a broken powerlaw distribution" << std::endl;
            // determine the momentum of the break from temperature
            electrons_thermal = kariba::Thermal(n_bins_e);
            electrons_thermal.set_temp_kev(electron_temperature);
            electrons_thermal.set_p();
            electrons_thermal.set_norm(electron_number_density);
            electrons_thermal.set_ndens();
            double p_break = electrons_thermal.av_p();

            electrons_bpl = kariba::Bknpower(n_bins_e);
            electrons_bpl.set_cutoff_type(cutoff_type);
            electrons_bpl.set_pspec1(-2.);
            electrons_bpl.set_pspec2(index_injected_electrons);
            electrons_bpl.set_p(0.1 * p_break, p_break, radiation_energy_density, magnetic_field, factor_break_electrons, radius, factor_max_energy_electrons);
            electrons_bpl.set_norm(electron_number_density);
            electrons_bpl.set_ndens();
            electrons_bpl.cooling_steadystate(radiation_energy_density, electron_number_density, magnetic_field, radius, factor_break_electrons);
            if (verbosity_level > 2)
                electrons_bpl.test();

            gsl_spline_init(spline_electrons, electrons_bpl.get_gamma().data(),
                            electrons_bpl.get_gdens().data(), n_bins_e);
            gsl_spline_init(spline_electrons_derivative,
                            electrons_bpl.get_gamma().data(),
                            electrons_bpl.get_pdensp2_diff_logp().data(), n_bins_e);
        }
        else if (fraction_nonthermal_electrons == 1.)
        {
            if (verbosity_level > 1)
                std::cout << "initializing a powerlaw distribution" << std::endl;
            // Powerlaw spectrum only
            // determine the min momentum from temperature
            electrons_thermal = kariba::Thermal(n_bins_e);
            electrons_thermal.set_temp_kev(electron_temperature);
            electrons_thermal.set_p();
            electrons_thermal.set_norm(electron_number_density);
            electrons_thermal.set_ndens();
            double p_min = electrons_thermal.av_p();

            electrons_pl = kariba::Powerlaw(n_bins_e);
            electrons_pl.set_cutoff_type(cutoff_type);
            electrons_pl.set_pspec(index_injected_electrons);
            electrons_pl.set_p(p_min, radiation_energy_density, magnetic_field, factor_break_electrons, radius, factor_max_energy_electrons);
            electrons_pl.set_norm(electron_number_density);
            electrons_pl.set_ndens();
            electrons_pl.cooling_steadystate(radiation_energy_density, electron_number_density, magnetic_field, radius, factor_break_electrons);
            if (verbosity_level > 2)
                electrons_pl.test();

            gsl_spline_init(spline_electrons, electrons_pl.get_gamma().data(),
                            electrons_pl.get_gdens().data(), n_bins_e);
            gsl_spline_init(spline_electrons_derivative,
                            electrons_pl.get_gamma().data(),
                            electrons_pl.get_pdensp2_diff_logp().data(), n_bins_e);
        }
        else
        {
            throw std::out_of_range("fraction_nonthermal_electrons has to be <=1!");
        }
        // same for protons...

        if (profile_time)
            computation_times[0] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tstart).count();
    }

    void RadiationZone::compute_radiation()
    {
        if (verbosity_level > 1)
            std::cout << "Computing radiation using default energy array" << std::endl;
        int N = 210;
        double lo = -19; // 1e-19 erg
        double hi = 2;   // 1e2 erg

        std::vector<double> v(N);
        double d = (hi - lo) / (N - 1);

        for (int i = 0; i < N; i++)
            v[i] = std::pow(10.0, lo + d * i);
        compute_radiation(v);
    }

    void RadiationZone::compute_radiation(std::vector<double> obs_energy_grid)
    {
        if (verbosity_level > 1)
            std::cout << "Computing radiation using given energy array" << std::endl;

        double gmin = get_electron_gamma_grid()[0];
        double gmax = get_electron_gamma_grid()[n_bins_e - 1];

        double syn_min, syn_max;
        // calculate emission of each zone
        // note: the syn_en array is used for the seed photon fields in the IC
        // part, so it needs to include also the black body and disk part. This
        // is why the maximum frequency is taken as the maximum of the two scale
        // frequencies.
        syn_min = 0.1 * std::pow(gmin, 2.) * karcst::charg * magnetic_field /
                  (2. * karcst::pi * karcst::emgm * karcst::cee);
        // syn_min = std::max(syn_min, 1.6e-18 / karcst::herg); // don't go below 1e-8eV

        syn_max = 50. * std::pow(gmax, 2.) * karcst::charg * magnetic_field /
                  (2. * karcst::pi * karcst::emgm * karcst::cee);
        // compare if external targets extend to higher energies
        if (target_component_names.size() > 0)
        {
            for (size_t i = 0; i < target_component_names.size(); i++)
            {
                // this simply assumes that the energy grids are not extensively larger than 
                // the spectra's relevant flux ranges
                size_t Nbins_i = target_component_radiation_energy_grid[i].size();
                if (target_component_radiation_energy_grid[i][Nbins_i-1] > syn_max * karcst::herg) {
                    syn_max = target_component_radiation_energy_grid[i][Nbins_i-1] / karcst::herg;
                }
                if(verbosity_level > 2) {
                    std::cout << "syn_max:" << syn_max << "Hz, target_component_radiation_energy_grid[i][Nbins_i-1]:" << 
                    target_component_radiation_energy_grid[i][Nbins_i-1] << "erg = " << 
                    target_component_radiation_energy_grid[i][Nbins_i-1]/ karcst::herg << "Hz " << std::endl;
                }
            }
        }

        size_t nsyn = (size_t)(std::log10(syn_max) - std::log10(syn_min)) * syn_res;
        photon_energy_grid_electron_cyclosyn = std::vector<double>(nsyn, 1e-100);
        photon_observed_luminosity_electron_cyclosyn = std::vector<double>(nsyn, 1e-100);

        auto tstart = std::chrono::steady_clock::now();
        kariba::Cyclosyn Syncro(nsyn);
        if(verbosity_level > 2) std::cout << "syn_min:" << syn_min << "Hz, syn_max:" << syn_max << "Hz, nsyn:" << nsyn << std::endl;
        // std::cout << "after" << std::endl;
        Syncro.set_frequency(syn_min, syn_max);
        // std::cout << "after freq" << std::endl;

        // calculate cyclosynchrotron spectrum
        // Set up the calculation by reading in magnetic
        // field,beaming,volume,counterjet presence
        Syncro.set_bfield(magnetic_field);
        Syncro.set_beaming(std::max(1e-10, theta_obs), beta_bulk, doppler_factor_bulk);
        Syncro.set_geometry(geometry, radius, height);
        Syncro.set_counterjet(include_counterjet);
        // std::cout << "before" << std::endl;

        // read out syn self-abs timescale
        Syncro.cycsyn_spectrum(gmin, gmax, spline_electrons, spline_electrons_accel, spline_electrons_derivative, spline_electrons_derivative_accel);
        std::vector<double> syn_abs_rate = Syncro.get_cyclosyn_absorption_rate();
        cyclosyn_selfabsorption_rate.resize(syn_abs_rate.size());
        for (size_t k = 0; k < syn_abs_rate.size(); k++) {
            cyclosyn_selfabsorption_rate[k] = 1/ syn_abs_rate[k];
        }
        cyclosyn_energy = Syncro.get_energy();

        if (include_counterjet)
        {
            sum_jet_and_counterjet(nsyn, Syncro.get_energy_obs(), Syncro.get_nphot_obs(),
                                   photon_energy_grid_electron_cyclosyn, photon_observed_luminosity_electron_cyclosyn);
        }
        else
        {
            sum_jet_only(nsyn, Syncro.get_energy_obs(), Syncro.get_nphot_obs(),
                         photon_energy_grid_electron_cyclosyn, photon_observed_luminosity_electron_cyclosyn);
        }
        if (profile_time)
            computation_times[1] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tstart).count();
        if (verbosity_level > 2)
            Syncro.test();

        photon_energy_grid_total = obs_energy_grid;
        photon_observed_luminosity_total = std::vector<double>(obs_energy_grid.size(), 1e-100);
        add_emission_on_interpolated_grid(
            photon_energy_grid_electron_cyclosyn, photon_observed_luminosity_electron_cyclosyn,
            photon_energy_grid_total, photon_observed_luminosity_total);

        std::vector<double> syn_energy_density(nsyn, 1e-100);
        for (size_t i = 0; i < nsyn; i++){
            syn_energy_density[i] = pow(Syncro.get_energy()[i], 2.) * Syncro.get_nphot()[i]/
                       (karcst::cee * karcst::herg * Syncro.get_energy()[i] * karcst::pi * radius * radius);
        }
        remove_target_photon_field("cyclosynchrotron");
        add_target_photon_field(Syncro.get_energy(), syn_energy_density, "cyclosynchrotron");

        // scattered photon energy grid
        // find max of absorbed synch spectrum
        double temp_lum = 0.;
        size_t ind_max = 0;
        for (size_t i = 0; i < nsyn; i++)
        {
            if (photon_observed_luminosity_electron_cyclosyn[i] > temp_lum)
            {
                temp_lum = photon_observed_luminosity_electron_cyclosyn[i];
                ind_max = i;
            }
        }
        // set min frequency somewhere below
        double com_min = 0.1 * photon_energy_grid_electron_cyclosyn[ind_max] / karcst::herg;

        // this can be done more clever: com_max = min( 10 * min(gamma^2 * Etarget/h, gamma*me*c^2/h), grid_max),
        // Etarget = max(syn_max, other target energy eg BB temperature * 10)
        // grid_max = cut off, e.g. if I only have data until keV, no need to calc higher
        double E_el_max = gmax * karcst::emerg;
        double Thomson_max = gmax * gmax * syn_max;
        double com_max = 100 * std::min(E_el_max / karcst::herg, Thomson_max);
        com_max = std::min(com_max, obs_energy_grid[obs_energy_grid.size() - 1] / karcst::herg);
        size_t ncom = (size_t)(std::log10(com_max) - std::log10(com_min)) * com_res;

        photon_energy_grid_electron_compton = std::vector<double>(ncom, 1e-100);
        photon_observed_luminosity_electron_compton = std::vector<double>(ncom, 1e-100);
        if (profile_time)
            tstart = std::chrono::steady_clock::now();
        kariba::Compton InvCompton(ncom, nsyn);
        InvCompton.set_frequency(com_min, com_max);
        for (size_t i = 0; i < ncom; i++)
        { // needed to make sure the energy grid is initialised even when skipping Compton calulcation
            photon_energy_grid_electron_compton[i] = InvCompton.get_energy()[i];
        }

        // std::cout << "com min " << com_min << "Hz, com max " << com_max << "Hz " << std::endl;
        // std::cout << "doppler_factor_bulk " << doppler_factor_bulk << std::endl;

        // readout of the target fields needed to compute Compton cooling time
        // if (compton_calculation_necessary()) {

        // Set up the calculation by reading in/calculating
        // beaming,volume,counterjet presence,tau
        // std::cout << "set beaming" << std::endl;
        InvCompton.set_beaming(std::max(1e-10, theta_obs), beta_bulk, doppler_factor_bulk);
        InvCompton.set_geometry(geometry, radius, height);
        InvCompton.set_counterjet(include_counterjet);
        // std::cout << "set tau" << std::endl;
        InvCompton.set_tau(electron_number_density, electron_temperature);
        // std::cout << "set niter?" << std::endl;
        // Multiple scatters only if ypar and tau are large enough
        if (InvCompton.get_ypar() > 1.e-2 && InvCompton.get_tau() > 5.e-2)
        {
            // std::cout << "set niter(15)" << std::endl;
            InvCompton.set_niter(15);
        }

        // std::cout << "set syn seed" << std::endl;
        // Cyclosynchrotron photons are always considered in the scattering
        // InvCompton.cyclosyn_seed(Syncro.get_energy(), Syncro.get_nphot());
        InvCompton.set_target_energy_array(Syncro.get_energy());

        // std::cout << "set bb seed" << std::endl;
        for (size_t i = 0; i < target_component_names.size(); i++)
        {
            InvCompton.add_target_energy_density(
                target_component_radiation_energy_grid[i], 
                target_component_radiation_energy_density[i]);
            // InvCompton.bb_seed_kev(std::pow(doppler_factor_bulk, 2) * target_vector_blackbody[i].energy_density,
            //                        doppler_factor_bulk * target_vector_blackbody[i].temperature);
            // use here a doppler factor ^2 to roughly include anisotropy
        }

        // std::cout << "before Compton switch" << std::endl;
        // if(additional_target_field_energy_density.size()>0){
        //     std::vector<double> target_extra(nsyn, 1e-100);
        //     std::vector<double> syn_energy = Syncro.get_energy();
        //     add_emission_on_interpolated_grid(
        //         additional_target_field_energy, additional_target_field_energy_density,
        //         syn_energy, target_extra);
        //     InvCompton.add_seed(Syncro.get_energy(), target_extra);
        //     // is here a doppler factor ^2 or gamma^2 needed?
        // }
        total_target_radiation_energy_grid = InvCompton.get_target_energy();
        total_target_radiation_energy_density = std::vector<double>(total_target_radiation_energy_grid.size());
        for (size_t i = 0; i < total_target_radiation_energy_grid.size(); i++){
            // convert to energy density
            total_target_radiation_energy_density[i] = InvCompton.get_target_diff_spec()[i] * std::pow(total_target_radiation_energy_grid[i], 2);
        }

        if (compton_calculation_necessary()) //legacy: compton switch
        {
            // this is the comp. expensive call
            InvCompton.compton_spectrum(gmin, gmax, spline_electrons, spline_electrons_accel);

            if (verbosity_level>2) std::cout << "compton sum" << std::endl;
            if (include_counterjet)
            {
                sum_jet_and_counterjet(ncom, InvCompton.get_energy_obs(), InvCompton.get_nphot_obs(),
                                       photon_energy_grid_electron_compton, photon_observed_luminosity_electron_compton);
            }
            else
            {
                sum_jet_only(ncom, InvCompton.get_energy_obs(), InvCompton.get_nphot_obs(),
                             photon_energy_grid_electron_compton, photon_observed_luminosity_electron_compton);
            }

            // std::cout << "add com " << photon_energy_grid_electron_compton.size() << "," <<
            //      photon_observed_luminosity_electron_compton.size() << "," <<
            //      photon_energy_grid_total.size() << "," <<
            //      photon_observed_luminosity_total.size() << "," << std::endl;
            add_emission_on_interpolated_grid(
                photon_energy_grid_electron_compton, photon_observed_luminosity_electron_compton,
                photon_energy_grid_total, photon_observed_luminosity_total);
        }
        if (profile_time)
            computation_times[2] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tstart).count();
    }

    std::vector<double> RadiationZone::get_electron_momentum_grid()
    {
        const std::vector<double> &vec =
            (fraction_nonthermal_electrons == 0) ? electrons_thermal.get_p() : (fraction_nonthermal_electrons < 0.5) ? electrons_mixed.get_p()
                                                                           : (fraction_nonthermal_electrons < 1.)    ? electrons_bpl.get_p()
                                                                           : (fraction_nonthermal_electrons == 1)    ? electrons_pl.get_p()
                                                                                                                     : throw std::out_of_range("fraction_nonthermal_electrons has to be <=1!");
        return vec;
    }
    std::vector<double> RadiationZone::get_electron_gamma_grid()
    {
        const std::vector<double> &vec =
            (fraction_nonthermal_electrons == 0) ? electrons_thermal.get_gamma() : (fraction_nonthermal_electrons < 0.5) ? electrons_mixed.get_gamma()
                                                                               : (fraction_nonthermal_electrons < 1.)    ? electrons_bpl.get_gamma()
                                                                               : (fraction_nonthermal_electrons == 1)    ? electrons_pl.get_gamma()
                                                                                                                         : throw std::out_of_range("fraction_nonthermal_electrons has to be <=1!");
        return vec;
    }
    std::vector<double> RadiationZone::get_electron_momentum_number_density()
    {
        std::vector<double> mom = get_electron_momentum_grid();
        std::vector<double> vec =
            (fraction_nonthermal_electrons == 0) ? electrons_thermal.get_pdens() : (fraction_nonthermal_electrons < 0.5) ? electrons_mixed.get_pdens()
                                                                               : (fraction_nonthermal_electrons < 1.)    ? electrons_bpl.get_pdens()
                                                                               : (fraction_nonthermal_electrons == 1)    ? electrons_pl.get_pdens()
                                                                                                                         : throw std::out_of_range("fraction_nonthermal_electrons has to be <=1!");
        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] *= mom[i];
        }
        return vec;
    }
    std::vector<double> RadiationZone::get_electron_gamma_number_density()
    {
        std::vector<double> gamma = get_electron_gamma_grid();
        std::vector<double> vec =
            (fraction_nonthermal_electrons == 0) ? electrons_thermal.get_gdens() : (fraction_nonthermal_electrons < 0.5) ? electrons_mixed.get_gdens()
                                                                               : (fraction_nonthermal_electrons < 1.)    ? electrons_bpl.get_gdens()
                                                                               : (fraction_nonthermal_electrons == 1)    ? electrons_pl.get_gdens()
                                                                                                                         : throw std::out_of_range("fraction_nonthermal_electrons has to be <=1!");
        for (size_t i = 0; i < vec.size(); i++)
        {
            vec[i] *= gamma[i];
        }
        return vec;
    }

    std::vector<double> RadiationZone::get_total_photon_target_energy()
    {
        return total_target_radiation_energy_grid;
    }
    std::vector<double> RadiationZone::get_total_photon_target_energy_density()
    {
        return total_target_radiation_energy_density;
    }


    std::vector<double> RadiationZone::get_photon_target_energy_density(std::string name)
    {
        bool target_in_list = false;
        for (size_t i = 0; i < target_component_names.size(); i++)
        {
            if (name == target_component_names[i])
            {
                target_in_list = true;
                return target_component_radiation_energy_density[i];
            }
        }
        if (!target_in_list)
        {
            throw std::out_of_range("photon target field called " + name + " not found!");
        }
        return std::vector<double>();
    }
    std::vector<double> RadiationZone::get_photon_target_energy(std::string name)
    {
        bool target_in_list = false;
        for (size_t i = 0; i < target_component_names.size(); i++)
        {
            if (name == target_component_names[i])
            {
                target_in_list = true;
                return target_component_radiation_energy_grid[i];
            }
        }
        if (!target_in_list)
        {
            throw std::out_of_range("photon target field called " + name + " not found!");
        }
        return std::vector<double>();
    }


    std::vector<double> RadiationZone::get_observed_photon_energy_grid_electron_cyclosyn()
    {
        return photon_energy_grid_electron_cyclosyn;
    }
    std::vector<double> RadiationZone::get_observed_photon_luminosity_electron_cyclosyn()
    {
        return photon_observed_luminosity_electron_cyclosyn;
    }
    std::vector<double> RadiationZone::get_observed_photon_flux_electron_cyclosyn()
    {
        std::vector<double> flux(photon_observed_luminosity_electron_cyclosyn.size(), 1e-100);
        for (size_t i = 0; i < flux.size(); i++)
        {
            flux[i] = kariba_luminosity_to_number_flux(photon_observed_luminosity_electron_cyclosyn[i], redshift, distance);
        }
        return flux;
    }

    std::vector<double> RadiationZone::get_observed_photon_energy_grid_electron_compton()
    {
        return photon_energy_grid_electron_compton;
    }
    std::vector<double> RadiationZone::get_observed_photon_luminosity_electron_compton()
    {
        return photon_observed_luminosity_electron_compton;
    }
    std::vector<double> RadiationZone::get_observed_photon_flux_electron_compton()
    {
        std::vector<double> flux(photon_observed_luminosity_electron_compton.size(), 1e-100);
        for (size_t i = 0; i < flux.size(); i++)
        {
            flux[i] = kariba_luminosity_to_number_flux(photon_observed_luminosity_electron_compton[i], redshift, distance);
        }
        return flux;
    }

    std::vector<double> RadiationZone::get_observed_photon_energy_grid_total()
    {
        return photon_energy_grid_total;
    }
    std::vector<double> RadiationZone::get_observed_photon_luminosity_total()
    {
        return photon_observed_luminosity_total;
    }
    std::vector<double> RadiationZone::get_observed_photon_flux_total()
    {
        std::vector<double> flux(photon_observed_luminosity_total.size(), 1e-100);
        for (size_t i = 0; i < flux.size(); i++)
        {
            flux[i] = kariba_luminosity_to_number_flux(photon_observed_luminosity_total[i], redshift, distance);
        }
        return flux;
    }

    // This function takes the observed arrays of the Cyclosyn and Compton classes
    // for jet/counterjet sums up the contributions of both and stores them in one
    // array of observed frequencies and one of comoving luminosities
    void RadiationZone::sum_jet_and_counterjet(size_t size, const std::vector<double> &input_en,
                                               const std::vector<double> &input_lum, std::vector<double> &en,
                                               std::vector<double> &lum)
    {
        double en_cj_min, en_j_min, en_cj_max, en_j_max, einc;
        std::vector<double> en_j(size, 1e-100);
        std::vector<double> lum_j(size, 1e-100);
        en_j_min = input_en[0];
        en_j_max = input_en[size - 1];

        std::vector<double> en_cj(size, 1e-100);
        std::vector<double> lum_cj(size, 1e-100);
        en_cj_min = input_en[size];
        en_cj_max = input_en[2 * size - 1];

        einc = (std::log10(en_j_max) - std::log10(en_cj_min)) / static_cast<double>(size - 1);

        for (size_t i = 0; i < size; i++)
        {
            en[i] = std::pow(10., std::log10(en_cj_min) + static_cast<double>(i) * einc);
            en_j[i] = input_en[i];
            en_cj[i] = input_en[i + size];
            lum_j[i] = std::max(input_lum[i], 1.e-50);
            lum_cj[i] = std::max(input_lum[i + size], 1.e-50);
        }

        gsl_interp_accel *acc_j = gsl_interp_accel_alloc();
        gsl_spline *spline_j = gsl_spline_alloc(gsl_interp_akima, size);
        gsl_spline_init(spline_j, en_j.data(), lum_j.data(), size);

        gsl_interp_accel *acc_cj = gsl_interp_accel_alloc();
        gsl_spline *spline_cj = gsl_spline_alloc(gsl_interp_akima, size);
        gsl_spline_init(spline_cj, en_cj.data(), lum_cj.data(), size);

        for (size_t i = 0; i < size; i++)
        {
            if (i == 0)
            {
                lum[i] = lum_cj[i];
            }
            else if (i == size - 1)
            {
                lum[i] = lum_j[i];
            }
            else if (en[i] < en_j_min)
            {
                lum[i] = gsl_spline_eval(spline_cj, en[i] * 1.0000001, acc_cj);
            }
            else if (en[i] < en_cj_max)
            {
                lum[i] =
                    gsl_spline_eval(spline_j, en[i], acc_j) + gsl_spline_eval(spline_cj, en[i], acc_cj);
            }
            else
            {
                lum[i] = gsl_spline_eval(spline_j, en[i] * 0.999999,
                                         acc_j); // note: the factor 0.999 is to avoid
                                                 // occasional gsl interpolation errors
                                                 // due to numerical inaccuracies
            }
        }

        gsl_spline_free(spline_j), gsl_interp_accel_free(acc_j);
        gsl_spline_free(spline_cj), gsl_interp_accel_free(acc_cj);

        return;
    }

    // This function takes the observed arrays of the Cyclosyn and Compton classes
    // for jet/counterjet sums up the contributions of both and stores them in one
    // array of observed frequencies and one of comoving luminosities
    void RadiationZone::sum_jet_only(size_t size, const std::vector<double> &input_en,
                                     const std::vector<double> &input_lum, std::vector<double> &en,
                                     std::vector<double> &lum)
    {
        double en_cj_min, en_j_min, en_cj_max, en_j_max, einc;
        std::vector<double> en_j(size, 1e-100);
        std::vector<double> lum_j(size, 1e-100);
        en_j_min = input_en[0];
        en_j_max = input_en[size - 1];

        einc = (std::log10(en_j_max) - std::log10(en_j_min)) / static_cast<double>(size - 1);
        for (size_t i = 0; i < size; i++)
        {
            en[i] = std::pow(10., std::log10(en_j_min) + static_cast<double>(i) * einc);
            en_j[i] = input_en[i];
            lum_j[i] = std::max(input_lum[i], 1.e-50);
        }

        gsl_interp_accel *acc_j = gsl_interp_accel_alloc();
        gsl_spline *spline_j = gsl_spline_alloc(gsl_interp_akima, size);
        gsl_spline_init(spline_j, en_j.data(), lum_j.data(), size);

        for (size_t i = 0; i < size; i++)
        {
            if (i == 0)
            {
                lum[i] = lum_j[i];
            }
            else if (i == size - 1)
            {
                lum[i] = lum_j[i];
            }
            else
            {
                lum[i] = gsl_spline_eval(spline_j, en[i] * 0.999999, acc_j);
                // note: the factor 0.999 is to avoid
                // occasional gsl interpolation errors
                // due to numerical inaccuracies
            }
        }

        gsl_spline_free(spline_j), gsl_interp_accel_free(acc_j);

        return;
    }

    std::vector<double> RadiationZone::get_timescale_electron_cyclosyn(std::vector<double> momentum)
    {
        double gamma = 1;
        double U_B = std::pow(magnetic_field, 2.) / (8. * karcst::pi);
        double pdot_B = (4. * karcst::sigtom * karcst::cee * U_B) /
                        (3. * karcst::emgm * std::pow(karcst::cee, 2.));

        std::vector<double> t_syn(momentum.size(), 1e100);
        for (size_t i = 0; i < momentum.size(); i++)
        {
            gamma = std::pow(1 + std::pow(momentum[i] / karcst::emgm / karcst::cee, 2), 0.5);
            t_syn[i] = 1 / pdot_B / gamma;
        }

        return t_syn;
    }

    std::vector<double> RadiationZone::get_timescale_electron_compton_thomson(std::vector<double> momentum)
    {
        double gamma = 1;
        double U_rad = 1e-100;

        if (total_target_radiation_energy_grid.size() < 2 ||
            total_target_radiation_energy_density.size() !=
                total_target_radiation_energy_grid.size())
        {
            std::cout << "Invalid target radiation field in RadiationZone, egridsize=" << total_target_radiation_energy_grid.size() << ", densitygridsize=" << total_target_radiation_energy_density.size() << std::endl;
            throw std::runtime_error(
                "Invalid target radiation field in RadiationZone (Compton not calculated?)");
        }

        double dlogE = std::log(total_target_radiation_energy_grid[1] / total_target_radiation_energy_grid[0]);
        for (size_t i = 0; i < total_target_radiation_energy_density.size(); i++)
        {
            U_rad += total_target_radiation_energy_density[i] * dlogE;
        }

        double pdot_rad = (4. * karcst::sigtom * karcst::cee * U_rad) /
                          (3. * karcst::emgm * std::pow(karcst::cee, 2.));

        std::vector<double> t_com(momentum.size(), 1e100);
        for (size_t i = 0; i < momentum.size(); i++)
        {
            gamma = std::pow(1 + std::pow(momentum[i] / karcst::emgm / karcst::cee, 2), 0.5);
            t_com[i] = 1 / pdot_rad / gamma;
        }

        return t_com;
    }

    std::vector<double> RadiationZone::get_timescale_electron_compton(std::vector<double> momentum)
    {
        double fac_gamma2p, gamma, alpha, integral, betagamma, dn_dalpha;

        if (total_target_radiation_energy_grid.size() < 2 ||
            total_target_radiation_energy_density.size() !=
                total_target_radiation_energy_grid.size())
        {
            std::cout << "Invalid target radiation field in RadiationZone, egridsize=" << total_target_radiation_energy_grid.size() << ", densitygridsize=" << total_target_radiation_energy_density.size() << std::endl;
            throw std::runtime_error(
                "Invalid target radiation field in RadiationZone (Compton not calculated?)");
        }

        double dlogEph = std::log(total_target_radiation_energy_grid[1] / total_target_radiation_energy_grid[0]);
        double fac_comp = 0.5 * karcst::pi * std::pow(karcst::re0, 2) * karcst::cee;
        std::vector<double> t_com(momentum.size(), 1e100);
        for (size_t i = 0; i < momentum.size(); i++)
        {
            // dp = gamma/p dgamma -> t = p/(dp/dt) = p^2/gamma / (dgamma/dt)
            gamma = std::pow(1 + std::pow(momentum[i] / karcst::emgm / karcst::cee, 2), 0.5);
            fac_gamma2p = std::pow(momentum[i] / karcst::emgm / karcst::cee, 2) / gamma;
            integral = 0;
            for (size_t j = 0; j < total_target_radiation_energy_grid.size(); j++)
            {
                alpha = total_target_radiation_energy_grid[j]/karcst::emerg;
                betagamma = momentum[i] / (karcst::emgm * karcst::cee);
                dn_dalpha = karcst::emerg * total_target_radiation_energy_density[j]* std::pow(total_target_radiation_energy_grid[j], -2);// need target density here...
                integral += dlogEph * fac_comp * Fic(gamma, alpha) * dn_dalpha / alpha / betagamma / gamma; 
            }

            t_com[i] = fac_gamma2p / integral;
        }

        return t_com;
    }

    std::vector<double> RadiationZone::get_timescale_electron_adiabatic(std::vector<double> momentum)
    {
        std::vector<double> t_adi(momentum.size(), 1e100);
        for (size_t i = 0; i < momentum.size(); i++)
        {
            t_adi[i] = 1 / (factor_break_electrons * karcst::cee / radius);
        }

        return t_adi;
    }

    std::vector<double> RadiationZone::get_timescale_electron_acceleration(std::vector<double> momentum)
    {
        std::vector<double> t_acc(momentum.size(), 1e100);
        for (size_t i = 0; i < momentum.size(); i++)
        {
            t_acc[i] = 4 * momentum[i] / (3 * factor_max_energy_electrons * karcst::charg * magnetic_field);
        }

        return t_acc;
    }

    std::vector<double> RadiationZone::get_timescale_photon_cyclosyn_selfabsorption(std::vector<double> momentum)
    {
        size_t n = cyclosyn_energy.size();

        std::vector<double> log_energy(n), log_rate(n);
        for (size_t i = 0; i < n; i++) {
            log_energy[i] = std::log10(cyclosyn_energy[i]);
            log_rate[i]   = std::log10(cyclosyn_selfabsorption_rate[i]);
        }

        gsl_interp_accel* acc = gsl_interp_accel_alloc();
        gsl_spline* spline = gsl_spline_alloc(gsl_interp_steffen, n);
        gsl_spline_init(spline, log_energy.data(), log_rate.data(), n);

        double pmin = cyclosyn_energy.front() / karcst::cee;
        double pmax = cyclosyn_energy.back() / karcst::cee;

        std::vector<double> result(momentum.size());
        for (size_t i = 0; i < momentum.size(); i++) {
            if (momentum[i] < pmin || momentum[i] > pmax) {
                result[i] = 1e100;
            } else {
                double log_en = std::log10(momentum[i] * karcst::cee);  // convert back to energy
                result[i] = std::pow(10., gsl_spline_eval(spline, log_en, acc));
            }
        }

        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);

        return result;
    }

    std::vector<double> RadiationZone::get_timescale_photon_escape(std::vector<double> momentum)
    {
        double t_escape = radius / karcst::cee;
        if(geometry == "cylinder") {
            t_escape /= 2;
        } else {
            t_escape /= 3;
        }
        return std::vector<double>(momentum.size(), t_escape);
    }

    std::vector<double> RadiationZone::get_computation_times()
    {
        return computation_times;
    }

    double RadiationZone::get_electron_max_momentum()
    {
        std::vector<double> momentum = get_electron_momentum_grid();
        std::vector<double> t_acc = get_timescale_electron_acceleration(momentum);
        std::vector<double> t_adi = get_timescale_electron_adiabatic(momentum);
        std::vector<double> t_syn = get_timescale_electron_cyclosyn(momentum);
        std::vector<double> t_com = get_timescale_electron_compton(momentum);

        size_t imax = momentum.size() - 1; // default: highest momentum

        for (size_t i = 0; i < momentum.size(); ++i)
        {
            double t_loss = 1.0 / (1.0 / t_adi[i] + 1.0 / t_syn[i] + 1.0 / t_com[i]);

            if (t_acc[i] > t_loss)
            {
                imax = (i == 0) ? 0 : i - 1;
                break;
            }
        }
        return momentum[imax];
    }
    double RadiationZone::get_electron_break_momentum()
    {
        std::vector<double> momentum = get_electron_momentum_grid();
        std::vector<double> t_adi = get_timescale_electron_adiabatic(momentum);
        std::vector<double> t_syn = get_timescale_electron_cyclosyn(momentum);
        std::vector<double> t_com = get_timescale_electron_compton(momentum);

        size_t ibreak = momentum.size() - 1;

        for (size_t i = 0; i < momentum.size(); ++i)
        {
            double t_rad = 1.0 / (1.0 / t_syn[i] + 1.0 / t_com[i]);

            if (t_adi[i] > t_rad)
            {
                ibreak = (i == 0) ? 0 : i - 1;
                break;
            }
        }
        return momentum[ibreak];
    }

    bool RadiationZone::compton_calculation_necessary()
    {
        if (force_compton_calculation)
            return compton_switch;
        else
        {
            double Lumnorm, Ub, Usyn, Lsyn, Lcom;
            // this is a rough estimate comparing the power in syn. to the power in IC:
            // P_syn/IC ~ V * n_e * sigma_T * c * < beta^2 gamma^2> * u_B/target
            // with u_target ~ P_syn * R/c / V + doppler^2 * u'_rad
            // and u_B = B^2/8pi
            // then P_com/P_syn ~ n_e * sigma_T * R * <gamma^2> + doppler^2 * u'_rad / u_B
            double average_gamma_squared =
                (fraction_nonthermal_electrons == 0) ? electrons_thermal.av_gammasq() : (fraction_nonthermal_electrons < 0.5) ? electrons_mixed.av_gammasq()
                                                                                    : (fraction_nonthermal_electrons < 1.)    ? electrons_bpl.av_gammasq()
                                                                                    : (fraction_nonthermal_electrons == 1)    ? electrons_pl.av_gammasq()
                                                                                                                              : throw std::out_of_range("fraction_nonthermal_electrons has to be <=1!");
            Ub = std::pow(magnetic_field, 2.) / (8. * karcst::pi);
            double P_com_P_syn = electron_number_density * karcst::sigtom * radius * average_gamma_squared;
            P_com_P_syn += std::pow(doppler_factor_bulk, 2.) * radiation_energy_density / Ub;


            Lumnorm = karcst::pi * std::pow(radius, 2.) * height * std::pow(doppler_factor_bulk, 4.) *
                      electron_number_density * karcst::sigtom * karcst::cee * average_gamma_squared;
            
            if (verbosity_level > 2)
            {
                std::cout << "P_com / P_syn=" << P_com_P_syn << ", compton?=" << (P_com_P_syn > compton_threshold) << std::endl;
            }
            return (P_com_P_syn > compton_threshold);
        }
    }

    // dp = gamma/p dgamma, alpha = E_photon / (m_e*c^2)
    double RadiationZone::Fic(double gamma, double alpha)
    {
        double gamma_bar = gamma + sqrt(gamma * gamma - 1.0 + 1.0e-10);
        double za = alpha * gamma_bar;
        double zb = alpha / gamma_bar;

        double y = gamma * (Fic_1(za) - Fic_1(zb)) - alpha * (Fic_2(za) - Fic_2(zb));

        return y > 0.0 ? y : 0.0;
    }

    double RadiationZone::Fic_1(double z)
    {
        return z < 0.1 ? 8.0 / 9.0 * z * z * z - 14.0 / 5.0 * z * z * z * z + 196.0 / 25.0 * z * z * z * z * z : (z + 6.0 + 3.0 / z) * log(1.0 + 2.0 * z) - (22.0 * z * z * z / 3.0 + 24.0 * z * z + 18.0 * z + 4.0) / (1.0 + 2.0 * z) / (1.0 + 2.0 * z) - 2.0 + 2.0 * polylog(-2.0 * z);
    }

    double RadiationZone::Fic_2(double z)
    {
        return z < 0.1 ? 4.0 * z * z / 3.0 - 128.0 * z * z * z / 45.0 + 7.0 * z * z * z * z : (z + 31.0 / 6.0 + 5.0 / z + 3.0 / 2.0 / z / z) * log(1.0 + 2.0 * z) - (22.0 * z * z * z / 3.0 + 28.0 * z * z + 103.0 * z / 3.0 + 17.0 + 3.0 / z) / (1.0 + 2.0 * z) / (1.0 + 2.0 * z) - 2.0 + polylog(-2.0 * z);
    }

    double RadiationZone::polylog(double z)
    {
        double y = 0.0;

        if (z >= 0.99)
            throw std::invalid_argument("in function RadiationZone::polylog_2, z>0.9 out of range\n");
        else if (z >= -0.1 and z < 0.9)
            y = polylog_smallz(z);
        else if (z < -0.1 and z > -20.0)
            y = polylog_rk4(z);
        else
            y = polylog_largez(z);
        return y;
    }

    double RadiationZone::polylog_smallz(double z)
    {
        if (z < -0.1 or z > 0.9)
        {
            throw std::invalid_argument("in function RadiationZone::polylog_smallz, argument z too large - requires -0.1<z<0.9 \n");
        }

        return z * (1.0 + z * (1.0 / 4.0 + z * (1.0 / 9.0 + z * (1.0 / 16.0 + z * (1.0 / 25.0 + z / 36.0)))));
    }

    double RadiationZone::polylog_rk4(double z) // applies when z<-0.1
    {
        double x = 0.1;
        double y = 0.0;
        const int N_INTG = 20;
        const double diff_x = (-z - 0.1) / N_INTG;

        if (z > -0.1)
        {
            throw std::invalid_argument("in function RadiationZone::polylog_rk4, argument z should be z<-0.1\n");
        }

        for (int i = 0; i < N_INTG; i++)
        {
            double k1 = -diff_x * log(1.0 + x) / x;
            double k2 = -diff_x * log(1.0 + x + 0.5 * diff_x) / (x + 0.5 * diff_x);
            double k4 = -diff_x * log(1.0 + x + 1.0 * diff_x) / (x + 1.0 * diff_x);
            y += (1.0 / 6.0) * (k1 + k4) + (2.0 / 3.0) * k2;
            x += diff_x;
        }

        return -0.097605235 + y; // the number = Li2[-0.1]
    }

    double RadiationZone::polylog_largez(double z)
    {
        double x = -z;
        return -3.1415926 * 3.1415927 / 6.0 - 0.5 * log(x) * log(x) + 1.0 / x - 1.0 / 4.0 / x / x + 1.0 / 9.0 / x / x / x;
    }

} // namespace bhjet
