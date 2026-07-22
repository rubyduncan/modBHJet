#include "BHJet.hpp"
#include "JetDynamics.hpp"
#include "RadiationZone.hpp"
#include "utils.hpp"
#include <iostream>
#include <chrono>
#include "kariba/constants.hpp"
namespace karcst = kariba::constants;

namespace bhjet
{

    void BHJet::init_jet_dynamics(std::shared_ptr<JetDynamics> jet_dynamics_)
    {
        if (verbosity_level > 1)
            std::cout << "Initialising jet dynamics object" << std::endl;
        jet_dynamics = jet_dynamics_;
    }
    std::shared_ptr<JetDynamics> BHJet::get_jet_dynamics()
    {
        return jet_dynamics;
    }
    void BHJet::add_target_photon_field(std::shared_ptr<TargetPhotonField> target)
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: adding target photon field: " << target->name << std::endl;
        target_list.push_back(target);
    }

    void BHJet::remove_target_photon_field(std::shared_ptr<TargetPhotonField> target)
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: removing target photon field: " << target->name << std::endl;

        for (size_t i=0;  i<target_list.size(); i++){
            if(target_list[i]->name == target->name) {
                target_list.erase(target_list.begin() + i);
                std::cout << "BHJet: removed target photon field "<< i << ": " << target->name << std::endl;
            }
        }
    }

    void BHJet::clear_target_photon_fields()
    {
        if (verbosity_level > 1)
            std::cout << "BHJet: clearing all target photon fields " << std::endl;
        target_list = std::vector<std::shared_ptr<TargetPhotonField>>();
    }
    

    void BHJet::compute_full_jet(
        std::vector<double> photon_energy_grid)
    {
        if (verbosity_level > 1)
            std::cout << "Computing full jet (dynamics + radiation)" << std::endl;
        auto tstart = std::chrono::steady_clock::now();

        // compute jet dynamics
        jet_dynamics->compute_jet_dynamics();
        // after compute_jet_dynamics() to update n_zones first
        computation_times = std::vector<double>(jet_dynamics->n_zones + 1, 0.0);

        if (profile_time)
            computation_times[0] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - tstart).count();


        // init zones and compute radiation
        observed_photon_energy_grid = photon_energy_grid;
        size_t n_bins_phot = observed_photon_energy_grid.size();
        observed_photon_flux_total = std::vector<double>(n_bins_phot, 0.0);
        radiation_zones = std::vector<RadiationZone>(jet_dynamics->n_zones);
        std::vector<double> z_center_grid = jet_dynamics->get_z_center_grid();
        std::vector<double> radius_grid = jet_dynamics->get_radius_grid();
        std::vector<double> z_height_grid = jet_dynamics->get_z_height_grid();
        std::vector<double> beta_gamma_grid = jet_dynamics->get_beta_gamma_grid();
        std::vector<double> magnetic_field_grid = jet_dynamics->get_magnetic_field_grid();
        std::vector<double> electron_density_grid = jet_dynamics->get_electron_density_grid();
        std::vector<double> proton_density_grid = jet_dynamics->get_proton_density_grid();
        std::vector<double> electron_temperature_grid = jet_dynamics->get_electron_temperature_grid();
        std::vector<double> proton_temperature_grid = jet_dynamics->get_proton_temperature_grid();
        std::vector<double> fraction_nonthermal_electrons_grid = jet_dynamics->get_fraction_nonthermal_electrons_grid();
        std::vector<double> fraction_nonthermal_protons_grid = jet_dynamics->get_fraction_nonthermal_protons_grid();
        std::vector<double> factor_break_electrons_grid = jet_dynamics->get_factor_break_electrons_grid();
        std::vector<double> factor_break_protons_grid = jet_dynamics->get_factor_break_protons_grid();
        std::vector<double> factor_max_energy_electrons_grid = jet_dynamics->get_factor_max_energy_electrons_grid();
        std::vector<double> factor_max_energy_protons_grid = jet_dynamics->get_factor_max_energy_protons_grid();
        std::vector<double> index_injected_electrons_grid = jet_dynamics->get_index_injected_electrons_grid();
        std::vector<double> index_injected_protons_grid = jet_dynamics->get_index_injected_protons_grid();
        std::vector<int> cutoff_type_grid = jet_dynamics->get_cutoff_type_grid();
        bool compton_switch = true;
        // double compton_threshold = 1e-5;
        // Each iteration writes only to its own RadiationZone.
#pragma omp parallel for schedule(dynamic)
        for (std::ptrdiff_t zone = 0; zone < static_cast<std::ptrdiff_t>(jet_dynamics->n_zones); zone++)
        {
            size_t i = static_cast<size_t>(zone);
            auto zone_start = std::chrono::steady_clock::now();
            // always do compton emission in first two zones
            bool force_compton = (i < 2);
            if (verbosity_level > 1)
#pragma omp critical(bhjet_zone_logging)
                std::cout << "Computing zone " << i << std::endl;
            radiation_zones[i] = RadiationZone(
                magnetic_field_grid[i],
                radius_grid[i],
                z_height_grid[i],
                "cylinder",
                beta_gamma_grid[i],
                theta_obs, distance, redshift,
                electron_density_grid[i],
                proton_density_grid[i],
                electron_temperature_grid[i],
                proton_temperature_grid[i],
                fraction_nonthermal_electrons_grid[i],
                fraction_nonthermal_protons_grid[i],
                factor_break_electrons_grid[i],
                factor_break_protons_grid[i],
                factor_max_energy_electrons_grid[i],
                factor_max_energy_protons_grid[i],
                index_injected_electrons_grid[i],
                index_injected_protons_grid[i],
                include_counterjet, force_compton, compton_switch,
                compton_threshold, profile_time, verbosity_level,
                cutoff_type_grid[i]);
            
            if (verbosity_level > 2)
#pragma omp critical(bhjet_zone_logging)
                std::cout << "Adding targets in zone " << i << std::endl;
            for (size_t b = 0; b < target_list.size(); b++)
            {
                auto [energies, energy_densities] = target_list[b]->get_target_energy_grid_and_density(
                    z_center_grid[i], beta_gamma_grid[i], theta_obs,
                    observed_photon_energy_grid[0], observed_photon_energy_grid[observed_photon_energy_grid.size()-1]
                );
                radiation_zones[i].add_target_photon_field(energies, energy_densities, target_list[b]->name);
            }

            radiation_zones[i].compute_particles();
            radiation_zones[i].compute_radiation(observed_photon_energy_grid);

            if (profile_time)
                computation_times[i + 1] = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - zone_start).count();
        }

        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            // Keep total emission serial to preserve the accumulation order.
            add_emission_on_interpolated_grid(
                radiation_zones[i].get_observed_photon_energy_grid_total(),
                radiation_zones[i].get_observed_photon_flux_total(),
                observed_photon_energy_grid, observed_photon_flux_total);
        }
        // black bodies
        for (size_t b = 0; b < target_list.size(); b++)
        {
            if(target_list[b]->add_to_total_flux){
                target_list[b]->update_observed_flux();
                add_emission_on_interpolated_grid(
                    target_list[b]->get_observed_energy(), target_list[b]->get_observed_number_flux(),
                    observed_photon_energy_grid, observed_photon_flux_total);
            }
        }
    }

    std::vector<double> BHJet::get_observed_photon_integrated_flux_total(double z_min, double z_max)
    {
        std::vector<double> photon_flux_obs_cum(observed_photon_energy_grid.size(), 0.0);
        std::vector<double> z_values = jet_dynamics->get_z_min_grid();
        std::vector<double> h_values = jet_dynamics->get_z_height_grid();
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            if ((z_values[i] >= z_min) && (z_values[i] + h_values[i] <= z_max))
            {

                add_emission_on_interpolated_grid(
                    radiation_zones[i].get_observed_photon_energy_grid_total(),
                    radiation_zones[i].get_observed_photon_flux_total(),
                    observed_photon_energy_grid, photon_flux_obs_cum);
            }
        }
        return photon_flux_obs_cum;
    }

    std::vector<double> BHJet::get_observed_photon_integrated_flux_electron_cyclosyn(double z_min, double z_max)
    {
        std::vector<double> photon_flux_obs_cum(observed_photon_energy_grid.size(), 0.0);
        std::vector<double> z_values = jet_dynamics->get_z_min_grid();
        std::vector<double> h_values = jet_dynamics->get_z_height_grid();
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            if ((z_values[i] >= z_min) && (z_values[i] + h_values[i] <= z_max))
            {

                add_emission_on_interpolated_grid(
                    radiation_zones[i].get_observed_photon_energy_grid_electron_cyclosyn(),
                    radiation_zones[i].get_observed_photon_flux_electron_cyclosyn(),
                    observed_photon_energy_grid, photon_flux_obs_cum);
            }
        }
        return photon_flux_obs_cum;
    }

    std::vector<double> BHJet::get_observed_photon_integrated_flux_electron_compton(double z_min, double z_max)
    {
        std::vector<double> photon_flux_obs_cum(observed_photon_energy_grid.size(), 0.0);
        std::vector<double> z_values = jet_dynamics->get_z_min_grid();
        std::vector<double> h_values = jet_dynamics->get_z_height_grid();
        for (size_t i = 0; i < jet_dynamics->n_zones; i++)
        {
            if ((z_values[i] >= z_min) && (z_values[i] + h_values[i] <= z_max))
            {

                add_emission_on_interpolated_grid(
                    radiation_zones[i].get_observed_photon_energy_grid_electron_compton(),
                    radiation_zones[i].get_observed_photon_flux_electron_compton(),
                    observed_photon_energy_grid, photon_flux_obs_cum);
            }
        }
        return photon_flux_obs_cum;
    }

    std::vector<double> BHJet::get_observed_photon_energy_grid()
    {
        return observed_photon_energy_grid;
    }
    std::vector<double> BHJet::get_observed_photon_flux_total()
    {
        return observed_photon_flux_total;
    }
    std::vector<double> BHJet::get_observed_target_photon_energy(std::string name)
    {
        for (size_t b = 0; b < target_list.size(); b++)
        {
            if (name == target_list[b]->name)
            {
                return target_list[b]->get_observed_energy();
            }
        }
        throw std::out_of_range("target photon field called " + name + " not found!");
    }
    std::vector<double> BHJet::get_observed_target_photon_flux(std::string name)
    {
        for (size_t b = 0; b < target_list.size(); b++)
        {
            if (name == target_list[b]->name)
            {
                return target_list[b]->get_observed_energy_flux();
            }
        }
        throw std::out_of_range("target photon field called " + name + " not found!");
    }


    std::vector<double> BHJet::get_computation_times()
    {
        return computation_times;
    }

} // namespace bhjet
