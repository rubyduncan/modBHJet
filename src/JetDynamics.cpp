#include "JetDynamics.hpp"
#include <iostream>

namespace bhjet
{

    void JetDynamics::compute_jet_dynamics()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: Computing jet dynamics (empty)" << std::endl;
    }

    std::string JetDynamics::info() const
    {
        return "Base JetDynamics";
    }

    void JetDynamics::reinit_grid_arrays()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: Reset grid arrays" << std::endl;

        z_min_grid = std::vector<double>(n_zones, 0.0);
        z_height_grid = std::vector<double>(n_zones, 0.0);
        z_center_grid = std::vector<double>(n_zones, 0.0);
        radius_grid = std::vector<double>(n_zones, 0.0);
        gamma_grid = std::vector<double>(n_zones, 0.0);
        beta_grid = std::vector<double>(n_zones, 0.0);
        beta_gamma_grid = std::vector<double>(n_zones, 0.0);
        magnetic_field_grid = std::vector<double>(n_zones, 0.0);
        electron_density_grid = std::vector<double>(n_zones, 0.0);
        proton_density_grid = std::vector<double>(n_zones, 0.0);
        electron_temperature_grid = std::vector<double>(n_zones, 0.0);
        proton_temperature_grid = std::vector<double>(n_zones, 0.0);
        fraction_nonthermal_electrons_grid = std::vector<double>(n_zones, 0.0);
        fraction_nonthermal_protons_grid = std::vector<double>(n_zones, 0.0);
        factor_break_electrons_grid = std::vector<double>(n_zones, 0.0);
        factor_break_protons_grid = std::vector<double>(n_zones, 0.0);
        factor_max_energy_electrons_grid = std::vector<double>(n_zones, 0.0);
        factor_max_energy_protons_grid = std::vector<double>(n_zones, 0.0);
        index_injected_electrons_grid = std::vector<double>(n_zones, 0.0);
        index_injected_protons_grid = std::vector<double>(n_zones, 0.0);
        cutoff_type_grid = std::vector<int>(n_zones, defaults::DEFAULT_CUTOFF_TYPE);
    }

    std::vector<double> JetDynamics::get_z_min_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_z_min_grid" << std::endl;
        return z_min_grid;
    }
    std::vector<double> JetDynamics::get_z_height_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_z_height_grid" << std::endl;
        return z_height_grid;
    }
    std::vector<double> JetDynamics::get_z_center_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_z_center_grid" << std::endl;
        return z_center_grid;
    }
    std::vector<double> JetDynamics::get_radius_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_radius_grid" << std::endl;
        return radius_grid;
    }
    std::vector<double> JetDynamics::get_gamma_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_gamma_grid" << std::endl;
        return gamma_grid;
    }
    std::vector<double> JetDynamics::get_beta_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_beta_grid" << std::endl;
        return beta_grid;
    }
    std::vector<double> JetDynamics::get_beta_gamma_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_beta_gamma_grid" << std::endl;
        return beta_gamma_grid;
    }
    std::vector<double> JetDynamics::get_magnetic_field_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_magnetic_field_grid" << std::endl;
        return magnetic_field_grid;
    }
    std::vector<double> JetDynamics::get_electron_density_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_electron_density_grid" << std::endl;
        return electron_density_grid;
    }
    std::vector<double> JetDynamics::get_proton_density_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_proton_density_grid" << std::endl;
        return proton_density_grid;
    }
    std::vector<double> JetDynamics::get_electron_temperature_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_electron_temperature_grid" << std::endl;
        return electron_temperature_grid;
    }
    std::vector<double> JetDynamics::get_proton_temperature_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_proton_temperature_grid" << std::endl;
        return proton_temperature_grid;
    }
    std::vector<double> JetDynamics::get_fraction_nonthermal_electrons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_fraction_nonthermal_electrons_grid" << std::endl;
        return fraction_nonthermal_electrons_grid;
    }
    std::vector<double> JetDynamics::get_fraction_nonthermal_protons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_fraction_nonthermal_protons_grid" << std::endl;
        return fraction_nonthermal_protons_grid;
    }
    std::vector<double> JetDynamics::get_factor_break_electrons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_factor_break_electrons_grid" << std::endl;
        return factor_break_electrons_grid;
    }
    std::vector<double> JetDynamics::get_factor_break_protons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_factor_break_protons_grid" << std::endl;
        return factor_break_protons_grid;
    }
    std::vector<double> JetDynamics::get_factor_max_energy_electrons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_factor_max_energy_electrons_grid" << std::endl;
        return factor_max_energy_electrons_grid;
    }
    std::vector<double> JetDynamics::get_factor_max_energy_protons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_factor_max_energy_protons_grid" << std::endl;
        return factor_max_energy_protons_grid;
    }
    std::vector<double> JetDynamics::get_index_injected_electrons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_index_injected_electrons_grid" << std::endl;
        return index_injected_electrons_grid;
    }
    std::vector<double> JetDynamics::get_index_injected_protons_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_index_injected_protons_grid" << std::endl;
        return index_injected_protons_grid;
    }
    std::vector<int> JetDynamics::get_cutoff_type_grid()
    {
        if (verbosity_level > 1)
            std::cout << "JetDynamics: get_cutoff_type_grid" << std::endl;
        return cutoff_type_grid;
    }

} // namespace bhjet
