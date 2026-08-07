#include "IsoJet.hpp"

#include <array>
#include <iostream>
#include "kariba/Thermal.hpp"
#include "kariba/constants.hpp"

namespace karcst = kariba::constants;

namespace bhjet
{

    void IsoJet::compute_jet_dynamics()
    {
        if (verbosity_level > 1)
            std::cout << "Computing IsoJet dynamics" << std::endl;

        // clean up existing gsl interpolation
        if (verbosity_level > 2)
            std::cout << "cleaning up spline memory" << std::endl;
        if (spline_speed)
            gsl_spline_free(spline_speed);


        // init basic internal params
        eddington_luminosity = 1.25e38 * mass_bh;
        r_g = karcst::gconst * mass_bh * karcst::msun / (karcst::cee * karcst::cee);

        jet_dyn.min = z_jet_launching * r_g;
        jet_dyn.max = z_max_calculation * r_g;
        jet_dyn.h0 = 2. * r_initial * r_g + z_jet_launching * r_g;
        jet_dyn.r0 = r_initial * r_g;
        jet_dyn.beta0 = sqrt(4. / 3. * (4. / 3. - 1.) /
                             (4. / 3. + 1.)); // set initial jet speed for relativistic fluid, g=4/3
        jet_dyn.gam0 = 1. / sqrt(1. - (std::pow(jet_dyn.beta0,
                                                2.))); // set corresponding lorentz factor
        jet_dyn.Rg = r_g;

        // complex way to calculate integral over thermal distribution using kariba
        // Dummy particle distribution, needed for average lorentz factor in
        // equipartition function The number density is just set to unity, the
        // normalisation is not needed to calculate the average Lorenz factor of the
        // thermal distribution anyway
        kariba::Thermal dummy_elec(70); // use 70 bins for integral
        dummy_elec.set_temp_kev(electron_temperature_jet_base);
        dummy_elec.set_p();
        dummy_elec.set_norm(1.);
        dummy_elec.set_ndens();

        nozzle_ener.pbeta = plasma_beta_jet_base;
        nozzle_ener.Nj = jet_power_eddington * eddington_luminosity;

        nozzle_ener.av_gamma = dummy_elec.av_gamma();

        // fill the spline for the speed profile (formerly "velprof_iso")
        calc_velocity_profile_iso_jet();
        // fill the other values of nozzle_ener (formerly "equipartition")
        calc_nozzle_energetics_equipartition();

        // grid building
        // double dlgz = 0.1;
        double zmin = jet_dyn.min;
        double zmax = jet_dyn.max;
        double zdiss = z_dissipation*r_g;

        if (zdiss < zmin*pow(10, dlgz*1.5)) {
            zdiss = zmin*pow(10, dlgz/2);
            if (verbosity_level > 1)
                std::cout << "z_dissipation<z_jet_launching, set z_dissipation=z_jet_launching*10^(dlgz/2)" << std::endl;
            
        } else if (zdiss > zmax*pow(10, -dlgz*1.5)) {
            zdiss = zmax*pow(10, -dlgz/2);
            if (verbosity_level > 1)
                std::cout << "z_dissipation>z_max_calculation, set z_dissipation=z_max_calculation*10^(-dlgz/2)" << std::endl;
        }

        n_zones= static_cast<size_t>(std::log10(zmax/zmin)/dlgz);
        size_t N1 = static_cast<size_t>(std::log10(zdiss*pow(10, -dlgz/2)/zmin)/dlgz);
        double dlgz1 = std::log10(zdiss*pow(10, -dlgz/2)/zmin)/static_cast<double>(N1);
        size_t N2 = n_zones - N1 - 1;
        double dlgz2 = std::log10(zmax/(zdiss*pow(10, dlgz/2)))/static_cast<double>(N2);
        
        reinit_grid_arrays();
        // std::cout<< "Nzones=" << n_zones << ", N1=" << N1 << std::endl;

        for (size_t i = 0; i < n_zones; i++)
        {
            // std::cout<< "zone " << i << std::endl;
            // calculate the size of the next zone (formerly "jetgrid")
            // fills z_min_grid and z_height_grid
            // calc_grid_next_zone(i, cut, zcut);
            if(i<N1) {
                z_min_grid[i] = zmin*pow(10, dlgz1 * i);
                z_center_grid[i] = zmin*pow(10, dlgz1 * (i+0.5));
                z_height_grid[i] = zmin*(pow(10, dlgz1 * (i+1)) - pow(10, dlgz1 * i));
            } else if(i==N1){
                z_min_grid[i] = zdiss*pow(10, -dlgz/2);
                z_center_grid[i] = zdiss;
                z_height_grid[i] = zdiss*(pow(10, dlgz/2) - pow(10, -dlgz/2));
            } else {
                z_min_grid[i] = zdiss*pow(10, dlgz/2 + dlgz2 * (i-N1-1));
                z_center_grid[i] = zdiss*pow(10, dlgz/2 +  dlgz2 * (i-N1-1+0.5));
                z_height_grid[i] = zdiss*(pow(10, dlgz/2 + dlgz2 * (i-N1-1+1)) - pow(10, dlgz/2 + dlgz2 * (i-N1-1)));
            }

            // calculate the parameters of the cell (formerly "isojetpars")
            calc_zone_properties(i);
            // isojetpars(z, jet_dyn, nozzle_ener, tshift, zone, spline_speed, acc_speed);
        }
    }

    
    // Velocity profile function: interpolate 1-d hydro jet velocity tables
    // profile returns x axis in either units of initial jet radius 
    // and y axis in units of gamma*beta, where
    // gamma is the jet bulk lorentz factor and beta its speed in units of c For
    // information on the 2 velocity profiles below see Crumley et al. 2016
    void IsoJet::calc_velocity_profile_iso_jet()
    {
        if (verbosity_level > 1)
            std::cout << "IsoJet::calc_velocity_profile_iso_jet" << std::endl;

        // Tabulated velocity for 1D quasi-isothermal Bernoulli eq.
        const size_t n_bins_speed = 54;

        if (verbosity_level > 2)
            std::cout << "allocating spline memory" << std::endl;
        spline_speed_accel = gsl_interp_accel_alloc();
        spline_speed = gsl_spline_alloc(gsl_interp_steffen, n_bins_speed);


        std::array<double, n_bins_speed> gbx_vel_iso = {
            1.0,         1.00001,     1.00005,     1.00023,     1.00101,     1.00456,     1.02053,
            1.09237,     1.41567,     2.87053,     6.26251,     14.3691,     34.0825,     82.8831,
            205.572,     518.283,     1325.11,     3429.81,     8975.1,      23719.5,     63255.,
            170099.,     460962.,     1.25824e+06, 3.4578e+06,  9.5633e+06,  2.66095e+07, 7.44655e+07,
            2.09528e+08, 5.92636e+08, 1.6846e+09,  4.81148e+09, 1.38055e+10, 3.9787e+10,  1.15153e+11,
            3.34651e+11, 9.76408e+11, 2.85981e+12, 8.40728e+12, 2.4805e+13,  7.34419e+13, 2.18185e+14,
            6.50341e+14, 1.94471e+15, 5.83352e+15, 1.75523e+16, 5.29701e+16, 1.60321e+17, 4.86616e+17,
            1.48111e+18, 4.52032e+18, 1.38326e+19, 4.24394e+19, 1.0e+20};

        std::array<double, n_bins_speed> gby_vel_iso = {
            0.485071, 1.05031, 1.05032, 1.05039, 1.05067, 1.05193, 1.05751, 1.08105, 1.16389,
            1.35278,  1.52406, 1.68077, 1.82495, 1.95888, 2.08429, 2.20255, 2.3147,  2.42158,
            2.52386,  2.62205, 2.71662, 2.80793, 2.89628, 2.98195, 3.06516, 3.14611, 3.22497,
            3.30189,  3.37702, 3.45046, 3.52232, 3.59270, 3.66169, 3.72937, 3.79580, 3.86106,
            3.92520,  3.98827, 4.05033, 4.11143, 4.17160, 4.23089, 4.28933, 4.34696, 4.40381,
            4.45992,  4.51530, 4.56999, 4.62402, 4.67740, 4.73015, 4.78230, 4.83388, 4.87281};

        gsl_spline_init(spline_speed, gbx_vel_iso.data(), gby_vel_iso.data(), n_bins_speed);
    }

    // Equipartition functions: calculate bfield,lepton number density,proton number
    // density at the base for given jet power, jet base radius, initial speed,
    // initial plasma beta, accounting for 1 or 2 jets.
    void IsoJet::calc_nozzle_energetics_equipartition()
    {
        if (verbosity_level > 1)
            std::cout << "IsoJet::calc_nozzle_energetics_equipartition" << std::endl;

        double eq_fac, dyn_fac;    // the two numbers that change equipartition are
                                   // the jet dynamics and equipartition assumptions
        size_t npsw = 1; // hardcoded in bhjet anyways
        if (npsw == 0) {
            // 0: no protons
            eq_fac = nozzle_ener.av_gamma * karcst::emerg * (1. + 1. / nozzle_ener.pbeta);
            dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.beta0 * jet_dyn.gam0 * karcst::cee;
            nozzle_ener.lepdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
            nozzle_ener.protdens = 0;
            nozzle_ener.bfield = std::pow(8. * karcst::pi * nozzle_ener.av_gamma * nozzle_ener.lepdens * karcst::emgm *
                                    std::pow(karcst::cee, 2.) / nozzle_ener.pbeta,
                                1. / 2.);
        } else if (npsw == 1) {
            // 1: Up = Ue+Ub
            eq_fac = 2. * nozzle_ener.av_gamma * karcst::emerg * (1. + 1. / nozzle_ener.pbeta);
            dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.beta0 * jet_dyn.gam0 * karcst::cee;
            nozzle_ener.lepdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
            nozzle_ener.protdens =
                (1. + 1. / nozzle_ener.pbeta) * nozzle_ener.av_gamma * nozzle_ener.lepdens * (karcst::emgm / karcst::pmgm);
            nozzle_ener.bfield = std::pow(8. * karcst::pi * nozzle_ener.av_gamma * nozzle_ener.lepdens * karcst::emgm *
                                    std::pow(karcst::cee, 2.) / nozzle_ener.pbeta,
                                1. / 2.);
        } else if (npsw == 2) {
            // 2: ne = np
            eq_fac = karcst::emerg * (karcst::pmgm / karcst::emgm + nozzle_ener.av_gamma * (1. + 1. / nozzle_ener.pbeta));
            dyn_fac = 2. * karcst::pi * std::pow(jet_dyn.r0, 2.) * jet_dyn.beta0 * jet_dyn.gam0 * karcst::cee;
            nozzle_ener.lepdens = nozzle_ener.Nj / (eq_fac * dyn_fac);
            nozzle_ener.protdens = nozzle_ener.lepdens;
            nozzle_ener.bfield = std::pow(8. * karcst::pi * nozzle_ener.av_gamma * nozzle_ener.lepdens * karcst::emgm *
                                    std::pow(karcst::cee, 2.) / nozzle_ener.pbeta,
                                1. / 2.);
        }
        nozzle_ener.eta = nozzle_ener.lepdens / nozzle_ener.protdens;
        nozzle_ener.sig0 = std::pow(nozzle_ener.bfield, 2.) /
                (4. * karcst::pi * nozzle_ener.protdens * karcst::pmgm * std::pow(karcst::cee, 2.));
    }

    // isojetpars
    void IsoJet::calc_zone_properties(size_t i)
    {
        if (verbosity_level > 1)
            std::cout << "IsoJet::calc_zone_properties " << i << std::endl;

        // Note on z_eval: for the old agnjet, the velocity is given as a funct ion
        // of jet launching point min (see Crumley et al. 2016)
        double z_in_r0 = (std::max(z_center_grid[i] - jet_dyn.h0, 0.) + jet_dyn.r0) / jet_dyn.r0;
        double gb, mj;
        double gb0 = jet_dyn.gam0 * jet_dyn.beta0;
        if (z_center_grid[i] < jet_dyn.h0)
        {
            gb = gb0;
        }
        else
        {
            gb = gsl_spline_eval(spline_speed, z_in_r0, spline_speed_accel);
        }
        // std::cout << "z_min_grid[i]= " << z_min_grid[i] << std::endl;
        // std::cout << "z_center_grid[i]= " << z_center_grid[i] << std::endl;
        // std::cout << "z_height_grid[i]= " << z_height_grid[i] << std::endl;

        mj = gb / gb0;

        beta_gamma_grid[i] = gb;
        gamma_grid[i] = std::sqrt(std::pow(gb, 2.) + 1.);
        beta_grid[i] = std::sqrt((std::pow(gamma_grid[i], 2.) - 1.) / std::pow(gamma_grid[i], 2.));

        radius_grid[i] = jet_dyn.r0 + std::max(z_center_grid[i] - jet_dyn.h0, 0.) / mj;

        electron_density_grid[i] = nozzle_ener.lepdens * (
            std::pow(jet_dyn.r0 / radius_grid[i], 2.) / mj);
        // Todo: fill here proton density and temperature arrays
        proton_density_grid[i] = nozzle_ener.protdens * (
            std::pow(jet_dyn.r0 / radius_grid[i], 2.) / mj); // check this
        proton_temperature_grid[i] = electron_temperature_jet_base; // check this

        if (z_center_grid[i] < jet_dyn.h0) {
            magnetic_field_grid[i] = nozzle_ener.bfield * (
                jet_dyn.r0 / radius_grid[i]) / std::pow(mj, 0.5);
            electron_temperature_grid[i] = electron_temperature_jet_base;
        } else {
            magnetic_field_grid[i] = nozzle_ener.bfield * (
                jet_dyn.r0 / radius_grid[i]) / std::pow(mj, 0.5 + 1. / 6.);
            electron_temperature_grid[i] = electron_temperature_jet_base / std::pow(mj, 1. / 3.);
        }
        
        // same as for bljet
        if (z_center_grid[i] < z_dissipation * r_g)
        {
            fraction_nonthermal_electrons_grid[i] = 0.0;
            fraction_nonthermal_protons_grid[i] = 0.0;
        }
        else
        {
            fraction_nonthermal_electrons_grid[i] = fraction_nonthermal_electrons;
            fraction_nonthermal_protons_grid[i] = fraction_nonthermal_protons;
        }
        factor_break_electrons_grid[i] = factor_break_electrons;
        factor_break_protons_grid[i] = factor_break_protons;
        factor_max_energy_electrons_grid[i] = factor_max_energy_electrons;
        factor_max_energy_protons_grid[i] = factor_max_energy_protons;
        index_injected_electrons_grid[i] = index_injected_electrons;
        index_injected_protons_grid[i] = index_injected_protons;
        cutoff_type_grid[i] = cutoff_type;
    }

} // namespace bhjet
