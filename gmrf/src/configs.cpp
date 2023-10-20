//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//
#include "gmrf_configs.h"
#include <iostream>
#include <toml++/toml.h>
#include "res_loc.h"

GmrfParams::GmrfParams(): GmrfParams(RES_LOC "/configs/gmrf_conf.toml") {
}

void GmrfParams::exportConfig(std::string path) {
    std::string source_str = "";
    // General
    source_str += "[general]\n";
    source_str += "size_x = " + std::to_string(size_x) + "\n";
    source_str += "size_y = " + std::to_string(size_y) + "\n";
    source_str += "N_X = " + std::to_string(N_X) + "\n";
    source_str += "N_Y = " + std::to_string(N_Y) + "\n";
    // Field
    source_str += "\n[field_params]\n";
    source_str += "prior_mean = " + std::to_string(prior_mean) + "\n";
    source_str += "prior_sigma = " + std::to_string(prior_sigma) + "\n";
    source_str += "sigma_zeta = " + std::to_string(sigma_zeta) + "\n";
    source_str += "sigma_p = " + std::to_string(sigma_p_per_m) + "\n";
    source_str += "sigma_s = " + std::to_string(sigma_s) + "\n";
    source_str += "sigma_flow = " + std::to_string(sigma_flow) + "\n";
    // Flow Field
    source_str += "\n[flow_field_params]\n";
    source_str += "prior_mean_x = " + std::to_string(ff_prior_mean_x) + "\n";
    source_str += "prior_mean_y = " + std::to_string(ff_prior_mean_y) + "\n";
    source_str += "prior_sigma = " + std::to_string(ff_prior_sigma) + "\n";
    source_str += "sigma_zeta = " + std::to_string(ff_sigma_zeta) + "\n";
    source_str += "sigma_p = " + std::to_string(ff_sigma_p_per_m) + "\n";
    source_str += "sigma_s = " + std::to_string(ff_sigma_s) + "\n";
    source_str += "sigma_c = " + std::to_string(ff_sigma_c) + "\n";
    source_str += "sigma_o = " + std::to_string(ff_sigma_o) + "\n";
    // Energy
    source_str += "\n[energies]\n";
    source_str += "field_neighbors = " + std::to_string(field_neighbors) + "\n";
    source_str += "flow_neighbors = " + std::to_string(flow_neighbors) + "\n";
    source_str += "flow_conservation = " + std::to_string(flow_conservation) + "\n";
    source_str += "obstacle_flow = " + std::to_string(obstacle_flow) + "\n";
    source_str += "flow_field = " + std::to_string(flow_field) + "\n";
    // Source loc
    source_str += "\n[source_loc]\n";
    source_str += "all_neighbors = " + std::to_string(all_neighbors) + "\n";

    std::ofstream out(path + "/configs/gmrf_conf.toml");
    out << source_str;
    out.close();

}

GmrfParams::GmrfParams(std::string path) {
    toml::table tbl;
    try
    {
        tbl = toml::parse_file(path);
//        std::cout << "Loading GMRF config: " << std::endl;
//        std::cout << tbl << "\n" << "\n";

        // Load general params
        size_x = tbl["general"]["size_x"].value_or(10.0f);
        size_y = tbl["general"]["size_y"].value_or(10.0f);
        N_X = tbl["general"]["N_X"].value_or(10);
        N_Y = tbl["general"]["N_Y"].value_or(10);

        // Load field params
        prior_mean = tbl["field_params"]["prior_mean"].value_or(0.0f);
        prior_sigma = tbl["field_params"]["prior_sigma"].value_or(0.0f);
        sigma_zeta = tbl["field_params"]["sigma_zeta"].value_or(0.0f);
        sigma_p_per_m = tbl["field_params"]["sigma_p"].value_or(0.0f);
        sigma_s = tbl["field_params"]["sigma_s"].value_or(0.0f);
        sigma_flow = tbl["field_params"]["sigma_flow"].value_or(0.0f);
        // sigma_p = sigma_p_per_m * size_x / (N_X - 1);
        // Load flow field params
        ff_prior_mean_x = tbl["flow_field_params"]["prior_mean_x"].value_or(0.0f);
        ff_prior_mean_y = tbl["flow_field_params"]["prior_mean_y"].value_or(0.0f);
        ff_prior_sigma = tbl["flow_field_params"]["prior_sigma"].value_or(0.0f);
        ff_sigma_zeta = tbl["flow_field_params"]["sigma_zeta"].value_or(0.0f);
        ff_sigma_p_per_m = tbl["flow_field_params"]["sigma_p"].value_or(0.0f);
        ff_sigma_s = tbl["flow_field_params"]["sigma_s"].value_or(0.0f);
        ff_sigma_c = tbl["flow_field_params"]["sigma_c"].value_or(0.0f);
        ff_sigma_o = tbl["flow_field_params"]["sigma_o"].value_or(0.0f);
        // ff_sigma_p = ff_sigma_p_per_m * size_x / (N_X - 1);

        // Load Energy config
        field_neighbors = tbl["energies"]["field_neighbors"].value_or(true);
        flow_neighbors = tbl["energies"]["flow_neighbors"].value_or(true);
        flow_conservation = tbl["energies"]["flow_conservation"].value_or(true);
        obstacle_flow = tbl["energies"]["obstacle_flow"].value_or(true);
        flow_field = tbl["energies"]["flow_field"].value_or(true);

        // Load source loc config
        all_neighbors = tbl["source_loc"]["all_neighbors"].value_or(true);
    }
    catch (const toml::parse_error& err)
    {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}

Viz_Params::Viz_Params() {
    toml::table tbl;
    try
    {
        tbl = toml::parse_file(RES_LOC "/configs/viz_conf.toml");
        std::cout << "Loading Viz config: " << std::endl;
        std::cout << tbl << "\n" << "\n";
        res_x = tbl["viz_params"]["res_x"].value_or(1920);
        margin = tbl["viz_params"]["margin"].value_or(0.01f);
        infobox = tbl["viz_params"]["infobox"].value_or(0.1f);
        draw_map = tbl["viz_params"]["draw_map"].value_or(true);
        draw_uncertainty = tbl["viz_params"]["draw_uncertainty"].value_or(true);
        draw_gt = tbl["viz_params"]["draw_groundtruth"].value_or(true);
    }
    catch (const toml::parse_error& err)
    {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}
