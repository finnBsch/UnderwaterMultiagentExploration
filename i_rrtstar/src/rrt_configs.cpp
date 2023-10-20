//
// Created by finn on 6/15/23.
//
#include "rrt_configs.h"
#include <iostream>
#include <toml++/toml.h>
#include <string>
#include "res_loc.h"


template<typename T> T read_toml_value(toml::node_view<toml::node> node_view){
    auto optional_v = node_view.value<T>();
    if(optional_v){
        return optional_v.value();
    }
    else{
        std::cerr << "TOML field not populated or incorrect datatype." <<
                  std::endl;
        return T();
    }
}

RRTParams::RRTParams(): RRTParams(RES_LOC "/configs/rrt_conf.toml") {
}

void RRTParams::exportConfig(std::string path) {
    std::string source_str = "";
    // General
    source_str += "[general]\n";
    source_str += "batch_size = " + std::to_string(batch_size) + "\n";
    source_str += "draw_all = " + std::to_string(draw_all) + "\n";
    source_str += "num_nodes = " + std::to_string(num_nodes) + "\n";
    source_str += "expand_radius = " + std::to_string(expand_radius) + "\n";
    source_str += "rewire_radius = " + std::to_string(rewire_radius) + "\n";
    source_str += "size_x = " + std::to_string(size_x) + "\n";
    source_str += "size_y = " + std::to_string(size_y) + "\n";
    source_str += "offset_replan = " + std::to_string(offset_replan) + "\n";
    source_str += "min_tail_length = " + std::to_string(min_tail_length) + "\n";
    source_str += "maximize_util = " + std::to_string(maximize_util) + "\n";

    // Information
    source_str += "\n[information]\n";
    source_str += "gamma = " + std::to_string(info_p.gamma) + "\n";
    source_str += "alpha = " + std::to_string(info_p.alpha) + "\n";
    source_str += "weight_mean = " + std::to_string(info_p.weight_mean) + "\n";
    source_str += "weight_source = " + std::to_string(info_p.weight_source) +
                  "\n";
    source_str += "discount_info = " + std::to_string(info_p.discount_info) +
                  "\n";

    // Path Cost
    source_str += "\n[path_cost]\n";
    source_str += "weight_distance = " + std::to_string(cost_p.weight_distance) +
                  "\n";
    source_str += "weight_steering = " + std::to_string(cost_p.weight_steering) +
                  "\n";

    // Export toml file
    std::ofstream out(path + "/configs/rrt_conf.toml");
    out << source_str;
    out.close();
}

RRTParams::RRTParams(std::string path) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(path);
//        std::cout << "Loading RRT config: " << std::endl;
//        std::cout << tbl << "\n" << "\n";

        // Load RRT params
        batch_size = read_toml_value<int>(tbl["general"]["batch_size"]);
        draw_all = read_toml_value<bool>(tbl["general"]["draw_all"]);
        num_nodes = read_toml_value<int>(tbl["general"]["num_nodes"]);
        expand_radius = read_toml_value<float>(tbl["general"]["expand_radius"]);
        rewire_radius = read_toml_value<float>(tbl["general"]["rewire_radius"]);
        size_x = read_toml_value<float>(tbl["general"]["size_x"]);
        size_y = read_toml_value<float>(tbl["general"]["size_y"]);
        offset_replan = read_toml_value<int>(tbl["general"]["offset_replan"]);
        min_tail_length = read_toml_value<int>(tbl["general"]["min_tail_length"]);
        maximize_util = read_toml_value<bool>(tbl["general"]["maximize_util"]);

        // Load Information params
        info_p.gamma = read_toml_value<float>
                (tbl["information"]["gamma"]);
        info_p.alpha = read_toml_value<float>
                (tbl["information"]["alpha"]);
        info_p.weight_mean = read_toml_value<float>
                (tbl["information"]["weight_mean"]);
        info_p.weight_source = read_toml_value<float>
                (tbl["information"]["weight_source"]);
        info_p.discount_info = read_toml_value<float>
                (tbl["information"]["discount_info"]);

        // Load PathCost params
        cost_p.weight_distance = read_toml_value<float>
                (tbl["path_cost"]["weight_distance"]);
        cost_p.weight_steering = read_toml_value<float>
                (tbl["path_cost"]["weight_steering"]);
    }
    catch (const toml::parse_error& err)
    {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
    rewire_radius_sq = powf(rewire_radius, 2);
}

