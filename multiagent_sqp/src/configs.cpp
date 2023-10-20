//
// Created by finn on 5/28/23.
//
#include "configs.h"
#include <iostream>
#include <toml++/toml.h>
#include <string>
#include "bicycle_dynamics.h"
#include "hippocampus_dynamics.h"
#include "bluerov_dynamics.h"

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

/**
 * Constructor for SQPParams. Loads the config file and populates the
 * parameters. If originally loading just one agent, not all weights will be
 * loaded.
 */
SQPParams::SQPParams():SQPParams(RES_LOC "/configs/sqp_conf.toml") {
}

/**
 * Export config to toml file.
 * @param path
 */
void SQPParams::exportConfig(std::string path) {
    std::string source_str = "";
    // General Params
    source_str += "[general]\n";
    source_str += "N_horizon = " + std::to_string(N_horizon) + "\n";
    source_str += "num_agents = " + std::to_string(num_agents) + "\n";

    source_str += "agent_type = \"" + convert(agent_type) + "\"\n";
    source_str += "sqp_iters = " + std::to_string(sqp_iters) + "\n";
    source_str += "alpha_damping = " + std::to_string(alpha_damping) + "\n";
    source_str += "dt = " + std::to_string(dt) + "\n";
    source_str += "jit = " + std::to_string(jit) + "\n";
    source_str += "use_generated_code = " + std::to_string(use_generated_code) + "\n";
    source_str += "generate_code = " + std::to_string(generate_code) + "\n";
    source_str += "viz_all = " + std::to_string(viz_all) + "\n";

    // Modules
    source_str += "\n";
    source_str += "[modules]\n";
    source_str += "information_cost = " + std::to_string(information_cost) + "\n";
    source_str += "path_tracking = " + std::to_string(path_tracking) + "\n";
    source_str += "corridor_constraint = " + std::to_string(corridor_constraint) + "\n";
    source_str += "communication_constraint = " + std::to_string(communication_constraint) + "\n";
    source_str += "distribution_cost = " + std::to_string(distribution_cost) + "\n";
    source_str += "agent_collision_constraint = " + std::to_string(agent_collision_constraint) + "\n";
    source_str += "nearest_obstacle_constraint = " + std::to_string(nearest_obstacle_constraint) + "\n";

    // Information
    source_str += "\n";
    source_str += "[information]\n";
    source_str += "linear_weight = " + std::to_string(information.linear_weight) + "\n";

    // Path Tracking
    source_str += "\n";
    source_str += "[tracking]\n";
    source_str += "progress_weight = " + std::to_string(tracking.progress_weight) + "\n";
    source_str += "lag_weight = " + std::to_string(tracking.lag_weight) + "\n";
    source_str += "contouring_weight = " + std::to_string(tracking.contouring_weight) + "\n";
    source_str += "max_progress = " + std::to_string(tracking.max_progress) + "\n";

    // Nearest Obstacle
    source_str += "\n";
    source_str += "[nearest_obstacle]\n";
    source_str += "slack_weight = " + std::to_string(nearest_obstacle.slack_weight) + "\n";

    // Collision
    source_str += "\n";
    source_str += "[collision]\n";
    source_str += "slack_weight = " + std::to_string(collision.slack_weight) + "\n";
    source_str += "collision_distance = " + std::to_string(collision.collision_distance) + "\n";

    // Connectivity
    source_str += "\n";
    source_str += "[connectivity]\n";
    source_str += "slack_weight = " + std::to_string(connectivity.slack_weight) + "\n";
    source_str += "r_c = " + std::to_string(connectivity.r_c)+ "\n";
    source_str += "fiedler_eps = " + std::to_string(connectivity.fiedler_eps)+
            "\n";

    // Distribution
    source_str += "\n";
    source_str += "[distribution]\n";
    source_str += "lag_weight = " + std::to_string(distribution.lag_weight) +
            "\n";
    source_str += "distribution_weight = " + std::to_string(distribution.distribution_weight) +
            "\n";

    // Corridor
    source_str += "\n";
    source_str += "[corridor]\n";
    source_str += "slack_weight = " + std::to_string(corridor.slack_weight) + "\n";
    source_str += "lag_weight = " + std::to_string(corridor.lag_weight) + "\n";

    // Export
    std::ofstream out(path + "/configs/sqp_conf.toml");
    out << source_str;
    out.close();
}

void SQPParams::setNumAgents(int num) {
    num_agents = num;
    if(num == 1){
        communication_constraint = false;
        distribution_cost = false;
        agent_collision_constraint = false;
    }
    else{
        communication_constraint = true;
        distribution_cost = true;
        agent_collision_constraint = true;
    }
}

SQPParams::SQPParams(std::string path) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(path);
//        std::cout << "Loading SQP config: " << std::endl;
//        std::cout << tbl << "\n" << "\n";

        // Load general params
        N_horizon = read_toml_value<int>(tbl["general"]["N_horizon"]);
        num_agents = read_toml_value<int>(tbl["general"]["num_agents"]);

        agent_type = convert(read_toml_value<std::string>
                                     (tbl["general"]["agent_type"]));
        num_inputs = getNumInputs(agent_type);
        num_states = getNumStates(agent_type);

        sqp_iters = read_toml_value<int>(tbl["general"]["sqp_iters"]);
        alpha_damping = read_toml_value<double>(tbl["general"]["alpha_damping"]);
        dt = read_toml_value<double>(tbl["general"]["dt"]);
        jit = read_toml_value<bool>(tbl["general"]["jit"]);
        use_generated_code = read_toml_value<bool>(tbl["general"]["use_generated_code"]);
        generate_code = read_toml_value<bool>(tbl["general"]["generate_code"]);
        viz_all = read_toml_value<bool>(tbl["general"]["viz_all"]);
        information_cost = read_toml_value<bool>
                (tbl["modules"]["information_cost"]);
        if(information_cost){
            information.linear_weight =
                    read_toml_value<double>(tbl["information"]["linear_weight"]);
        }
        path_tracking = read_toml_value<bool>(tbl["modules"]["path_tracking"]);
        corridor_constraint = read_toml_value<bool>
                (tbl["modules"]["corridor_constraint"]);
        if(num_agents > 1) {
            communication_constraint = read_toml_value<bool>
                    (tbl["modules"]["communication_constraint"]);
            distribution_cost = read_toml_value<bool>
                    (tbl["modules"]["distribution_cost"]);
            agent_collision_constraint = read_toml_value<bool>
                    (tbl["modules"]["agent_collision_constraint"]);

        }
        else {
            communication_constraint = false;
            distribution_cost = false;
            agent_collision_constraint = false;
        }
        // Connectivity
        connectivity.slack_weight = read_toml_value<double>
                (tbl["connectivity"]["slack_weight"]);
        connectivity.r_c = read_toml_value<double>
                (tbl["connectivity"]["r_c"]);
        connectivity.fiedler_eps = read_toml_value<double>
                (tbl["connectivity"]["fiedler_eps"]);
        // Path Tracking
        tracking.progress_weight = read_toml_value<double>
                (tbl["tracking"]["progress_weight"]);
        tracking.lag_weight = read_toml_value<double>
                (tbl["tracking"]["lag_weight"]);
        tracking.contouring_weight = read_toml_value<double>
                (tbl["tracking"]["contouring_weight"]);
        tracking.max_progress = read_toml_value<double>
                (tbl["tracking"]["max_progress"]);

        // Distribution
        distribution.lag_weight = read_toml_value<double>
                (tbl["distribution"]["lag_weight"]);
        distribution.distribution_weight = read_toml_value<double>
                (tbl["distribution"]["distribution_weight"]);
        distribution.debug_draw = read_toml_value<bool>
                (tbl["distribution"]["debug_draw"]);

        // Collision
        collision.slack_weight = read_toml_value<double>
                (tbl["collision"]["slack_weight"]);
        collision.collision_distance = read_toml_value<double>
                (tbl["collision"]["collision_distance"]);
        // Corridor
        corridor.slack_weight = read_toml_value<double>
                (tbl["corridor"]["slack_weight"]);
        corridor.lag_weight = read_toml_value<double>
                (tbl["corridor"]["lag_weight"]);
        corridor.debug_draw = read_toml_value<bool>
                (tbl["corridor"]["debug_draw"]);

        // Nearest obstacle
        nearest_obstacle.slack_weight = read_toml_value<double>
                (tbl["nearest_obstacle"]["slack_weight"]);

    }
    catch (const toml::parse_error& err)
    {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
    for(int i = 0; i < num_agents; i++){
        colors.push_back(new sf::Color(rand() % 255, rand() % 255, rand() %
                                                                   255));
    }
}

BicycleParams::BicycleParams() {
    toml::table tbl;
    try {
        tbl = toml::parse_file(RES_LOC "/configs/bicycle_conf.toml");
//        std::cout << "Loading Dynamics config: " << std::endl;
//        std::cout << tbl << "\n" << "\n";

        // Load general params
        yaw_rate_max = read_toml_value<float>(tbl["yaw_rate_max"]);
        accel_max = read_toml_value<float>(tbl["accel_max"]);
        speed_max = read_toml_value<float>(tbl["speed_max"]);
    }
    catch (const toml::parse_error& err)
    {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}

void BicycleParams::exportConfig(std::string path) {
    std::string source_str = "";
    source_str += "[general]\n";
    source_str += "yaw_rate_max = " + std::to_string(yaw_rate_max) + "\n";
    source_str += "accel_max = " + std::to_string(accel_max) + "\n";
    source_str += "speed_max = " + std::to_string(speed_max) + "\n";

    // Export
    std::ofstream out(path + "/configs/bicycle_conf.toml");
    out << source_str;
    out.close();
}

Dynamics *createAgent(SQPParams *params, int id) {
        switch(params->agent_type){
            case bicycle:
                return new BicycleDynamics(params, id);
            case hippocampus:
                return new HippocampusDynamics(params, id);
            case bluerov:
                return new BlueROVDynamics(params, id);
            default:
                throw std::runtime_error("Invalid agent type");
        }
}

int getNumStates(AgentType agent_type) {
    switch(agent_type){
        case bicycle:
            return BicycleDynamics::getNumStates();
        case hippocampus:
            return HippocampusDynamics::getNumStates();
        case bluerov:
            return BlueROVDynamics::getNumStates();
        default:
            throw std::runtime_error("Invalid agent type");
    }}

int getNumInputs(AgentType agent_type) {
    switch(agent_type){
        case bicycle:
            return BicycleDynamics::getNumInputs();
        case hippocampus:
            return HippocampusDynamics::getNumInputs();
        case bluerov:
            return BlueROVDynamics::getNumInputs();
        default:
            throw std::runtime_error("Invalid agent type");
    }
}

AgentType convert(const std::string &str) {
    if(str == "bicycle"){
        return bicycle;
    }
    else if(str == "hippocampus"){
        return hippocampus;
    }
    else if(str == "bluerov"){
        return bluerov;
    }
    else{
        throw std::runtime_error("Invalid agent type \"" + str + "\"");
    }
}

std::string convert(const AgentType &agent_type) {
    switch(agent_type){
        case bicycle:
            return "bicycle";
        case hippocampus:
            return "hippocampus";
        case bluerov:
            return "bluerov";
        default:
            throw std::runtime_error("Invalid agent type");
    }
}

HippocampusParams::HippocampusParams() {
    toml::table tbl;
    try {
        tbl = toml::parse_file(RES_LOC "/configs/hippocampus_conf.toml");
//        std::cout << "Loading Hippocampus config: " << std::endl;
//        std::cout << tbl << "\n" << "\n";

        // Load general params
        yaw_rate_max = read_toml_value<double>(tbl["yaw_rate_max"]);
        accel_max = read_toml_value<double>(tbl["accel_max"]);
        speed_max = read_toml_value<double>(tbl["speed_max"]);
        damping_factor = read_toml_value<double>(tbl["damping_factor"]);
    }
    catch (const toml::parse_error &err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}

void HippocampusParams::exportConfig(std::string path) {
    std::string source_str = "";
    source_str += "[general]\n";
    source_str += "yaw_rate_max = " + std::to_string(yaw_rate_max) + "\n";
    source_str += "accel_max = " + std::to_string(accel_max) + "\n";
    source_str += "speed_max = " + std::to_string(speed_max) + "\n";

    // Export
    std::ofstream out(path + "/configs/hippocampus_conf.toml");
    out << source_str;
    out.close();
}


BlueROVParams::BlueROVParams() {
    toml::table tbl;
    try {
        tbl = toml::parse_file(RES_LOC "/configs/bluerov_conf.toml");
//        std::cout << "Loading BlueROV config: " << std::endl;
//        std::cout << tbl << "\n" << "\n";

        // Load general params
        yaw_rate_max = read_toml_value<double>(tbl["yaw_rate_max"]);
        accel_max = read_toml_value<double>(tbl["accel_max"]);
        speed_max = read_toml_value<double>(tbl["speed_max"]);
        damping_factor = read_toml_value<double>(tbl["damping_factor"]);
    }
    catch (const toml::parse_error &err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
}

void BlueROVParams::exportConfig(std::string path) {
    std::string source_str = "";
    source_str += "[general]\n";
    source_str += "yaw_rate_max = " + std::to_string(yaw_rate_max) + "\n";
    source_str += "accel_max = " + std::to_string(accel_max) + "\n";
    source_str += "speed_max = " + std::to_string(speed_max) + "\n";
    source_str += "damping_factor = " + std::to_string(damping_factor) + "\n";

    // Export
    std::ofstream out(path + "/configs/bluerov_conf.toml");
    out << source_str;
    out.close();
}
