//
// Created by finn on 5/28/23.
//

#ifndef SQP_CONFIGS_H
#define SQP_CONFIGS_H
#include <cmath>
#include "res_loc.h"
#include "string"
#include <stdexcept>
#include <SFML/Graphics.hpp>

enum AgentType {
    bicycle = 0,
    hippocampus = 1,
    bluerov = 2
};

struct SQPParams;
class Dynamics;

Dynamics* createAgent(SQPParams* params, int id);
int getNumStates(AgentType agent_type);
int getNumInputs(AgentType agent_type);


AgentType convert(const std::string& str);
std::string convert(const AgentType& agent_type);

struct InformationParams{
    double linear_weight;
};

struct CollisionParams{
    double slack_weight;
    double collision_distance;
};

struct ConnectivityParams{
    double slack_weight;
    double r_c;  // Communication radius
    double fiedler_eps;  // Epsilon for fiedler eigenvalue
};

struct TrackingParams{
    double progress_weight;
    double lag_weight;
    double contouring_weight;
    double max_progress;
};

struct NearestObstacleParams{
    double slack_weight;
};

struct CorridorParams{
    double slack_weight;
    double lag_weight;
    bool debug_draw;
};

struct DistributionParams{
    double lag_weight;
    double distribution_weight;
    bool debug_draw;
};

struct SQPParams{
    double size_x;
    double size_y;
    int N_horizon;
    AgentType agent_type;
    int num_agents;
    int num_inputs;
    int num_states;
    double alpha_damping; // 0 equals no damping
    double dt;
    bool jit;
    bool use_generated_code;
    bool generate_code;

    bool viz_all;
    bool information_cost;
    bool path_tracking;
    bool corridor_constraint;
    bool nearest_obstacle_constraint;
    bool communication_constraint;
    bool distribution_cost;
    bool agent_collision_constraint;

    int sqp_iters;
    InformationParams information;
    TrackingParams tracking;
    CorridorParams corridor;
    CollisionParams collision;
    ConnectivityParams connectivity;
    DistributionParams distribution;
    NearestObstacleParams nearest_obstacle;

    std::vector<sf::Color*> colors;

    SQPParams();
    SQPParams(std::string path);

    void setNumAgents(int num);
    void exportConfig(std::string path);
};


struct AgentParams {
    double speed_max;
};

struct BicycleParams: AgentParams{
    double yaw_rate_max = 1.5;
    double accel_max = 1.0;
    double speed_max = 2.0;

    BicycleParams();
    void exportConfig(std::string path);
};

struct HippocampusParams : AgentParams{
    double yaw_rate_max;
    double accel_max;
    double speed_max;
    double damping_factor;

    HippocampusParams();
    void exportConfig(std::string path);
};

struct BlueROVParams : AgentParams{
    double yaw_rate_max;
    double accel_max;
    double speed_max;
    double damping_factor;

    BlueROVParams();
    void exportConfig(std::string path);
};
#endif //SQP_CONFIGS_H
