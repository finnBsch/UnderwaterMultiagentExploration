//
// Created by finn on 6/25/23.
//

#ifndef MASTERMAIN_NEAREST_OBSTACLE_AVOIDANCE_H
#define MASTERMAIN_NEAREST_OBSTACLE_AVOIDANCE_H
#include "sqp_module.h"
#include "configs.h"
#include "path.h"
class NearestObstacleAvoidance : public SQPModule{
private:
    SQPParams* params;
    int id;
    std::vector<SQPObstacle*>* obstacles;

    Function grad_fn;
    Function dist_fn;
    Function grad_to_path_fn;

    // Parameters
    std::array<int, 2> sl_obs_x0;
    std::array<int, 2> sl_obs_y0;
    std::array<int, 2> sl_obs_x1;
    std::array<int, 2> sl_obs_y1;
    std::array<int, 2> sl_dist;
    std::array<int, 2> sl_dir_t;


    // Variables
    std::array<int, 2> sl_slack;

public:
    NearestObstacleAvoidance(SQPParams* params, int id, std::vector<SQPObstacle*>* obstacles);
    NearestObstacleAvoidance(SQPParams* params, int id,
                             std::vector<SQPObstacle*>* obstacles, sf::Color*
                             color);

    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u, OptiWrapper *opti_obj, DM &states,
                       MX &ca_param_vec, DM &param_vec) override;

    int addToParams(int start_idx) noexcept override;

    int addToVariables(int start_idx) noexcept override;

    void iterationUpdate(OptiWrapper *opti_obj) override;

    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

    void buildGradients();
};

#endif //MASTERMAIN_NEAREST_OBSTACLE_AVOIDANCE_H
