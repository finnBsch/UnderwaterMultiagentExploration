//
// Created by finn on 5/28/23.
//

#ifndef SQP_COLLISION_CONSTRAINT_H
#define SQP_COLLISION_CONSTRAINT_H
#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include "bicycle_dynamics.h"
#include "sqp_module.h"
#include "configs.h"

using namespace casadi;
namespace E = Eigen;


/**
 * Implementation of collision constraint for circles for one timestep.
 * For now, just first order gradients, cannot be used with full SQP approach.
 * Will create 'n * (n-1)/2' constraints for all agents pairwise.
 * Given a collision radius r, we can derive for two agents:
 * (x1 - x0)^2 + (y1 - y0)^2 >= r^2
 *
 */
class CollisionConstraint: public SQPModule{
private:
    int time_id;
    SQPParams* params;

    // Variables
    std::array<int, 2> sl_slack;
public:
    int addToParams(int start_idx) noexcept;
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u,
                       OptiWrapper *opti_obj, DM &states, MX &ca_param_vec,
                       DM &param_vec);
    void iterationUpdate(OptiWrapper *opti_obj) {}
    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {}
    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {}
    void draw(sf::RenderTarget& target, sf::RenderStates states) const{}
    CollisionConstraint(int time_id, SQPParams* params);
//    bool checkFeasible(Opti *opti_obj, DM &states);

    int addToVariables(int start_idx) noexcept override;

};
#endif //SQP_COLLISION_CONSTRAINT_H
