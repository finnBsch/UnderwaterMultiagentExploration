//
// Created by finn on 5/29/23.
//

#ifndef SQP_PATH_TRACKING_H
#define SQP_PATH_TRACKING_H
#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include "sqp_module.h"
#include "configs.h"
#include "path.h"
using namespace casadi;


/**
 * Implements the path-tracking cost for one agent for piecewise linear
 * paths.
 * See https://github.com/alexliniger/MPCC and https://arxiv.org/pdf/1711.07300.pdf
 */
class PathTracking : public SQPModule{
private:
    Path path;
    SQPParams* params;
    int id;
    MX theta;  // Approximated progress variable
    MX theta0;  // Approximated progress variable
    MX delta_theta;  // Approximated progress variable
    DM current_theta;
    MX theta_0_constr;
    MX x0;
    MX y0;
    MX alpha;
    MX v;  // progress velocity
    // coordinates of *current* theta state. Function of linearization point
    MX xbar;
    std::vector<sf::CircleShape> lin_pts;

public:
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u,
                       OptiWrapper *opti_obj, DM &states, MX &ca_param_vec,
                       DM &param_vec) override;
    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;
    void iterationUpdate(OptiWrapper *opti_obj) override;
    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;
    PathTracking(SQPParams* params, int id);
    bool checkFeasible(Opti *opti_obj, DM &states){
        return true;
    }
    void draw(sf::RenderTarget& target, sf::RenderStates states) const;
};
#endif //SQP_PATH_TRACKING_H
