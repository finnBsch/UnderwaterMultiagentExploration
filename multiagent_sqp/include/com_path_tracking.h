//
// Created by finn on 5/31/23.
//

#ifndef SQP_COM_PATH_TRACKING_H
#define SQP_COM_PATH_TRACKING_H
#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include "sqp_module.h"
#include "configs.h"
#include "path.h"
using namespace casadi;


class COMPathTracking : public SQPModule{
private:
    SplinePath* path;
    SQPParams* params;

    std::array<int, 2> sl_x0;
    std::array<int, 2> sl_y0;
    std::array<int, 2> sl_theta0;
    std::array<int, 2> sl_alpha;
    std::array<int, 2> sl_theta;
    std::array<int, 2> sl_xbar;
    std::array<int, 2> sl_theta0_constr;

    // Variables
    std::array<int, 2> sl_v;
    std::array<int, 2> sl_delta_theta;
//    MX theta;  // Approximated progress variable
//    MX theta0;  // Approximated progress variable
    DM current_theta;
//    MX theta_0_constr;
//    MX x0;
//    MX y0;
//    MX alpha;
    // coordinates of *current* theta state. Function of linearization point
//    MX xbar;
    std::vector<sf::CircleShape> lin_pts;
    sf::CircleShape ct_pt;
    Function Hess;
    Function Jac;

public:
    double getTheta() const;
    int addToParams(int start_idx) noexcept;
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u,
                       OptiWrapper *opti_obj, DM &states, MX &ca_param_vec,
                       DM &param_vec) override;
    void deriveGradients();
    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;
    void iterationUpdate(OptiWrapper *opti_obj) override;
    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;
    COMPathTracking(SQPParams* params, SplinePath* path);
    bool checkFeasible(Opti *opti_obj, DM &states){
        return true;
    }
    void draw(sf::RenderTarget& target, sf::RenderStates states) const;

    int addToVariables(int start_idx) noexcept override;

};
#endif //SQP_COM_PATH_TRACKING_H
