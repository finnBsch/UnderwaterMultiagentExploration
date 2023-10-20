//
// Created by finn on 6/5/23.
//

#ifndef SQP_CORRIDOR_CONSTRAINT_H
#define SQP_CORRIDOR_CONSTRAINT_H
#include "sqp_module.h"
#include "configs.h"
#include "path.h"

using namespace casadi;

class CorridorConstraint : public SQPModule{
private:
    std::vector<sf::CircleShape*> circles_0;
    std::vector<sf::CircleShape*> circles_1;
    std::vector<sf::CircleShape*> circles_c;

    std::array<int, 2> sl_theta;
    std::array<int, 2> sl_theta0;
    std::array<int, 2> sl_x0;
    std::array<int, 2> sl_y0;
    std::array<int, 2> sl_alpha;
    std::array<int, 2> sl_xbar;
    std::array<int, 2> sl_fl0;
    std::array<int, 2> sl_fr0;
    std::array<int, 2> sl_dalpha_l;
    std::array<int, 2> sl_dalpha_r;
    std::array<int, 2> sl_maxtheta;

    // Variables
    std::array<int, 2> sl_slack;
    std::array<int, 2> sl_delta_theta;
    std::array<int, 2> sl_v;

    SQPParams* params;
    DM current_theta;
    MX max_theta;
    SplinePath* path;
    Function Jac_Constr_L;
    Function Jac_Constr_R;
    Function Hess;
    Function Jac_Cost;
    // DEBUG
    Function e_cost;
public:
    CorridorConstraint(SQPParams* params, SplinePath* path);
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u, OptiWrapper *opti_obj, DM &states,
                       MX &ca_param_vec, DM &param_vec) override;

    int addToParams(int start_idx) noexcept override;

    void iterationUpdate(OptiWrapper *opti_obj) override;

    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

    void buildGradients();

    double getTheta() const;

    int addToVariables(int start_idx) noexcept override;
};
#endif //SQP_CORRIDOR_CONSTRAINT_H
