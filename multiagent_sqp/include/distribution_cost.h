//
// Created by finn on 6/19/23.
//

#ifndef MASTERMAIN_DISTRIBUTION_COST_H
#define MASTERMAIN_DISTRIBUTION_COST_H
#include "sqp_module.h"
#include "configs.h"
#include "path.h"

class DistributionCost : public SQPModule {
private:
    std::vector<sf::CircleShape*> circles;
    SQPParams* params;
    SplinePath* path;
    // Slices
    std::array<int, 2> sl_theta;
    std::array<int, 2> sl_theta0;
    std::array<int, 2> sl_x0;
    std::array<int, 2> sl_y0;
    std::array<int, 2> sl_alpha;
    std::array<int, 2> sl_xbar;


    // Variables
    std::array<int, 2> sl_delta_theta;
    std::array<int, 2> sl_v;

    Function Jac_Cost;
    Function Jac_Contouring;
    Function Hess_Cost;
    Function Hess_Contouring;

    DM current_theta;
//    MX delta_theta;



public:
    DistributionCost(SQPParams* params, SplinePath* path);

    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u, OptiWrapper *opti_obj, DM &states,
                       MX &ca_param_vec, DM &param_vec) override;

    int addToParams(int start_idx) noexcept override;

    void iterationUpdate(OptiWrapper *opti_obj) override;

    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void draw(sf::RenderTarget &target, sf::RenderStates states) const
    override;

    void buildCasadiGradients();

    int addToVariables(int start_idx) noexcept override;
};

#endif //MASTERMAIN_DISTRIBUTION_COST_H
