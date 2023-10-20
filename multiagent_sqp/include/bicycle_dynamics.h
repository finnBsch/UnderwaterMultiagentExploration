//
// Created by finn on 5/18/23.
//

#ifndef CASADI_TEST_SIMPLE_AGENT_H
#define CASADI_TEST_SIMPLE_AGENT_H

#include "dynamics.h"

using namespace casadi;
using e_state = Eigen::Matrix<double, 4, 1>;
using e_u = Eigen::Matrix<double, 2, 1>;



class BicycleDynamics : public Dynamics {
private:
    std::array<double, 4> init_state = {0, 0, 0, 0};
    const int num_states = 4;  // x, y, phi, v
    const int num_inputs = 2;
    SQPParams* params;
    MX mu; // Lagrangian for dynamics.
    DM mu_vals;
    MX dmu; // Lagrangian for dynamics.
    MX x_init; // Initial state
public:
    BicycleParams a_params;
    MX getDeriv(const MX& state, const MX& u);
    e_state getDeriv(const e_state& state, const e_u& u);
    BicycleDynamics(SQPParams* params, int id);
    // SQP Builder
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u,
                       OptiWrapper *opti_obj, DM &states, MX &ca_param_vec, DM &param_vec) override;
    void iterationUpdate(OptiWrapper *opti_obj) override;
    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;
    void draw(sf::RenderTarget& target, sf::RenderStates states) const override {};
    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override{};
    void initializeNextSol(Opti* opti_obj, Eigen::Matrix<double, Eigen::Dynamic,
                           1> init_state);

    int addToParams(int start_idx) noexcept override{return 0;};
    static int getNumStates(){return 4;};
    static int getNumInputs(){return 2;}

    int addToVariables(int start_idx) noexcept override;

};



#endif //CASADI_TEST_SIMPLE_AGENT_H
