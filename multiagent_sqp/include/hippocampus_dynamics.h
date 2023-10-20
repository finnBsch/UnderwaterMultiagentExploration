//
// Created by finn on 6/16/23.
//

#ifndef MASTERMAIN_HIPPOCAMPUS_DYNAMICS_H
#define MASTERMAIN_HIPPOCAMPUS_DYNAMICS_H
#include "dynamics.h"

class HippocampusDynamics : public Dynamics {
private:
    static const int num_states = 5; // x, y, phi, v_x, v_y
    static const int num_inputs = 2; // d_phi, a
    SQPParams* params;
    Function dynamics_grad; // Gradient of dynamics
    Function dynamics; // Dynamics
public:
    HippocampusParams a_params; // Parameters for dynamics
    HippocampusDynamics(SQPParams *params, int id);
    static int getNumStates(){return num_states;};
    static int getNumInputs(){return num_inputs;};

    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u, OptiWrapper *opti_obj, DM &states,
                       MX &ca_param_vec, DM &param_vec) override;

    int addToParams(int start_idx) noexcept override;

    void iterationUpdate(OptiWrapper *opti_obj) override;

    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

    void buildCasadiGradients();

    int addToVariables(int start_idx) noexcept override;

};


#endif //MASTERMAIN_HIPPOCAMPUS_DYNAMICS_H
