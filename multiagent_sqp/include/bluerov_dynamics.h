//
// Created by finn on 6/19/23.
//

#ifndef MASTERMAIN_BLUEROV_DYNAMICS_H
#define MASTERMAIN_BLUEROV_DYNAMICS_H

#include "dynamics.h"
using namespace casadi;

class BlueROVDynamics : public Dynamics {
private:
    static const int num_states = 5;
    static const int num_inputs = 3;
    SQPParams* params;
    Function dynamics_grad;
    Function dynamics;
public:
    BlueROVParams a_params;
    BlueROVDynamics(SQPParams *params, int id);
    static int getNumStates(){return num_states;};
    static int getNumInputs(){return num_inputs;};

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


#endif //MASTERMAIN_BLUEROV_DYNAMICS_H
