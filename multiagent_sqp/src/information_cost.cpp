//
// Created by finn on 5/29/23.
//
#include "information_cost.h"

InformationCost::InformationCost(SQPParams* params,  int id, InfoField* field)
        : params(params), id(id), field(field) {

}

void InformationCost::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                    MX &delta_z, MX &delta_u,
                                    OptiWrapper *opti_obj, DM &states,
                                    MX &ca_param_vec, DM &param_vec) {

    auto gradientx = ca_param_vec(Slice(sl_gradient_x[0], sl_gradient_x[1]));
    auto gradienty = ca_param_vec(Slice(sl_gradient_y[0], sl_gradient_y[1]));
    auto value = ca_param_vec(Slice(sl_value[0], sl_value[1]));
    // adding the cost terms:
    Slice state_slice(id*params->num_states, (id+1)*params->num_states);
    Slice all;
    auto dz = delta_z(state_slice, all);
    for(int i = 0; i < params->N_horizon; i++) {
        MX c1 = gradientx(i) * dz(0, i) + gradienty(i) * dz(1, i);
        cost_fun += c1 * params->information.linear_weight;
    }
    linearize(opti_obj, states, param_vec);
}

void InformationCost::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    Slice state_slice(id*params->num_states, (id+1)*params->num_states);
    Slice all;
    auto states_sl = states(state_slice, all);
    for(int i = 0; i < params->N_horizon; i++){
        auto grads = field->getGradient(states_sl(0, i).scalar(), states_sl(1, i).scalar());
//        std::cout << "Dynamics " << id << ", x: " << states_sl(0, i).scalar() << ", y: " << states_sl(1, i).scalar() << " gradient: " << grads[0] << ", " << grads[1] << std::endl;
        param_vec(sl_gradient_x[0] + i) = -grads[0];
        param_vec(sl_gradient_y[0] + i) = -grads[1];
    }
}

int InformationCost::addToParams(int start_idx) noexcept {
    // gradient x and y
    sl_gradient_x = {start_idx, start_idx + params->N_horizon};
    sl_gradient_y = {start_idx + params->N_horizon,
                     start_idx + params->N_horizon * 2};
    sl_value = {start_idx + params->N_horizon * 2,
                start_idx + params->N_horizon * 3};
    return params->N_horizon * 3;
}

int InformationCost::addToVariables(int start_idx) noexcept {
    return 0;
}
