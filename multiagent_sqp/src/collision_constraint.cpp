//
// Created by finn on 5/28/23.
//
#include "collision_constraint.h"

// TODO Pass params struct
CollisionConstraint::CollisionConstraint(int time_id, SQPParams* params):
time_id(time_id), params(params) {

}

void CollisionConstraint::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                        MX &delta_z, MX &delta_u,
                                        OptiWrapper *opti_obj, DM &states,
                                        MX &ca_param_vec, DM &param_vec) {
    auto slack = opti_obj->ca_variables(Slice(sl_slack[0], sl_slack[1]));

    // For agent 0: N - 1 constraints
    // For agent 1: N - 2 ...
    int counter = 0;
    for(int i = 0; i < params->num_agents - 1; i++){
        auto x0 = ca_states(i * params->num_states, time_id);
        auto dx0 = delta_z(i * params->num_states, time_id);
        auto y0 = ca_states(i * params->num_states + 1, time_id);
        auto dy0 = delta_z(i * params->num_states + 1, time_id);
        for(int j = i + 1; j < params->num_agents; j++){
            auto x1 = ca_states(j * params->num_states, time_id);
            auto dx1 = delta_z(j * params->num_states, time_id);
            auto y1 = ca_states(j * params->num_states + 1, time_id);
            auto dy1 = delta_z(j * params->num_states + 1, time_id);
            MX grad_term = dx0 * ( 2 * x0 - 2 * x1)
                    + dx1*(-2*x0 + 2*x1) + dy0*(2*y0 - 2*y1) + dy1*(-2*y0 + 2*y1);
            MX func_eval = pow(x1 - x0, 2) + pow(y1 - y0, 2);
            opti_obj->subject_to(func_eval + grad_term >= params->collision.collision_distance
            - pow
            (slack
            (counter), 2));
            cost_fun += slack(counter) * params->collision.slack_weight;
            opti_obj->subject_to(slack >= 0);
            counter += 1;
        }
    }
}

//bool CollisionConstraint::checkFeasible(Opti *opti_obj, DM &states) {
//    auto slack_val = opti_obj->value(slack);
//    int counter = 0;
//    for(int i = 0; i < params->num_agents - 1; i++){
//        auto x0 = states(i * params->num_states, time_id);
//        auto y0 = states(i * params->num_states + 1, time_id);
//        for(int j = i + 1; j < params->num_agents; j++){
//
//            auto x1 = states(j * params->num_states, time_id);
//            auto y1 = states(j * params->num_states + 1, time_id);
//            double func_eval = (pow(x1 - x0, 2) + pow(y1 - y0, 2)).scalar();
//            if(func_eval + 0.01 < 0.16 - pow(slack_val(counter).scalar(), 2)){
//                return false;
//            }
//            counter += 1;
//        }
//    }
//    return true;
//}

int CollisionConstraint::addToParams(int start_idx) noexcept {
    return 0;
}

int CollisionConstraint::addToVariables(int start_idx) noexcept {
    sl_slack[0] = start_idx;
    sl_slack[1] = start_idx + params->num_agents * (params->num_agents - 1);

    return params->num_agents * (params->num_agents - 1);
}
