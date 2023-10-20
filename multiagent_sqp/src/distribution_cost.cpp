//
// Created by finn on 6/19/23.
//
#include "distribution_cost.h"

void DistributionCost::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                     MX &delta_z, MX &delta_u, OptiWrapper *opti_obj,
                                     DM &states, MX &ca_param_vec,
                                     DM &param_vec) {
    buildCasadiGradients();
    auto x0 = reshape(ca_param_vec(Slice(sl_x0[0], sl_x0[1])),
                      params->num_agents, params->N_horizon);
    auto y0 = reshape(ca_param_vec(Slice(sl_y0[0], sl_y0[1])),
                        params->num_agents, params->N_horizon);
    auto alpha = reshape(ca_param_vec(Slice(sl_alpha[0], sl_alpha[1])),
                         params->num_agents, params->N_horizon);
    auto xbar = reshape(ca_param_vec(Slice(sl_xbar[0], sl_xbar[1])),
                        params->num_agents, params->N_horizon);
    auto theta = reshape(ca_param_vec(Slice(sl_theta[0], sl_theta[1])),
                         params->num_agents, params->N_horizon);
    auto theta0 = reshape(ca_param_vec(Slice(sl_theta0[0], sl_theta0[1])),
                          params->num_agents, params->N_horizon);

    current_theta = DM::zeros(params->num_agents * params->N_horizon, 1);
    current_theta += 0.2;
    auto delta_theta = reshape(opti_obj->ca_variables(Slice(sl_delta_theta[0],
                                          sl_delta_theta[1])),
                               params->num_agents, params->N_horizon);
    auto v = reshape(opti_obj->ca_variables(Slice(sl_v[0], sl_v[1])),
                     params->num_agents, params->N_horizon);
    Slice all;
//    opti_obj->subject_to(delta_theta(all, 0) >= 0);
    for(int t = 1; t < params->N_horizon; t++) {
        auto z = ca_states(all, t);
        auto dz = delta_z(all, t);

        auto decision_vars = vertcat(z, theta(all, t));
        auto delta_vars = vertcat(dz, delta_theta(all, t));
//    auto delta_vars = dz;
        // Inputs: states, theta_, theta_0, alpha_, x0_,
        //                                     y0_
        auto J_map = Jac_Cost({{"i0", z},
                               {"i1", theta(all, t)},
                               {"i2", theta0(all, t)},
                               {"i3", alpha(all, t)},
                               {"i4", x0(all, t)},
                               {"i5", y0(all, t)}});
        auto J_cont_map = Jac_Contouring({{"i0", z},
                                    {"i1", theta(all, t)},
                                    {"i2", theta0(all, t)},
                                    {"i3", alpha(all, t)},
                                    {"i4", x0(all, t)},
                                    {"i5", y0(all, t)}});
        auto J = J_map["o0"];
        auto J_cont = J_cont_map["o0"];

        auto H_map = Hess_Cost({{"i0", z},
                                {"i1", theta(all, t)},
                                {"i2", theta0(all, t)},
                                {"i3", alpha(all, t)},
                                {"i4", x0(all, t)},
                                {"i5", y0(all, t)}});
        auto H_cont_map = Hess_Contouring({{"i0", z},
                                           {"i1", theta(all, t)},
                                           {"i2", theta0(all, t)},
                                           {"i3", alpha(all, t)},
                                           {"i4", x0(all, t)},
                                           {"i5", y0(all, t)}});
        auto H = H_map["o0"];
        auto H_cont = H_cont_map["o0"];

        cost_fun += mtimes(J, delta_vars)
                    + mtimes(transpose(delta_vars),
                             mtimes(H,delta_vars));
        cost_fun += mtimes(J_cont, delta_theta(all, t))
                    + mtimes(transpose(delta_theta(all, t)),
                             mtimes(H_cont,delta_theta(all, t)));
        opti_obj->subject_to(delta_theta(all, t) + theta(all, t) ==
        delta_theta(all, t - 1) + theta(all, t - 1) + v(all, t - 1) *
                                                          params->dt);
        for(int i = 0; i < params->num_agents; i++){
            cost_fun += v(i, t - 1) * 50.0;
        }
        opti_obj->subject_to(0.2 < v(all, t - 1));

    }
    opti_obj->subject_to(delta_theta(all, 0) == 0);

}


void DistributionCost::iterationUpdate(OptiWrapper *opti_obj) {
    DM dthet = opti_obj->value(sl_delta_theta[0], sl_delta_theta[1]);
    current_theta = current_theta + params->alpha_damping*dthet;
    opti_obj->update_parameter(sl_theta[0], sl_theta[1], current_theta);
}

void DistributionCost::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    for (int i = 0; i < params->num_agents; i++) {
        for (int t = 0; t < params->N_horizon; t++) {
            int index = t * params->num_agents + i;
            auto seg_pts = path->linSegment(current_theta(index).scalar());
            if(params->distribution.debug_draw) {
                if (t == 0) {
                    circles[i]->setPosition(seg_pts[0], seg_pts[1]);
                }
            }
            param_vec(sl_xbar[0] + index, 0) = cos(seg_pts[3]) *
                    (current_theta(index)
                                                              - seg_pts[2]) +
                                           seg_pts[0];
            param_vec(sl_theta0[0] + index, 0) = seg_pts[2];
            param_vec(sl_x0[0] + index, 0) = seg_pts[0];
            param_vec(sl_y0[0] + index, 0) = seg_pts[1];
            param_vec(sl_alpha[0] + index, 0) = seg_pts[3];
//            std::cout << "Dist" << std::endl;
//            std::cout << "t : " << t << ", fl0: " << seg_pts[4] << ", fr0: " << seg_pts[5] <<
//                      ", theta0 " << seg_pts[2]<< std::endl;
        }
    }
    param_vec(Slice(sl_theta[0], sl_theta[1]), 0) = current_theta;
//    std::cout << "theta end: " << reshape(current_theta, params->num_agents,
//                                          params->N_horizon)(Slice(), -1) <<
//                                          std::endl;

}

void DistributionCost::warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    Slice first(0, params->N_horizon - 1);
    Slice second(1, params->N_horizon);

    Slice all;
    auto theta = reshape(current_theta,
                         params->num_agents, params->N_horizon);
    theta(all, first) = theta(all, second);
    theta(all, -1) = theta(all, -2) * 2 - theta(all, -3);
    for(int i  = 0; i < params->N_horizon; i++) {
        auto mean = sum1(theta(all, i)) / params->num_agents;
        theta(all, i) = theta(all, i) + (mean - theta(all, i)) * 0.2;
    }
    current_theta = reshape(theta, params->num_agents * params->N_horizon, 1);


}

void DistributionCost::draw(sf::RenderTarget &target,
                            sf::RenderStates states) const {
    if(params->distribution.debug_draw) {
        for (auto &circ: circles) {
            target.draw(*circ, states);
        }
    }
}

void DistributionCost::buildCasadiGradients() {
    MX states = MX::sym("states", params->num_agents * params->num_states);
    MX e_c = MX::zeros(params->num_agents, 1);
    MX e_l = MX::zeros(params->num_agents, 1);
    MX theta_ = MX::sym("theta_", params->num_agents, 1);
    MX theta_0 = MX::sym("theta_0", params->num_agents, 1);
    MX alpha_ = MX::sym("alpha_", params->num_agents, 1);
    MX x0_ = MX::sym("x0_", params->num_agents, 1);
    MX y0_ = MX::sym("y0_", params->num_agents, 1);

//    auto decision_vars = states;
    MX cost = 0;
    MX cost_2 = 0;
    auto decision_vars = vertcat(states, theta_);
    for(int i = 0; i < params->num_agents; i++){
        auto x = states(i * params->num_states);
        auto y = states(i * params->num_states + 1);
        MX x_ = cos(alpha_(i)) * (theta_(i) - theta_0(i)) + x0_(i);
        e_l(i) = -(x - x_) * cos(alpha_(i)) -
                (y - (y0_(i) + sin(alpha_(i))
                *(theta_ (i) - theta_0(i)))) * sin(alpha_(i));
        e_c(i) = (x - x_) * sin(alpha_(i)) - (y - (y0_(i) + sin(alpha_(i)) *
                (theta_(i) - theta_0(i)))) * cos(alpha_(i));

    }
    // Get the min and max using mellowmax
//    double alpha_min = -100.0;
//    double alpha_max = 100;
//    auto exp_e_c_max = exp(e_c * alpha_max);
//    auto exp_e_c_min = exp(e_c * alpha_min);
//    auto min_ec = 1.0/alpha_min * log(1.0 / params->num_agents * sum1
//            (exp_e_c_min));
//    auto max_ec = 1.0/alpha_max * log(1.0 / params->num_agents * sum1
//            (exp_e_c_max));
    // Get the min and max using p-Norm
//    auto min_ec = pow(sum1(pow(-e_c, params->num_agents)), 1.0/params->num_agents);
//    auto max_ec = pow(sum1(pow(e_c, params->num_agents)), 1.0/params->num_agents);
    // Using softmax
//    double alpha_min = -3.0;
//    double alpha_max = 3.0;
//    auto exp_e_c_max = exp(e_c * alpha_max);
//    auto exp_e_c_min = exp(e_c * alpha_min);
//    auto max_ec = sum1(e_c * exp_e_c_max) / sum1(exp_e_c_max);
//    auto min_ec = sum1(e_c * exp_e_c_min) / sum1(exp_e_c_min);
//    auto cost = -(max_ec - min_ec) * 1.0;
    MX cont_cost = 0;
    for(int i = 0; i < params->num_agents - 1; i++){
        for(int j = i + 1; j < params->num_agents; j++) {
            cost -= pow(e_c(i) - e_c(j), 2) * params->distribution.distribution_weight;
        }
    }
    for(int i = 0; i < params->num_agents; i++){
//        cost += params->distribution.lag_weight * pow(e_l(i), 2);
        cont_cost += params->distribution.lag_weight * pow(e_l(i), 2);
        cont_cost += params->distribution.lag_weight/10 * pow(e_c(i), 2);
    }
    Jac_Cost = Function("Jac_Cost", {states, theta_, theta_0, alpha_, x0_,
                                     y0_}, {jacobian(cost, decision_vars)}, Dict({{"compiler", "shell"},
                                                                                  {"jit", params->jit},
                                                                                  {"jit_options"
                                                                                   ".compiler", "gcc"},
                                                                                  {"jit_options"
                                                                                   ".flags", "-O3"}}));
    Hess_Cost = Function("Hess_Cost", {states, theta_, theta_0, alpha_, x0_,
                                       y0_}, {hessian(cost_2, decision_vars)}, Dict({{"compiler", "shell"},
                                                                                     {"jit", params->jit},
                                                                                     {"jit_options"
                                                                                      ".compiler", "gcc"},
                                                                                     {"jit_options"
                                                                                      ".flags", "-O3"}}));
    Jac_Contouring = Function("Jac_Contouring", {states, theta_, theta_0, alpha_, x0_,
                                     y0_}, {jacobian(cont_cost, theta_)}, Dict
                                     ({{"compiler", "shell"},
                                                          {"jit", params->jit},
                                                          {"jit_options"
                                                           ".compiler", "gcc"},
                                                          {"jit_options"
                                                           ".flags", "-O3"}}));
    Hess_Contouring = Function("Hess_Contouring", {states, theta_, theta_0, alpha_, x0_,
                                        y0_}, {hessian(cont_cost, theta_)}, Dict
                                         ({{"compiler", "shell"},
                                                             {"jit", params->jit},
                                                             {"jit_options"
                                                              ".compiler", "gcc"},
                                                             {"jit_options"
                                                              ".flags", "-O3"}}));

}

DistributionCost::DistributionCost(SQPParams *params, SplinePath *path):
path(path)
, params(params) {
    if(params->distribution.debug_draw) {
        for (int i = 0; i < params->num_agents; i++) {
            circles.push_back(new sf::CircleShape(0.3));
            circles[i]->setOrigin(0.3, 0.3);
            circles[i]->setFillColor(*params->colors[i]);
        }
    }

}

int DistributionCost::addToVariables(int start_idx) noexcept {
    sl_delta_theta[0] = start_idx;
    sl_delta_theta[1] = sl_delta_theta[0] + params->num_agents * params->N_horizon;

    sl_v[0] = sl_delta_theta[1];
    sl_v[1] = sl_v[0] + params->num_agents * params->N_horizon;

    return 2 * params->num_agents * params->N_horizon;
}

int DistributionCost::addToParams(int start_idx) noexcept {
    sl_theta[0] = start_idx;
    sl_theta[1] = start_idx + params->num_agents * params->N_horizon;
    sl_theta0[0] = sl_theta[1];
    sl_theta0[1] = sl_theta0[0] + params->num_agents * params->N_horizon;
    sl_x0[0] = sl_theta0[1];
    sl_x0[1] = sl_x0[0] + params->num_agents * params->N_horizon;
    sl_y0[0] = sl_x0[1];
    sl_y0[1] = sl_y0[0] + params->num_agents * params->N_horizon;
    sl_alpha[0] = sl_y0[1];
    sl_alpha[1] = sl_alpha[0] + params->num_agents * params->N_horizon;
    sl_xbar[0] = sl_alpha[1];
    sl_xbar[1] = sl_xbar[0] + params->num_agents * params->N_horizon;

    num_params = params->num_agents * params->N_horizon * 6;

    return num_params;
}