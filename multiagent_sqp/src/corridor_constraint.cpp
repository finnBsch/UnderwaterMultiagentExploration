#include "corridor_constraint.h"

void CorridorConstraint::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                       MX &delta_z, MX &delta_u, OptiWrapper *opti_obj,
                                       DM &states, MX &ca_param_vec,
                                       DM &param_vec) {
    buildGradients();
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
    auto fl0 = reshape(ca_param_vec(Slice(sl_fl0[0], sl_fl0[1])),
                       params->num_agents, params->N_horizon);
    auto fr0 = reshape(ca_param_vec(Slice(sl_fr0[0], sl_fr0[1])),
                          params->num_agents, params->N_horizon);
    auto dalpha_l = reshape(ca_param_vec(Slice(sl_dalpha_l[0], sl_dalpha_l[1])),
                                params->num_agents, params->N_horizon);
    auto dalpha_r = reshape(ca_param_vec(Slice(sl_dalpha_r[0], sl_dalpha_r[1])),
                                params->num_agents, params->N_horizon);

    auto maxtheta = ca_param_vec(Slice(sl_maxtheta[0], sl_maxtheta[1]), 0);
    current_theta = DM::zeros(params->num_agents * params->N_horizon, 1);
    current_theta += 0.2;
    auto slack = reshape(opti_obj->ca_variables(Slice(sl_slack[0], sl_slack[1])),
                         params->num_agents, params->N_horizon);
    auto delta_theta = reshape(opti_obj->ca_variables(Slice(sl_delta_theta[0],
                                             sl_delta_theta[1])),
                               params->num_agents, params->N_horizon);
    auto v = reshape(opti_obj->ca_variables(Slice(sl_v[0], sl_v[1])),
                     params->num_agents, params->N_horizon);


    Slice all;
//    opti_obj->subject_to(delta_theta(all, 0) >= 0);
    opti_obj->subject_to(delta_theta(all, 0) == 0);
    for(int t = 1; t < params->N_horizon; t++){
        opti_obj->subject_to(theta(all, t) + delta_theta(all, t) ==
                             delta_theta(all, t - 1) + theta(all, t - 1) + v(all, t - 1) * params->dt);
        opti_obj->subject_to(slack(all, t) >= 0);
        opti_obj->subject_to(0.2 < v(all, t - 1));
        for(int j = 0; j < params->num_agents; j++) {
            Slice state_slice(j*params->num_states, (j+1)*params->num_states);
            auto z = ca_states(state_slice, all);
            auto dz = delta_z(state_slice, all);
            auto dec_vec = vertcat(ca_states(state_slice, t), theta(j, t));
//            auto delta_vec = vertcat(delta_z(state_slice, t), delta_theta(j, t));
            auto delta_vec = delta_z(state_slice, t);

            auto JJ_map = Jac_Cost({{"i0", ca_states(state_slice, t)},
                                    {"i1", theta(j, t)},
                                    {"i2", alpha(j, t)},
                                    {"i3", x0(j, t)},
                                    {"i4", y0(j, t)},
                                    {"i5", theta0(j, t)}});
            auto JJ = JJ_map["o0"];
            auto HH_map = Hess({{"i0", ca_states(state_slice, t)},
                                {"i1", theta(j, t)},
                                {"i2", alpha(j, t)},
                                {"i3", x0(j, t)},
                                {"i4", y0(j, t)},
                                {"i5", theta0(j, t)}});
            auto HH = HH_map["o0"];
            cost_fun +=
                    mtimes(JJ, delta_theta(j, t)) + mtimes(transpose(delta_theta(j, t)),
                                                         mtimes
                                                                 (HH,
                                                                  delta_theta(j, t)));
            auto J_C_L_m = Jac_Constr_L({{"i0", ca_states(state_slice, t)},
                                         {"i1", theta(j, t)},
                                         {"i2", alpha(j, t)},
                                         {"i3", x0(j, t)},
                                         {"i4", y0(j, t)},
                                         {"i5", theta0(j, t)},
                                         {"i6", fl0(j, t)},
                                         {"i7", dalpha_l(j, t)}});
            auto J_C_L = J_C_L_m["o0"];

            auto J_C_R_m = Jac_Constr_R({{"i0", ca_states(state_slice, t)},
                                         {"i1", theta(j, t)},
                                         {"i2", alpha(j, t)},
                                         {"i3", x0(j, t)},
                                         {"i4", y0(j, t)},
                                         {"i5", theta0(j, t)},
                                         {"i6", fr0(j, t)},
                                         {"i7", dalpha_r(j, t)}});
            auto J_C_R = J_C_R_m["o0"];

            auto e_c = (z(0, t) - xbar(j, t)) * sin(alpha(j, t)) - (z(1, t) - (y0(j, t) +
                                                                         sin(alpha(j, t)) *
                                                                         (theta(j, t) -
                                                                          theta0(j, t)))) *
                                                             cos(alpha(j, t));
            auto f_l = fl0(j, t) + sin(dalpha_l(j, t)) * (theta(j, t) - theta0(j, t));
            auto f_r = fr0(j, t) + sin(dalpha_r(j, t)) * (theta(j, t) - theta0(j, t));
            opti_obj->subject_to(
                    e_c - f_r + mtimes(J_C_R, delta_vec) <= pow(slack(j, t),
                                                                1));
            opti_obj->subject_to(
                    f_l - e_c + mtimes(J_C_L, delta_vec) <= pow(slack(j, t),
                                                                1));
            opti_obj->subject_to(slack(j, t) >= 0.0);
//        opti_obj->subject_to(theta(t) + delta_theta(t) <= maxtheta);
            cost_fun += params->corridor.slack_weight * pow(slack(j, t), 2) ;
//        cost_fun += -10.0 * mtimes(J_contouring_cost, delta_vec);
        }
    }

}

int CorridorConstraint::addToParams(int start_idx) noexcept {
    sl_theta[0] = start_idx;
    sl_theta[1] = start_idx + params->N_horizon * params->num_agents;
    sl_theta0[0] = sl_theta[1];
    sl_theta0[1] = sl_theta0[0] + params->N_horizon * params->num_agents;
    sl_x0[0] = sl_theta0[1];
    sl_x0[1] = sl_x0[0] + params->N_horizon * params->num_agents;
    sl_y0[0] = sl_x0[1];
    sl_y0[1] = sl_y0[0] + params->N_horizon * params->num_agents;
    sl_alpha[0] = sl_y0[1];
    sl_alpha[1] = sl_alpha[0] + params->N_horizon * params->num_agents;
    sl_xbar[0] = sl_alpha[1];
    sl_xbar[1] = sl_xbar[0] + params->N_horizon * params->num_agents;
    sl_fl0[0] = sl_xbar[1];
    sl_fl0[1] = sl_fl0[0] + params->N_horizon * params->num_agents;
    sl_fr0[0] = sl_fl0[1];
    sl_fr0[1] = sl_fr0[0] + params->N_horizon * params->num_agents;
    sl_dalpha_l[0] = sl_fr0[1];
    sl_dalpha_l[1] = sl_dalpha_l[0] + params->N_horizon * params->num_agents;
    sl_dalpha_r[0] = sl_dalpha_l[1];
    sl_dalpha_r[1] = sl_dalpha_r[0] + params->N_horizon * params->num_agents;
    sl_maxtheta[0] = sl_dalpha_r[1];
    sl_maxtheta[1] = sl_maxtheta[0] + 1;




    num_params = params->N_horizon * params->num_agents // theta
                 + params->N_horizon * params->num_agents // theta0
                 + params->N_horizon * params->num_agents // x0
                 + params->N_horizon * params->num_agents // y0
                 + params->N_horizon * params->num_agents // alpha
                 + params->N_horizon * params->num_agents // xbar
                 + params->N_horizon * params->num_agents // fl0
                 + params->N_horizon * params->num_agents // fr0
                 + params->N_horizon * params->num_agents // dalpha_l
                 + params->N_horizon * params->num_agents // dalpha_r
                 + 1; // maxtheta
    return num_params;
}

void CorridorConstraint::iterationUpdate(OptiWrapper *opti_obj) {
    DM dthet = opti_obj->value(sl_delta_theta[0], sl_delta_theta[1]);
//    current_theta = current_theta + params->alpha_damping*dthet;
    current_theta = current_theta + dthet;
    opti_obj->update_parameter(sl_theta[0], sl_theta[1], current_theta);

}

void CorridorConstraint::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    for (int t = 0; t < params->N_horizon; t++) {
        for(int i = 0; i < params->num_agents; i++) {
            int index = t * params->num_agents + i;
            auto seg_pts = path->linSegment(current_theta(index).scalar());
            param_vec(sl_xbar[0] + index, 0) = cos(seg_pts[3]) * (current_theta(index)
                                                              - seg_pts[2]) +
                                           seg_pts[0];
            param_vec(sl_theta0[0] + index, 0) = seg_pts[2];
            param_vec(sl_x0[0] + index, 0) = seg_pts[0];
            param_vec(sl_y0[0] + index, 0) = seg_pts[1];
            param_vec(sl_alpha[0] + index, 0) = seg_pts[3];
            param_vec(sl_fl0[0] + index, 0) = seg_pts[4];
            param_vec(sl_fr0[0] + index, 0) = seg_pts[5];
            param_vec(sl_dalpha_l[0] + index, 0) = seg_pts[6];
            param_vec(sl_dalpha_r[0] + index, 0) = seg_pts[7];
//            std::cout << "t : " << t << ", fl0: " << seg_pts[4] << ", fr0: " << seg_pts[5] <<
//            ", theta0 " << seg_pts[2]<< std::endl;
            if(params->corridor.debug_draw) {
                if (t == 0) {
                    circles_0[i]->setPosition(seg_pts[0] + sin(seg_pts[3]) *
                                                           seg_pts[4],
                                              seg_pts[1] - cos(seg_pts[3]) *
                                                           seg_pts[4]);
                    circles_1[i]->setPosition(seg_pts[0] + sin(seg_pts[3]) *
                                                           seg_pts[5],
                                              seg_pts[1] - cos(seg_pts[3]) *
                                                           seg_pts[5]);
                    circles_c[i]->setPosition(seg_pts[0], seg_pts[1]);

                }
            }
        }
    }
    param_vec(Slice(sl_theta[0], sl_theta[1]), 0) = current_theta;
    param_vec(sl_maxtheta[0], 0) = path->getMaxTheta();
}

void CorridorConstraint::warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
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

void CorridorConstraint::draw(sf::RenderTarget &target,
                              sf::RenderStates states) const {
    if(params->corridor.debug_draw) {
        for (auto circ: circles_0) {
            target.draw(*circ, states);
        }
        for (auto circ: circles_1) {
            target.draw(*circ, states);
        }
        for (auto circ: circles_c) {
            target.draw(*circ, states);
        }
    }
}

CorridorConstraint::CorridorConstraint(SQPParams *params, SplinePath*
path): params
                                               (params), path(path){
    if(params->corridor.debug_draw) {
        for (int i = 0; i < params->num_agents; i++) {
            circles_0.push_back(new sf::CircleShape(0.2));
            circles_1.push_back(new sf::CircleShape(0.2));
            circles_c.push_back(new sf::CircleShape(0.25));
            circles_0[i]->setOrigin(0.2, 0.2);
            circles_1[i]->setOrigin(0.2, 0.2);
            circles_c[i]->setOrigin(0.2, 0.2);
            circles_0[i]->setFillColor(*params->colors[i]);
            circles_1[i]->setFillColor(*params->colors[i]);
            circles_c[i]->setFillColor(*params->colors[i]);

        }
    }

}

int CorridorConstraint::addToVariables(int start_idx) noexcept {
    sl_slack[0] = start_idx;
    sl_slack[1] = sl_slack[0] + params->N_horizon * params->num_agents;
    sl_delta_theta[0] = sl_slack[1];
    sl_delta_theta[1] = sl_delta_theta[0] + params->N_horizon * params->num_agents;
    sl_v[0] = sl_delta_theta[1];
    sl_v[1] = sl_v[0] + params->N_horizon * params->num_agents;

    return params->N_horizon * params->num_agents * 3;
}

void CorridorConstraint::buildGradients() {
    MX states = MX::sym("sts", params->num_states);
    auto x = states(0);
    auto y = states(1);
    // Create temporary variables for symbolic differentiation
    MX theta_ = MX::sym("theta_", 1, 1);
    MX theta_0 = MX::sym("theta_0", 1, 1);
    MX alpha_ = MX::sym("alpha_", 1, 1);
    MX fr0_ = MX::sym("fr0_", 1, 1);
    MX fl0_ = MX::sym("fl0_", 1, 1);
    MX d_alpha_r = MX::sym("d_alpha_r", 1, 1);
    MX d_alpha_l = MX::sym("d_alpha_l", 1, 1);
    MX x0_ = MX::sym("x0_", 1, 1);
    MX y0_ = MX::sym("y0_", 1, 1);

    MX x_ = cos(alpha_) * (theta_ - theta_0) + x0_;
    auto e_l = -(x - x_) * cos(alpha_) - (y - (y0_
                                                       + sin(alpha_) * (theta_ - theta_0))) *
                                             sin(alpha_);
    auto e_c = (x - x_) * sin(alpha_) - (y - (y0_
                                                      + sin(alpha_) * (theta_ - theta_0))) * cos(alpha_);

    auto f_r = fr0_ + sin(d_alpha_r) * (theta_ - theta_0);
    auto f_l = fl0_ + sin(d_alpha_l) * (theta_ - theta_0);
    auto cont_constr_r = e_c - f_r;
    auto cont_constr_l = f_l - e_c;
    auto cost = pow(e_l, 2) * params->corridor.lag_weight;
//    auto vars_deriv = vertcat(states, theta_);
    auto vars_deriv = states;
    auto Jac_Const_L = simplify(jacobian(cont_constr_l, vars_deriv));
    auto Jac_Const_R = simplify(jacobian(cont_constr_r, vars_deriv));
    auto Jac_ = simplify(jacobian(cost, theta_));

    auto Hess_ = simplify(hessian(cost, theta_));
    // Create function that take linearized states and returns Jacobian.
    Jac_Cost = Function("jac_con", {states, theta_, alpha_, x0_, y0_, theta_0},
                        {Jac_}, Dict({{"compiler", "shell"},
                                      {"jit", params->jit},
                                      {"jit_options"
                                       ".compiler", "gcc"},
                                      {"jit_options"
                                       ".flags", "-O3"}}));


    Jac_Constr_L = Function("jac_constr_l", {states, theta_, alpha_, x0_, y0_,
                                             theta_0, fl0_, d_alpha_l},
                            {Jac_Const_L}, Dict({{"compiler", "shell"},
                                                 {"jit", params->jit},
                                                 {"jit_options"
                                                  ".compiler", "gcc"},
                                                 {"jit_options"
                                                  ".flags", "-O3"}}));
    Jac_Constr_R = Function("jac_constr_r", {states, theta_, alpha_, x0_, y0_,
                                             theta_0, fr0_, d_alpha_r},
                            {Jac_Const_R}, Dict({{"compiler", "shell"},
                                                 {"jit", params->jit},
                                                 {"jit_options"
                                                  ".compiler", "gcc"},
                                                 {"jit_options"
                                                  ".flags", "-O3"}}));
    // Create function that take linearized states and returns Hessian.
    Hess = Function("hess_con", {states, theta_, alpha_, x0_, y0_, theta_0},
                    {simplify(Hess_)}, Dict({{"compiler", "shell"},
                                             {"jit", params->jit},
                                             {"jit_options"
                                              ".compiler", "gcc"},
                                             {"jit_options"
                                              ".flags", "-O3"}}));

}

double CorridorConstraint::getTheta() const {
    return current_theta(0).scalar();
}
