//
// Created by finn on 5/31/23.
//
#include "com_path_tracking.h"

COMPathTracking::COMPathTracking(SQPParams *params, SplinePath* path) : params
(params),
                                                      ct_pt(0.1), path(path) {
    if (params->viz_all) {
        for (int i = 0; i < params->N_horizon; i++) {
            lin_pts.emplace_back(0.05 + (params->N_horizon - i) * 0.005);
            lin_pts[i].setFillColor(sf::Color::Blue);
        }
        ct_pt.setFillColor(sf::Color::Blue);
    }
}


int COMPathTracking::addToVariables(int start_idx) noexcept {
    sl_v[0] = start_idx;
    sl_v[1] = start_idx + params->N_horizon;

    sl_delta_theta[0] = sl_v[1];
    sl_delta_theta[1] = sl_delta_theta[0] + params->N_horizon;

    return params->N_horizon * 2;
}

void COMPathTracking::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                    MX &delta_z, MX &delta_u,
                                    OptiWrapper *opti_obj, DM &states,
                                    MX &ca_param_vec, DM &param_vec) {
    deriveGradients();
    // Coords of lin. segment
    auto x0 = ca_param_vec(Slice(sl_x0[0], sl_x0[1]), 0); // origin x
    auto y0 = ca_param_vec(Slice(sl_y0[0], sl_y0[1]), 0); // origin y
    auto theta0 = ca_param_vec(Slice(sl_theta0[0], sl_theta0[1]), 0); // origin
    auto alpha = ca_param_vec(Slice(sl_alpha[0], sl_alpha[1]), 0); // angle of
    // segment
    auto theta = ca_param_vec(Slice(sl_theta[0], sl_theta[1]), 0); // angle of
    // segment
    auto xbar = ca_param_vec(Slice(sl_xbar[0], sl_xbar[1]), 0); // angle of
    // segment
    auto theta_0_constr = ca_param_vec(Slice(sl_theta0_constr[0],
                                             sl_theta0_constr[1]), 0); //
                                             // angle of segment
    current_theta = DM::zeros(params->N_horizon, 1); // DM to hold current theta
    // lin
    // . pt
    Slice all(0, params->num_states * params->num_agents);
    // Optim variables
    auto delta_theta = opti_obj->ca_variables(Slice(sl_delta_theta[0],
                                             sl_delta_theta[1]));

    auto v = opti_obj->ca_variables(Slice(sl_v[0], sl_v[1]));

//    opti_obj->subject_to(delta_theta(0) + theta(0) == theta_0_constr);
    double x_com = 0;
    double y_com = 0;
    for (int i = 0; i < params->num_agents; i++) {
        x_com += states(i * params->num_states, 0).scalar();
        y_com += states(i * params->num_states + 1, 0).scalar();
    }
    x_com /= params->num_agents;
    y_com /= params->num_agents;
    double theta_0_current = 0.0;
    current_theta(0) = theta_0_current;
//    std::cout << theta_0_current << "\n";
    for (int i = 1; i < params->N_horizon; i++) {
        current_theta(i) = current_theta(0) + params->dt * 0.1 * i;
        auto dec_vec = vertcat(ca_states(all, i), theta(i));
        auto delta_vec = vertcat(delta_z(all, i), delta_theta(i));
        // vars_deriv, alpha_, x0_, y0_
        auto JJ_map = Jac({{"i0", dec_vec},
                           {"i1", alpha(i)},
                           {"i2", x0(i)},
                           {"i3", y0(i)},
                           {"i4", theta0(i)}});
        auto JJ = JJ_map["o0"];
        auto HH_map = Hess({{"i0", dec_vec},
                            {"i1", alpha(i)},
                            {"i2", x0(i)},
                            {"i3", y0(i)},
                            {"i4", theta0(i)}});
        auto HH = HH_map["o0"];

        cost_fun += mtimes(JJ, delta_vec)
                + mtimes(transpose(delta_vec),
                                                         mtimes
                                                                 (HH,
                                                                  delta_vec));

        cost_fun -= params->tracking.progress_weight * v(i - 1) *
                params->dt;
        opti_obj->subject_to(delta_theta(i) + theta(i) == delta_theta(i - 1)
                                                          + theta(i - 1) +
                                                          v(i - 1) *
                                                          params->dt);
        opti_obj->subject_to(-.0 < v(i - 1) <= params->tracking.max_progress);
//        opti_obj->subject_to(delta_theta(i)+ theta(i) >= 0.0);

    }
    opti_obj->subject_to(delta_theta(0) == 0);
    param_vec(sl_theta0_constr[0]) = theta_0_current;
    linearize(opti_obj, states, param_vec);
}

void COMPathTracking::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
//    std::cout << "XXX" << current_theta << "\n";
    for(int i = 0; i < params->N_horizon - 1; i++){
        if(current_theta(i + 1).scalar() < current_theta(i).scalar()){
//            std::cout << "WRAP\n";
        }
    }
    param_vec(Slice(sl_theta[0], sl_theta[1]), 0) = current_theta;

    for (int i = 0; i < params->N_horizon; i++) {
        auto seg_pts = path->linSegment(current_theta(i, 0).scalar());
        param_vec(sl_xbar[0] + i, 0) = cos(seg_pts[3]) * (current_theta(i)
                                                        - seg_pts[2]) +
                                     seg_pts[0];
        param_vec(sl_theta0[0] + i, 0) = seg_pts[2];
        param_vec(sl_x0[0] + i, 0) = seg_pts[0];
        param_vec(sl_y0[0] + i, 0) = seg_pts[1];
        param_vec(sl_alpha[0] + i, 0) = seg_pts[3];
//        std::cout << "Step " << i << ", angle: " << seg_pts[3] * 180 / M_PI <<
//                  ", dtheta: " << current_theta(i) - seg_pts[2] << ", dxbar: "
//                  << cos
//                             (seg_pts[3]) * (current_theta(i)
//                                             - seg_pts[2]) << ", dybar: "
//                  << tan(seg_pts[3]) * (cos(seg_pts[3]) * (current_theta(i)
//                                                           - seg_pts[2]))
//                  << ", xbar: " << cos(seg_pts[3]) * (current_theta(i)
//                                                      - seg_pts[2]) + seg_pts[0]
//                  << "\n";

    }
}

void COMPathTracking::iterationUpdate(OptiWrapper *opti_obj) {
    DM dthet = opti_obj->value(sl_delta_theta[0], sl_delta_theta[1]);
    current_theta = current_theta + params->alpha_damping*dthet;
    opti_obj->update_parameter(sl_theta[0], sl_theta[1], current_theta);
}

void COMPathTracking::warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    Slice first(0, params->N_horizon - 1);
    Slice second(1, params->N_horizon);
    // TODO Fix this constraint!! Find nearest point on path function required.
    double x_com = 0;
    double y_com = 0;
    for (int i = 0; i < params->num_agents; i++) {
        x_com += states(i * params->num_states).scalar();
        y_com += states(i * params->num_states + 1).scalar();
    }
    x_com /= params->num_agents;
    y_com /= params->num_agents;
//    double theta_0_current = path->getTheta(x_com, y_com, current_theta(0)
//    .scalar() - 0.1, current_theta(0).scalar() + params->tracking.max_progress*
//    params->dt);
//    param_vec(sl_theta0_constr[0]) = theta_0_current;
//    std::cout << "MMM " << theta_0_current << "\n";
//    std::cout << current_theta << "\n";
//    if(current_theta(1).scalar() < theta_0_current){
//        current_theta(second) += theta_0_current - current_theta(1).scalar();
//    }
    current_theta(first) = current_theta(second);
//    current_theta(0) = theta_0_current;
    current_theta(params->N_horizon - 1) = current_theta(params->N_horizon - 2)
                                           * 2
                                           -
                                           current_theta(params->N_horizon - 3);
    if (params->viz_all) {
        ct_pt.setPosition(x_com, y_com);
        for (int i = 0; i < params->N_horizon; i++) {
            auto seg_pts = path->linSegment(current_theta(i, 0).scalar());
            double xbar_local = (cos(seg_pts[3]) * (current_theta(i)
                                                    - seg_pts[2]) +
                                 seg_pts[0]).scalar();
            double ybar_local = (sin(seg_pts[3]) * (current_theta(i)
                                                    - seg_pts[2]) +
                                 seg_pts[1]).scalar();
            lin_pts[i].setPosition((float) xbar_local, (float) ybar_local);
        }
    }
}

void
COMPathTracking::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for (const auto &c: lin_pts) {
        target.draw(c, states);
    }
//    target.draw(lin_pts[0], states);
    target.draw(ct_pt);
}

void COMPathTracking::deriveGradients() {
    MX states = MX::sym("sts", params->num_agents * params->num_states);
    auto x_com = sum1(states(Slice(0,
                                   (params->num_agents - 1)
                                   * params->num_states + 1,
                                   params->num_states), 0)) /
                 params->num_agents;
    auto y_com = sum1(states(
            Slice(1, (params->num_agents - 1) * params->num_states + 2,
                  params->num_states), 0)) / params->num_agents;
    // Create temporary variables for symbolic differentiation
    MX theta_ = MX::sym("theta_", 1, 1);
    MX theta_0 = MX::sym("theta_0", 1, 1);
    MX alpha_ = MX::sym("alpha_", 1, 1);
    MX x0_ = MX::sym("x0_", 1, 1);
    MX y0_ = MX::sym("y0_", 1, 1);

    MX x_ = cos(alpha_) * (theta_ - theta_0) + x0_;
    // Error functions
    auto e_l = -(x_com - x_) * cos(alpha_) - (y_com - (y0_
                                              + sin(alpha_) * (theta_ - theta_0))) *
                                             sin(alpha_);
    auto e_c = (x_com - x_) * sin(alpha_) - (y_com - (y0_
                                                      + sin(alpha_) * (theta_ - theta_0))) * cos(alpha_);
    auto cost = pow(e_l, 2) * params->tracking.lag_weight
                + pow(e_c, 2) * params->tracking.contouring_weight;
    // optimization variables
    auto vars_deriv = vertcat(states, theta_);
    auto Jac_ = simplify(jacobian(cost, vars_deriv));
    auto Hess_ = simplify(hessian(cost, vars_deriv));
//    std::cout << "Hess_" << Hess_ << "\n";
    // Create function that take linearized states and returns Jacobian.
    Jac = Function("jac_con", {vars_deriv, alpha_, x0_, y0_, theta_0}, {Jac_}, Dict({{"compiler", "shell"},
                                                                                     {"jit", params->jit},
                                                                                     {"jit_options"
                                                                                      ".compiler", "gcc"},
                                                                                     {"jit_options"
                                                                                      ".flags", "-O3"}}));
    // Create function that take linearized states and returns Hessian.
    Hess = Function("hess_con", {vars_deriv, alpha_, x0_, y0_, theta_0},
                    {simplify(Hess_)}, Dict({{"compiler", "shell"},
                                             {"jit", params->jit},
                                             {"jit_options"
                                              ".compiler", "gcc"},
                                             {"jit_options"
                                              ".flags", "-O3"}}));
}

int COMPathTracking::addToParams(int start_idx) noexcept {
    sl_x0[0] = start_idx;
    sl_x0[1] = start_idx + params->N_horizon;
    sl_y0[0] = start_idx + params->N_horizon;
    sl_y0[1] = start_idx + 2 * params->N_horizon;
    sl_theta0[0] = start_idx + 2 * params->N_horizon;
    sl_theta0[1] = start_idx + 3 * params->N_horizon;
    sl_alpha[0] = start_idx + 3 * params->N_horizon;
    sl_alpha[1] = start_idx + 4 * params->N_horizon;
    sl_theta[0] = start_idx + 4 * params->N_horizon;
    sl_theta[1] = start_idx + 5 * params->N_horizon;
    sl_xbar[0] = start_idx + 5 * params->N_horizon;
    sl_xbar[1] = start_idx + 6 * params->N_horizon;

    sl_theta0_constr[0] = start_idx + 6 * params->N_horizon;
    sl_theta0_constr[1] = start_idx + 6 * params->N_horizon + 1;
    num_params = params->N_horizon  // x0
            + params->N_horizon  // y0
            + params->N_horizon  // theta0
            + params->N_horizon  // alpha
            + params->N_horizon  // theta
            + params->N_horizon // xbar
            + 1 ;  // theta0constaint
    return num_params;
}

double COMPathTracking::getTheta() const {
    return current_theta(params->N_horizon - 1).scalar();
}
