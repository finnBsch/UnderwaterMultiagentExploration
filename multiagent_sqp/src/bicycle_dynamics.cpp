//
// Created by finn on 5/18/23.
//

#include "bicycle_dynamics.h"

#define x_ 0
#define y_ 1
#define phi_ 2
#define v_ 3

/*
 * State x = [x, y, angle, v_body]
 * Input u = [yaw_rate, angle]
 */
using namespace casadi;
MX BicycleDynamics::getDeriv(const MX &state, const MX &u) {
    return vertcat(state(v_) * cos(state(phi_)),
                   state(v_) * sin(state(phi_)),
                   u(0),
                   u(1));
}

BicycleDynamics::BicycleDynamics(SQPParams* params, int id):  params(params) {
    this->id = id;
}

e_state BicycleDynamics::getDeriv(const e_state &state, const e_u &u) {
    return e_state(state(3) * cos(state(2)), state(3) * sin(state(2)), u(
            0
            ), u(1));
}

int BicycleDynamics::addToVariables(int start_idx) noexcept {
    return 0;
}

void BicycleDynamics::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                                    MX &delta_u,
                                    OptiWrapper *opti_obj, DM &states, MX &ca_param_vec, DM &param_vec) {
    mu_vals = DM::zeros(num_states, params->N_horizon - 1);
//    mu = opti_obj->parameter(num_states, params->N_horizon - 1);
//    dmu = opti_obj->variable(num_states, params->N_horizon - 1);
//    x_init = opti_obj->parameter(num_states, 1);
    Slice all;
    Slice state_slice(id*num_states, (id+1)*num_states);
    Slice u_slice(id*num_inputs, (id+1)*num_inputs);
    // First get and reshape this shit 8),
    // ca_states has dim (num_agents*num_states + num_lagrange,
    // params->N_horizon)
    // We need the states in (id*num_states:(id+1)*num_states,
    // params->N_horizon)
    auto z = ca_states(state_slice, all);
    auto dz = delta_z(state_slice, all);
    auto u = ca_u(u_slice, all);
    auto du = delta_u(u_slice, all);
    // System dynamics first. We need: all eq.-Constraints g for all time
    // steps. We need the gradient of these w.r.t. the states. We need the
//    MX init_1 = z(x_, 0) - x_init(x_, 0) + dz(x_, 0);
//    MX init_2 = z(y_, 0) - x_init(y_, 0) + dz(y_, 0);
//    MX init_3 = z(phi_, 0) - x_init(phi_, 0) + dz(phi_, 0);
//    MX init_4 = z(v_, 0) - x_init(v_, 0) + dz(v_, 0);
    MX init_1 = dz(x_, 0);
    MX init_2 = dz(y_, 0);
    MX init_3 = dz(phi_, 0);
    MX init_4 = dz(v_, 0);
    opti_obj->subject_to(init_1 == 0);
    opti_obj->subject_to(init_2 == 0);
    opti_obj->subject_to(init_3 == 0);
    opti_obj->subject_to(init_4 == 0);
    for(int i = 0; i < params->N_horizon - 1; i++) {
        opti_obj->subject_to(0 <= z(x_, i + 1) + dz(x_, i + 1) <=
        params->size_x);
        opti_obj->subject_to(0 <= z(y_, i + 1) + dz(y_, i + 1) <=
        params->size_y);
        opti_obj->subject_to(-a_params.speed_max <= z(v_, i) + dz(v_ ,i)  <=
        a_params.speed_max);
        MX g1 = z(x_, i + 1) - z(x_, i) - params->dt*cos(z(phi_, i)) * z(v_, i)
                // g
                // dphik*params->dt*vk*sin(phik) - params->dt*dvk*cos(phik) -
                // dxk +
                // dxk1
                + dz(phi_, i) * params->dt * z(v_, i) * sin(z(phi_, i))
                - params->dt * dz(v_, i) * cos(z(phi_, i))
                - dz(x_, i) + dz(x_, i + 1);
        MX g2 = z(y_, i + 1) - z(y_, i) - params->dt*sin(z(phi_, i)) * z(v_, i)
                // g
                // - dphik*params->dt*vk*cos(phik) - params->dt*dvk*sin(phik) -
                // dyk +
                // dyk1
                - dz(phi_, i) * params->dt * z(v_, i) * cos(z(phi_, i))
                - params->dt * dz(v_, i) * sin(z(phi_, i))
                - dz(y_, i) + dz(y_, i + 1);
        MX g3 = z(phi_, i + 1) - z(phi_, i) - params->dt*u(0, i) // g
                //  - dphik + dphik1 - params->dt*du0
                - dz(phi_, i) + dz(phi_, i + 1) - params->dt * du(0, i);
        MX g4 = z(v_, i + 1) - z(v_, i) - params->dt * u(1, i) // g
                // - params->dt*du1 - dvk + dvk1
                - params->dt * du(1, i) - dz(v_, i) + dz(v_, i + 1);
        MX gu1 = u(0, i) + du(0, i);
        MX gu2 = u(1, i) + du(1, i);
        opti_obj->subject_to(g1 == 0);
        opti_obj->subject_to(g2 == 0);
        opti_obj->subject_to(g3 == 0);
        opti_obj->subject_to(g4 == 0);
        opti_obj->subject_to( -a_params.yaw_rate_max < gu1 < a_params.yaw_rate_max);
        opti_obj->subject_to(-a_params.accel_max < gu2 < a_params.accel_max);
        // dphik**2*params->dt*vk*xmuk*cos(phik)
//        MX c1 = pow(dz(phi_, i), 2) * params->dt * z(v_, i) * mu(x_, i) * cos(z
//                (phi_,i))
//                //+ 2*dphik*params->dt*dvk*xmuk*sin(phik)
//                + 2 * dz(phi_, i) * params->dt * dz(v_, i) * mu(x_, i) * sin(z
//                (phi_,
//                                                                        i))
//                //+ 2*dphik*params->dt*dxmuk*vk*sin(phik)
//                + 2 * dz(phi_, i) * params->dt * dmu(x_, i) * z(v_, i) * sin(z
//                (phi_,
//                                                                        i))
//                // - 2*params->dt*dvk*dxmuk*cos(phik)
//                - 2 * params->dt * dz(v_, i) * dmu(x_, i) * cos(z(phi_, i))
//                //  - 2*dxk*dxmuk + 2*dxk1*dxmuk
//                - 2 * dz(x_, i) * dmu(x_, i) + 2 * dz(x_, i + 1) * dmu(x_, i);
//        // dphik**2*params->dt*vk*ymuk*sin(phik)
//        MX c2 = pow(dz(phi_, i), 2) * params->dt * z(v_, i) * mu(y_, i) * sin(z
//                (phi_,i))
//                //  + 2*dphik*params->dt*dvk*ymuk*cos(phik)
//                + 2 * dz(phi_, i) * params->dt * dz(v_, i) * mu(y_, i) * cos(z
//                (phi_,
//                                                                        i))
//                //  + 2*dphik*params->dt*dymuk*vk*cos(phik)
//                + 2 * dz(phi_, i) * params->dt * dmu(y_, i) * z(v_, i) * cos(z
//                (phi_,
//                                                                        i))
//                // + 2*params->dt*dvk*dymuk*sin(phik)
//                + 2 * params->dt * dz(v_, i) * dmu(y_, i) * sin(z(phi_, i))
//                // - 2*dyk*dymuk + 2*dyk1*dymuk
//                - 2 * dz(y_, i) * dmu(y_, i) + 2 * dz(y_, i + 1) * dmu(y_, i);
//        // 2*dphimuk*(-dphik + dphik1 - params->dt*du0)
//        MX c3 = 2 * dmu(phi_, i)
//                * ( - dz(phi_, i) + dz(phi_, i + 1) - params->dt * du(0,i));
//        // 2*dvmuk*( - params->dt*du1 - dvk + dvk1)
//        MX c4 = 2 * dmu(v_, i) * ( - params->dt * du(1, i) - dz(v_, i) + dz(v_,
//                                                                          i +
//        1));
        // TEMPORARY
//        MX c_pos = pow(z(x_, i + 1) + dz(x_, i + 1) - (id)*2  +0.1, 2)
//                + pow(z(y_, i + 1) + dz(y_,i + 1) -  10.0, 2);
//        MX c_pos = z(y_, i + 1) + dz(y_, i + 1);
//        cost_fun += -(c1 + c2 + c3 + c4) * 0.5 + c_pos;
        cost_fun += 0.1 * pow(u(0, i) + du(0,i), 2) + 5.0* pow(u(1, i) + du
                (1,i), 2);

    }
}

void BicycleDynamics::iterationUpdate(OptiWrapper *opti_obj) {
//    auto sol_dmu = opti_obj->value(dmu);
//    mu_vals += sol_dmu;
}

void BicycleDynamics::initializeNextSol(Opti* opti_obj,
                                        Eigen::Matrix<double, Eigen::Dynamic, 1> init_state) {
    // TODO REMOVE
//    DM x_init_dm = DM::zeros(num_states, 1);
//    std::memcpy(x_init_dm.ptr(), init_state.data(), sizeof(double)*init_state
//    .size());
//    opti_obj->set_value(x_init, x_init_dm);
//    opti_obj->set_value(mu, 0);
}

void BicycleDynamics::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
//    opti_obj->set_value(mu, mu_vals);
}
