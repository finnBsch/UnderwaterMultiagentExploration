//
// Created by finn on 6/19/23.
//
#include "bluerov_dynamics.h"

#define x_ 0
#define y_ 1
#define phi_ 2
#define vx_ 3
#define vy_ 4

#define dphi_ 0
#define a_x_ 1
#define a_y_ 2


void BlueROVDynamics::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                    MX &delta_z, MX &delta_u, OptiWrapper *opti_obj,
                                    DM &states, MX &ca_param_vec,
                                    DM &param_vec) {
    Slice all;
    Slice state_slice(id * num_states, (id + 1) * num_states);
    Slice input_slice(id * num_inputs, (id + 1) * num_inputs);

    auto z = ca_states(state_slice, all);
    auto u = ca_u(input_slice, all);

    auto dz = delta_z(state_slice, all);
    auto du = delta_u(input_slice, all);
    opti_obj->subject_to(dz(all, 0) == 0);
    // Initial constraints
    opti_obj->subject_to(dz(all, 0) == 0);
    for(int i = 0; i < params->N_horizon - 1; i++){
        opti_obj->subject_to(-a_params.accel_max <= u(a_x_, i) + du(a_x_, i) <=
                             a_params
                                     .accel_max);
        opti_obj->subject_to(-a_params.accel_max <= u(a_y_, i) + du(a_y_, i) <=
                                a_params
                                .accel_max);


        opti_obj->subject_to(-a_params.yaw_rate_max <= u(dphi_, i) + du
                (dphi_, i) <= a_params.yaw_rate_max);
        auto z1 = dynamics(MXDict({{"i0", z(all,
                                            i)},
                                   {"i1", u(all, i)}}));
        auto grad = dynamics_grad(MXDict({{"i0", z(all,
                                                   i)},
                                          {"i1", u(all, i)},
                                          {"i2", z(all, i + 1)}}));
        auto dparam = vertcat(dz(all, i), du(all, i), dz(all, i + 1));
        MX xx = z(all, i + 1) - z1["o0"] + mtimes(grad["o0"],
                                                  dparam);
//        std::cout << "XX " << xx << "\n";
//        std::cout << xx.size() << "\n";
        opti_obj->subject_to(z(all, i + 1) - z1["o0"] + mtimes(grad["o0"],
                                                               dparam) == 0);
        // TODO Speed max

    }

}

int BlueROVDynamics::addToParams(int start_idx) noexcept {
    return 0;
}

void BlueROVDynamics::iterationUpdate(OptiWrapper *opti_obj) {

}

void BlueROVDynamics::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {

}

void BlueROVDynamics::warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {

}

void
BlueROVDynamics::draw(sf::RenderTarget &target, sf::RenderStates states) const {

}

BlueROVDynamics::BlueROVDynamics(SQPParams *params, int id): params(params){
    this->id = id;
    buildCasadiGradients();
}

void BlueROVDynamics::buildCasadiGradients() {
    MX z = MX::sym("z", num_states);
    MX u = MX::sym("u", num_inputs);
    MX z1 = MX::sym("z1", num_states);
    double dt = params->dt;

    auto d_phi = u(dphi_);
    auto a_x = u(a_x_);
    auto a_y = u(a_y_);

    auto x = z(x_);
    auto y = z(y_);
    auto phi = z(phi_);
    auto v_x = z(vx_);
    auto v_y = z(vy_);


    auto x_1 = x + v_x * dt;
    auto y_1 = y + v_y * dt;
    auto phi_1 = phi + d_phi * dt;
    auto v_x_1 = v_x + (a_x - v_x * a_params.damping_factor) * dt;
    auto v_y_1 = v_y + (a_y - v_y * a_params.damping_factor) * dt;
    auto dyn_out = vertcat(x_1, y_1, phi_1, v_x_1, v_y_1);
    dynamics = Function("dynamics", {z, u}, {dyn_out}, Dict({{"compiler", "shell"},
                                                             {"jit", params->jit},
                                                             {"jit_options"
                                                              ".compiler", "gcc"},
                                                             {"jit_options"
                                                              ".flags", "-O3"}}));
    auto dyn_error = z1 - dyn_out;
    auto all_params = vertcat(z, u, z1);
    dynamics_grad = Function("dynamics_grad", {z, u, z1}, {jacobian
                                                                   (dyn_error, all_params)}, Dict({{"compiler", "shell"},
                                                                                                   {"jit", params->jit},
                                                                                                   {"jit_options"
                                                                                                    ".compiler", "gcc"},
                                                                                                   {"jit_options"
                                                                                                    ".flags", "-O3"}}));
}

int BlueROVDynamics::addToVariables(int start_idx) noexcept {
    return 0;
}
