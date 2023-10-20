//
// Created by finn on 5/29/23.
//
#include "path_tracking.h"
#define x_ 0
#define y_ 1
#define phi_ 2
#define v_ 3

// TODO: Make all params one vector. Turn individual params into views of
//  that vector. Then, all opti_obj->set_value() calls will be merged into
//  one set_value, reducing the huge overhead.
void
PathTracking::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                            MX &delta_u,
                            OptiWrapper *opti_obj, DM &states, MX &ca_param_vec,
                            DM &param_vec) {

    // Coords of lin. segment
    x0 = opti_obj->parameter(params->N_horizon, 1); // origin x
    y0 = opti_obj->parameter(params->N_horizon, 1); // origin y
    alpha = opti_obj->parameter(params->N_horizon, 1); // angle
    theta0 = opti_obj->parameter(params->N_horizon, 1);  // theta of origin

    // Values at linearized theta
    theta = opti_obj->parameter(params->N_horizon, 1);
    xbar = opti_obj->parameter(params->N_horizon, 1);

    current_theta = DM::zeros(params->N_horizon, 1); // DM to hold current theta
    // lin
    // . pt

    theta_0_constr = opti_obj->parameter(1, 1);
    // Optim variables
    delta_theta = opti_obj->variable(params->N_horizon, 1);
    v = opti_obj->variable(params->N_horizon, 1);

    Slice state_slice(id*params->num_states, (id+1)*params->num_states);
    Slice all;
    auto z = ca_states(state_slice, all);
    auto dz = delta_z(state_slice, all);
//    MX e_c_sq = pow((z(x_, all) - xbar)*sin(alpha) - (z(y_, all) - y0
//            + (x0 -xbar)*tan(alpha))*cos(alpha), 2);  NOT required, constant
//            because doesn't contain delta_z
    for(int i = 1; i < params->N_horizon; i++){
        current_theta(i) = params->dt * 0.5 * i;
        // Contouring
        // dec^2/dx * delta_x
        cost_fun += params->tracking.contouring_weight*dz(x_, i)
                * (2*((z(x_, i)  - xbar(i))*sin(alpha(i)) - (z(y_, i) - y0(i) +
                (x0(i)- xbar(i))*tan(alpha(i)))*cos(alpha(i)))*sin(alpha(i)));
        // dec^2/dy * delta_y
        cost_fun += params->tracking.contouring_weight*dz(y_, i)
                    * (-2*((z(x_, i) - xbar(i))*sin(alpha(i)) - (z(y_, i) -
                    y0(i)
                    + (x0(i)- xbar(i))*tan(alpha(i)))*cos(alpha(i)))*cos
                    (alpha(i)));
        // dec^2/dtheta * delta_theta
        cost_fun += 0;
        // Quadratic terms:
        cost_fun += 0.5 * params->tracking.contouring_weight*(-pow(dz(x_, i), 2)
                *cos(2*alpha(i)) +
                pow(dz
                (x_,
                                                                          i),
                                                                  2) -
                2*dz(x_, i)*dz(y_, i)*sin
                (2*alpha(i)) + pow(dz(y_, i), 2)*cos(2*alpha(i)) + pow(dz(y_,
                                                                          i),
                                                                       2));
        // LAG
        // del^2/dx * delta_x
        cost_fun += params->tracking.lag_weight * dz(x_, i)
                    * (2*((z(x_, i)  - xbar(i))*cos(alpha(i)) + (z(y_, i) - y0
                    (i) + (x0(i)- xbar(i))*tan(alpha(i)))*sin(alpha(i)))*cos
                    (alpha(i)));
        // del^2/dy * delta_y
        cost_fun += params->tracking.lag_weight * dz(y_, i)
                    * (2*((z(x_, i)  - xbar(i))*cos(alpha(i)) + (z(y_, i) -
                    y0(i)
                    + (x0(i)- xbar(i))*tan(alpha(i)))*sin(alpha(i)))*sin
                    (alpha(i)));
        cost_fun += params->tracking.lag_weight * delta_theta(i)
                    *( -2*((z(x_, i) - xbar(i))*cos(alpha(i))
                    + (z(y_, i) - y0(i)
                    + (x0(i)- xbar(i))*tan(alpha(i)))*sin(alpha(i)))
                    *cos(alpha(i))/cos(alpha(i)));
        // quadratic
        cost_fun += 0.5 * params->tracking.lag_weight * ((-2*delta_theta(i)*
                (-delta_theta(i) + dz(x_, i)
                *cos(alpha(i))
                +dz(y_, i)*sin(alpha(i)))
                *pow(cos(alpha(i)), 2) +
                (-2*delta_theta(i)*dz(x_, i)*cos(alpha(i))
                - 2*delta_theta(i)*dz(y_, i)*sin(alpha(i))
                + pow(dz(x_, i), 2)*cos(2*alpha(i)) + pow(dz(x_, i), 2)
                + 2*dz(x_, i)*dz(y_, i)*sin(2*alpha(i))
                - pow(dz(y_, i), 2)*cos(2*alpha(i))
                +pow(dz(y_, i), 2))*pow(cos(alpha(i)),2))/pow(cos(alpha(i)),2));



        cost_fun -= params->tracking.progress_weight * v(i - 1) * params->dt;
        opti_obj->subject_to(delta_theta(i ) + theta(i) == delta_theta(i - 1)
        + theta(i -1)+ v(i - 1) * params->dt);
        opti_obj->subject_to(-1.0 < v(i - 1) <= 100.0);

    }
    current_theta(0) = 0;
    opti_obj->subject_to(delta_theta(0) + theta(0) == theta_0_constr);
    opti_obj->set_value(theta, current_theta);
    opti_obj->set_value(theta_0_constr, 0);
    linearize(opti_obj, states, param_vec);
}

PathTracking::PathTracking(SQPParams* params, int id):
                params(params), id(id) {
    if(params->viz_all){
        for(int i = 0; i < params->N_horizon; i++){
            lin_pts.emplace_back(0.05);
        }
    }
}

void PathTracking::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    std::cout << "XXX" << current_theta <<"\n";
    opti_obj->set_value(theta, current_theta);

    for(int i = 0; i < params->N_horizon; i++){
        auto seg_pts = path.linSegment(current_theta(i, 0).scalar());
        opti_obj->set_value(xbar(i), cos(seg_pts[3]) * (current_theta(i)
        -seg_pts[2]) + seg_pts[0]);
        std::cout << "Step " << i << ", angle: " << seg_pts[3] * 180/M_PI <<
        ", dtheta: " <<current_theta(i)-seg_pts[2] << ", dxbar: " << cos
        (seg_pts[3]) * (current_theta(i) -seg_pts[2]) << ", dybar: " << tan
        (seg_pts[3]) * (cos(seg_pts[3]) * (current_theta(i)-seg_pts[2]) )
        << ", xbar: " << cos(seg_pts[3]) * (current_theta(i)-seg_pts[2]) + seg_pts[0] << "\n";
        opti_obj->set_value(theta0(i), seg_pts[2]);
        opti_obj->set_value(x0(i), seg_pts[0]);
        opti_obj->set_value(y0(i), seg_pts[1]);
        opti_obj->set_value(alpha(i), seg_pts[3]);
    }
}

void PathTracking::iterationUpdate(OptiWrapper *opti_obj) {
    // Fetch the current values for delta_theta and update linearization point
    DM dthet = opti_obj->value(delta_theta);
    current_theta = current_theta + dthet;

}

void PathTracking::warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    Slice first(0, params->N_horizon - 1);
    Slice second(1, params->N_horizon);
    opti_obj->set_value(theta_0_constr, path.getTheta(states(id *
    params->num_states).scalar(),
                                                      states(id *
                                                      params->num_states +1
                                                      ).scalar()));
    current_theta(first) = current_theta(second);
    current_theta(params->N_horizon - 1) = current_theta(params->N_horizon - 2)
            * 2
            - current_theta(params->N_horizon - 3);
    if(params->viz_all){
        for(int i = 0; i < params->N_horizon; i++) {
            auto seg_pts = path.linSegment(current_theta(i, 0).scalar());
            double xbar_local = (cos(seg_pts[3]) * (current_theta(i)
                                                    - seg_pts[2]) +
                                 seg_pts[0]).scalar();
            double ybar_local = (sin(seg_pts[3]) * (current_theta(i)
                                                    - seg_pts[2]) +
                                 seg_pts[1])
                    .scalar();
            lin_pts[i].setPosition((float) xbar_local, (float) ybar_local);
        }
    }
}

void
PathTracking::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for(const auto & c : lin_pts) {
        target.draw(c, states);
    }
}

