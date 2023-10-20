
//
// Created by finn on 6/25/23.
//
#include "nearest_obstacle_avoidance.h"
#include "rrt_util.h"

void
NearestObstacleAvoidance::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                        MX &delta_z, MX &delta_u,
                                        OptiWrapper *opti_obj, DM &states,
                                        MX &ca_param_vec, DM &param_vec) {
    buildGradients();



    Slice state_slice(id * params->num_states, (id + 1) * params->num_states);
    Slice all;
    auto slack = opti_obj->ca_variables(Slice(sl_slack[0], sl_slack[1]));
    auto obs_x0 = ca_param_vec(Slice(sl_obs_x0[0], sl_obs_x0[1]));
    auto obs_x1 = ca_param_vec(Slice(sl_obs_x1[0], sl_obs_x1[1]));
    auto obs_y0 = ca_param_vec(Slice(sl_obs_y0[0], sl_obs_y0[1]));
    auto obs_y1 = ca_param_vec(Slice(sl_obs_y1[0], sl_obs_y1[1]));
    auto dist = ca_param_vec(Slice(sl_dist[0], sl_dist[1]));
    auto dir_t = ca_param_vec(Slice(sl_dir_t[0], sl_dir_t[1]));

    for(int i = 0; i < params->N_horizon; i++){
        auto J_map = grad_fn({{"i0", ca_states(state_slice, i)}, {"i1", obs_x0(i)}, {"i2", obs_y0(i)},
                              {"i3", obs_x1(i)}, {"i4", obs_y1(i)}});
        auto J = J_map["o0"];
        auto J_t_map = grad_to_path_fn({{"i0", ca_states(state_slice, i)}, {"i1", obs_x0(i)}, {"i2", obs_y0(i)},
                                    {"i3", obs_x1(i)}, {"i4", obs_y1(i)}});
        auto J_t = J_t_map["o0"];
        auto length_obs = (pow(obs_x1(i) - obs_x0(i), 2) + pow(obs_y1(i) - obs_y0(i), 2));

        auto t = ((ca_states(params->num_states * id , i) + delta_z(params->num_states * id , i) - obs_x0(i)) * (obs_x1(i) - obs_x0(i)) + (ca_states(params->num_states * id + 1, i ) + delta_z(params->num_states * id , i) - obs_y0(i)) * (obs_y1(i) - obs_y0(i)))/length_obs;
//        auto t_eval = ((ca_states(params->num_states * id , i)  - obs_x0(i)) * (obs_x1(i) - obs_x0(i)) + (ca_states(params->num_states * id + 1, i )- obs_y0(i)) * (obs_y1(i) - obs_y0(i)))/length_obs;
        auto val = if_else(0 <= t && t <= 1, 0.1, -100);
        auto delta = 0.0;
        auto val2 = if_else(-delta <= t && t <= 1.0 + delta, 1/(pow(dist(i), 2)+ 0.1),
                            0);
        opti_obj->subject_to(dist(i) + mtimes(J, delta_z(state_slice, i)) >=
        val);//  - pow(slack(i), 2));
        cost_fun += 400 * pow(slack(i), 2) + val2 * dir_t(i) * mtimes(J_t,
                                                                      delta_z
                                                                      (state_slice, i)) * 0;

    }
}

int NearestObstacleAvoidance::addToParams(int start_idx) noexcept {
    sl_obs_x0[0] = start_idx;
    sl_obs_x0[1] = start_idx + params->N_horizon;
    sl_obs_y0[0] = sl_obs_x0[1];
    sl_obs_y0[1] = sl_obs_y0[0] + params->N_horizon;
    sl_obs_x1[0] = sl_obs_y0[1];
    sl_obs_x1[1] = sl_obs_x1[0] + params->N_horizon;
    sl_obs_y1[0] = sl_obs_x1[1];
    sl_obs_y1[1] = sl_obs_y1[0] + params->N_horizon;
    sl_dist[0] = sl_obs_y1[1];
    sl_dist[1] = sl_dist[0] + params->N_horizon;
    sl_dir_t[0] = sl_dist[1];
    sl_dir_t[1] = sl_dir_t[0] + params->N_horizon;

    return 6 * params->N_horizon;
}

int NearestObstacleAvoidance::addToVariables(int start_idx) noexcept {
    sl_slack = {start_idx, start_idx + params->N_horizon};
    return params->N_horizon;
}

void NearestObstacleAvoidance::iterationUpdate(OptiWrapper *opti_obj) {

}

void NearestObstacleAvoidance::linearize(OptiWrapper *opti_obj, DM &states,
                                         DM &param_vec) {
    for(int i = 0; i < params->N_horizon; i++){
        auto x = states(id * params->num_states, i).scalar();
        auto y = states(id * params->num_states + 1, i).scalar();
        double min_dist = std::numeric_limits<double>::max();
        SQPObstacle* min_obs = nullptr;
        for(auto& obs: *obstacles){
            double t = ((x - obs->x0) * (obs->x1 - obs->x0) + (y - obs->y0) * (obs->y1 - obs->y0))/pow(obs->length, 2);
            if(0 <= t && t <= 1) {
                double px = obs->x0 + t * (obs->x1 - obs->x0);
                double py = obs->y0 + t * (obs->y1 - obs->y0);
//                std::cout << "px " << px << " py " << py << std::endl;
                double distance = (pow(px - x, 2) + pow(py - y, 2));
                if (distance < min_dist) {
                    min_dist = distance;
                    min_obs = obs;
                }
            }
        }
        if(!min_obs){
//            std::cout << "x " << x << " y " << y << " min_dist " << min_dist
//                      << std::endl;
            param_vec(sl_dist[0] + i) = 2000;
        }
        else {
//            std::cout << "x " << x << " y " << y << " min_dist " << min_dist
//                      << std::endl;
            param_vec(sl_obs_x0[0] + i) = min_obs->x0;
            param_vec(sl_obs_x1[0] + i) = min_obs->x1;
            param_vec(sl_obs_y0[0] + i) = min_obs->y0;
            param_vec(sl_obs_y1[0] + i) = min_obs->y1;
            param_vec(sl_dist[0] + i) = min_dist;
            min_obs->isactive = true;
            bool is_left = isLeft(min_obs->x0, min_obs->y0, min_obs->x1, min_obs->y1, x, y);
            if(is_left
            && !min_obs->left_block){
                param_vec(sl_dir_t[0] + i) = 0;
                min_obs->isactive = false;
            }
            else if(!is_left &&
            min_obs->left_block){
                param_vec(sl_dir_t[0] + i) = 0;
                min_obs->isactive = false;
            }
            else {
                param_vec(sl_dir_t[0] + i) = min_obs->dir_t;
            }
        }
    }

}

void NearestObstacleAvoidance::warmStart(OptiWrapper *opti_obj, DM &states,
                                         DM &param_vec) {

}

void NearestObstacleAvoidance::draw(sf::RenderTarget &target,
                                    sf::RenderStates states) const {

}

void NearestObstacleAvoidance::buildGradients() {
    MX states = MX::sym("sts", params->num_states);
    auto x = states(0);
    auto y = states(1);

    MX obs_x0 = MX::sym("obs_x0", 1);
    MX obs_y0 = MX::sym("obs_y0", 1);
    MX obs_x1 = MX::sym("obs_x1", 1);
    MX obs_y1 = MX::sym("obs_y1", 1);

    auto length_obs = pow(obs_x1 - obs_x0, 2) + pow(obs_y1 - obs_y0, 2);
    auto t = ((x - obs_x0) * (obs_x1 - obs_x0) + (y - obs_y0) * (obs_y1 - obs_y0))/length_obs;
    auto px = obs_x0 + t * (obs_x1 - obs_x0);
    auto py = obs_y0 + t * (obs_y1 - obs_y0);

    auto distance = (pow(px - x, 2) + pow(py - y, 2));

    auto d_distance_dx = jacobian(distance, states);

    dist_fn = Function("dist", {states, obs_x0, obs_y0, obs_x1, obs_y1}, {distance});
    grad_fn = Function("grad", {states, obs_x0, obs_y0, obs_x1, obs_y1}, {d_distance_dx});
    grad_to_path_fn = Function("grad_to_path", {states, obs_x0, obs_y0, obs_x1, obs_y1},
                               {jacobian(t, states)});
}

NearestObstacleAvoidance::NearestObstacleAvoidance(SQPParams *params, int id,
                                                   std::vector<SQPObstacle *> *obstacles):
                                                   params(params),
                                                   id(id),
                                                   obstacles(obstacles){

}
