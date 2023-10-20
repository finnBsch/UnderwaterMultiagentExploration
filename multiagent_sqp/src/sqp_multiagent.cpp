#include "../include/sqp_multiagent.h"
#include "chrono"
#include <limits.h>
using namespace std::chrono;

// TODO make all modules sqp_modules and store them in a vector while keeping
//  vectors for modules of specific class

SQPMultiAgent::SQPMultiAgent(SQPParams params_, SplinePath *path,
                             InfoField *info_field,
                             Eigen::Matrix<double, Eigen::Dynamic, 1>
                                     &initial_state,
                             std::vector<SQPObstacle*>* obstacles)
        : params(params_), path(path), info_field(info_field)
                        , opti(&params)
                        {
    if(params.viz_all){
        for(int i = 0; i < params.num_agents; i++){
            plans.emplace_back(sf::LineStrip, params.N_horizon);
            for(int j = 0; j < params.N_horizon; j++){
                plans[i][j] = sf::Vertex(sf::Vector2f(0.0f , 0.0f),
                                         sf::Color::White);
            }
        }
    }
    E::Matrix<double, E::Dynamic, 1> states;

    // Slicers
    sl_states = {0, params.num_agents * params.num_states * params.N_horizon};

    sl_u = {sl_states[1], sl_states[1] + params.num_agents * params.num_inputs * (params.N_horizon - 1)};

    states.resize(params.num_agents * params.num_states);
    states.setZero();
    states_vals = DM::zeros(params.num_agents*params.num_states,
                               params.N_horizon);
    eigen_states.resize(params.num_agents*params.num_states, params.N_horizon);
    eigen_states.col(0) = initial_state;
    std::memcpy(states_vals.ptr(), eigen_states.data(), sizeof(double)
                                                        *eigen_states.size());
    u_vals = DM::zeros(params.num_agents*params.num_inputs, params.N_horizon - 1);



    if(params.information_cost){
        for(int i = 0; i < params.num_agents; i++) {
            sqp_modules.push_back(new InformationCost(&params, i, info_field));
        }
    }
    // Add connectivity constraints
    for(int i = 1; i < params.N_horizon; i++){
        if(params.communication_constraint && params.num_agents > 1) {
            connectivity_constraints.push_back(
                    new ConnectivityConstraint(&params,
                                               0, i));
            sqp_modules.push_back(connectivity_constraints[i - 1]);
        }
        if(params.agent_collision_constraint && params.num_agents > 1) {
            collision_constraints.push_back(new
            CollisionConstraint(i, &params));
            sqp_modules.push_back(collision_constraints[i - 1]);
        }
    }
    // Distribution!!
    if(params.distribution_cost && params.num_agents > 1) {
        sqp_modules.push_back(
                new DistributionCost(&params, path));
    }
    if(params.corridor_constraint) {
        sqp_modules.push_back(new CorridorConstraint(&params, path));
    }
    // Add the actual agents
    for(int i = 0; i < params.num_agents; i++){
        agents.push_back(createAgent(&params, i));
        sqp_modules.push_back(agents[i]);
        if(params.nearest_obstacle_constraint){
            sqp_modules.push_back(new NearestObstacleAvoidance(&params, i, obstacles));
        }
    }
//    if(params.path_tracking){
//        sqp_modules.push_back(new PathTracking(&params, 0));
//    }
    if(params.path_tracking){
        path_tracking = new COMPathTracking(&params, path);
        sqp_modules.push_back(path_tracking);
    }
    int n_params = params.num_agents * params.N_horizon * params.num_states
            + params.num_agents * params.num_inputs * (params.N_horizon - 1);
    int n_variables = n_params;
    for(auto& sqpmod : sqp_modules){
        n_params += sqpmod->addToParams(n_params);
        n_variables += sqpmod->addToVariables(n_variables);
//        std::cout << n_params << "\n";
    }
    opti.initializeParams(n_params);
    opti.initializeVariables(n_variables);

    for(auto& sqpmod : sqp_modules){
        sqpmod->initializeSQP(opti.f, opti.ca_states, opti.ca_u, opti.delta_z,
                              opti.delta_u,
                              &opti,
                              states_vals, opti.ca_parameter_vector,
                              opti.parameter_vector);
    }

//    for(auto& cc: connectivity_constraints){
//        cc.initializeSQP(f, ca_states, ca_u, delta_z, delta_u, &opti,
//                         states_vals);
//    }
//    for(auto& colc: collision_constraints){
//        colc.initializeSQP(f, ca_states, ca_u, delta_z, delta_u, &opti,
//                           states_vals);
//    }
    opti.initializeSolver();

//    auto Func = opti.to_function("solver", {ca_parameter_vector, ca_states, ca_u}, {delta_z, delta_u}, Dict({{"compiler", "shell"},
//                                                                         {"jit", true},
//                                                                         {"jit_options"
//                                                                          ".compiler", "gcc"},
//                                                                         {"jit_options"
//                                                                          ".flags", "-O3"}}));
//    std::cout << Func <<"XXXXXXXXXX\n";

    //    Func(ca_parameter_vector);
//    auto x = opti.x();
//    auto p = opti.p();
//    auto f_h = hessian(opti.f(), opti.x());
//    f_Hess = Function("f_Hess", {opti.x(), opti.p()}, {f_h});

}

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &SQPMultiAgent::getSol() {
    return eigen_states;
}

void SQPMultiAgent::solve(const Eigen::Matrix<double, Eigen::Dynamic, 1>& init_states, bool reset) {
    auto start = high_resolution_clock::now();
    if(!reset){
        // Reuse linearization point. TODO Extrapolate
        // TODO Reuse previous u..
        eigen_states(Eigen::all, Eigen::seq(1, params.N_horizon - 2)) =
                eigen_states(Eigen::all, Eigen::seq(2, params.N_horizon - 1));
        eigen_states(Eigen::all, params.N_horizon - 1) = eigen_states.col(params.N_horizon - 2) *2 -  eigen_states.col(params.N_horizon - 3);
        eigen_states( Eigen::seq(0, params.num_agents * params.num_states - 1), 0) = init_states;
    }
    else{
        for(int j = 0; j < params.N_horizon; j++){
            eigen_states(Eigen::seq(0, params.num_agents * params.num_states - 1), j) = init_states;
        }
        opti.update_parameter(sl_u[0], sl_u[1], 0);
    }
//    for(int i = 0; i < params.num_agents; i++){
//        agents[i].initializeNextSol(&opti, init_states);
//    }
    std::memcpy(states_vals.ptr(), eigen_states.data(), sizeof(double)
                                                        *eigen_states.size());
    if(!reset) {
        for (auto &sqpmod: sqp_modules) {
            sqpmod->warmStart(&opti, states_vals, opti.parameter_vector);
        }
    }
    for(auto & sqpmod : sqp_modules){
        sqpmod->linearize(&opti, states_vals, opti.parameter_vector);
    }
    opti.update_parameter(sl_states[0],
                          sl_states[1], reshape(states_vals,
                                                params.num_agents * params.num_states *
                                                params.N_horizon, 1));
    opti.solve();
    for(int i = 0; i < params.sqp_iters; i++){
        // Not sure if opti.value on slices works!! Maybe better do this
        // using slices?
        auto du = reshape(opti.value(sl_u[0], sl_u[1]), params.num_agents *
        params
        .num_inputs, params.N_horizon - 1);
        auto dz = reshape(opti.value(sl_states[0], sl_states[1]), params
        .num_agents * params
        .num_states, params.N_horizon);
        u_vals = u_vals + params.alpha_damping*du;
        states_vals = states_vals + params.alpha_damping*dz;
        for(auto& sqpmod : sqp_modules){
            sqpmod->iterationUpdate(&opti);
        }
        for(auto & sqpmod : sqp_modules){
            sqpmod->linearize(&opti, states_vals, opti.parameter_vector);
        }
        opti.update_parameter(sl_states[0], sl_states[1], reshape
        (states_vals, params.num_agents * params.num_states * params.N_horizon, 1));
        opti.update_parameter(sl_u[0], sl_u[1], reshape(u_vals, params
        .num_agents * params.num_inputs * (params.N_horizon - 1), 1));
//        opti.set_initial(opti.delta_u, du);
//        opti.set_initial(opti.delta_z, dz);
        opti.solve();
    }
//    auto res = opti.advanced().res();
//    auto x_ = res["x"];
//    auto p_ = res["p"];
//    auto hess_m = f_Hess(DMDict{{"i0", x_},{"i1", p_}});
//    auto hess = hess_m["o0"];
//    std::cout << hess <<"\n";
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
//    std::cout << "Took " << duration.count() << "ms.\n";
    auto du = reshape(opti.value(sl_u[0], sl_u[1]), params.num_agents *
                                                    params
                                                            .num_inputs, params.N_horizon - 1);
    auto dz = reshape(opti.value(sl_states[0], sl_states[1]), params
                                                                      .num_agents * params
                                                                      .num_states, params.N_horizon);
    u_vals = u_vals + params.alpha_damping*du;
    states_vals = states_vals + params.alpha_damping*dz;
    std::memcpy( eigen_states.data(), states_vals.ptr(), sizeof(double)
                                                         *eigen_states.size());
    if(params.viz_all){
        for(int i = 0; i < params.num_agents; i++){
            for(int j = 0; j < params.N_horizon; j++){
                plans[i][j].position = sf::Vector2f(states_vals(i *
                        params.num_states, j).scalar(), states_vals(i *
                        params.num_states +
                        1, j).scalar());
            }
        }
    }
}

void
SQPMultiAgent::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    if(params.viz_all){
        for(const auto& plan : plans) {
            target.draw(plan, states);
        }
        for(const auto & sqp_mod : sqp_modules){
            sqp_mod->draw(target, states);
        }
    }
}

double SQPMultiAgent::getTheta() const {
    return path_tracking->getTheta();
}

double SQPMultiAgent::getMinTheta() const {
    double min_theta = MAXFLOAT;
    for(auto& corr : corridor_constraints){
        auto theta = corr->getTheta();
        if(theta < min_theta){
            min_theta = theta;
        }
    }
    return min_theta;
}

