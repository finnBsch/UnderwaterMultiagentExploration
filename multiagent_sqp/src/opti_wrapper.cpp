//
// Created by finn on 6/20/23.
//
#include "opti_wrapper.h"
#include <eigen3/Eigen/Dense>
#include <stdexcept>

OptiWrapper::OptiWrapper(SQPParams *params): f(1,1), opti("conic") {
//    opti = new Opti("conic");
    this->params = params;
}

void OptiWrapper::initializeParams(int num_params) {
    parameter_vector = DM::zeros(num_params);
//    ca_parameter_vector = opti->parameter(num_params);
    ca_parameter_vector = opti.parameter(num_params);
    ca_states = reshape(ca_parameter_vector(Slice(0, params->N_horizon *
                                                 params->num_states *
                                                 params->num_agents)), params->num_states * params->num_agents, params->N_horizon);
    ca_u = reshape(ca_parameter_vector(Slice(params->N_horizon *
                                         params->num_states *
                                         params->num_agents, params->N_horizon * params->num_states *
                                                             params->num_agents + (params->N_horizon - 1) * params->num_inputs *
                                                                                  params->num_agents)), params->num_inputs * params->num_agents,
                      params->N_horizon - 1);
}

void OptiWrapper::initializeVariables(int num_variables) {
    variables = DM::zeros(num_variables);
//    ca_variables = opti->variable(num_variables);
    ca_variables = opti.variable(num_variables);
//    opti->subject_to(ca_variables < 1000);
    opti.subject_to(ca_variables < 1000);
//    opti->subject_to(ca_variables > -1000);
    opti.subject_to(ca_variables > -1000);
    delta_z = reshape(ca_variables(Slice(0, params->N_horizon *
                                                 params->num_states *
                                                 params->num_agents)), params->num_states * params->num_agents, params->N_horizon);
    delta_u = reshape(ca_variables(Slice(params->N_horizon *
                                         params->num_states *
                                         params->num_agents, params->N_horizon * params->num_states *
                                                             params->num_agents + (params->N_horizon - 1) * params->num_inputs *
                                                                                  params->num_agents)), params->num_inputs * params->num_agents,
                      params->N_horizon - 1);
}

void OptiWrapper::initializeSolver() {
//    opti->minimize(f);
    opti.minimize(f);
    costHess = Function("costHess", {ca_parameter_vector},
                        {hessian(f, ca_variables)});
    //    opti->solver("qpoases");
//    opti->solver("ipopt");
//    opti->solver("osqp", Dict({{"error_on_fail", false}}),
    opti.solver("osqp", Dict({{"error_on_fail", false}}),
                 Dict({{"max_iter", 600},
                       {"polish",true},
                       {"verbose", false},
                       {"eps_rel", 1e-3},
                       {"adaptive_rho_interval", 50}}));
    if(params->generate_code){
//        solve_compiled = opti->to_function("solve_compiled",
        solve_compiled = opti.to_function("solve_compiled",
                                           {ca_parameter_vector},
                                           {ca_variables});
        solve_compiled.generate("solve_compiled.cpp", Dict({{"verbose",
                                                             false},
                                                            {"cpp", true},
                                                            {"with_header",true}}));
//        std::cout << solve_compiled << "\n";
//        std::cout << ca_parameter_vector.dim() << "\n";
//        std::cout << ca_variables.dim() << "\n";
        throw "done generating code.";
    }
    if(params->use_generated_code){

        solve_compiled = external("solve_compiled");

    }
}

// Overloading the

void OptiWrapper::subject_to(const MX &expr) {
//    opti->subject_to(expr);
    opti.subject_to(expr);
}

void OptiWrapper::subject_to(const std::vector<MX> &g) {
//    opti->subject_to(g);
    opti.subject_to(g);
}

native_DM OptiWrapper::value(const MX &x, const std::vector<MX> &values) const {
//    return opti->value(x, values);
    return opti.value(x, values);
}

native_DM OptiWrapper::value(const DM &x, const std::vector<MX> &values) const {
//    return opti->value(x, values);
    return opti.value(x, values);
}

native_DM OptiWrapper::value(const SX &x, const std::vector<MX> &values) const {
//    return opti->value(x, values);
    return opti.value(x, values);
}

void OptiWrapper::solve() {
    if(params->use_generated_code) {
        auto sol = solve_compiled({parameter_vector});
        solution = sol[0];
    }
    else {
//        auto hess = costHess(parameter_vector)[0];
//        DM dense_hess;
//        dense_hess = densify(hess);
//        int rows = dense_hess.rows();
//        int columns = dense_hess.columns();
//        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> hess_eigen;
//        hess_eigen.resize(rows, columns);
//        std::memcpy(hess_eigen.data(), dense_hess.ptr(), sizeof
//                                                                                 (double) *
//                                                                         hess_eigen.size());
//        Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(hess_eigen);
//        std::cout << hess_eigen(0, 0) << "\n";
//        std::cout << parameter_vector <<"\n";
//        for(int i = 0; i < hess_eigen.rows(); ++i){
//
//            auto E = eigensolver.eigenvalues().col(0)[i];
//            std::cout << E << "," ;
//        }
//        std::cout << std::endl;
        set_value(); // Set the parameters
        try {
//            opti->solve_limited();
            opti.solve_limited();
//            auto status = opti->stats()["return_status"].to_string();
            auto status = opti.stats()["return_status"].to_string();
            if(!(status == "solved")) {
                if (status == "solved inaccurate") {
//                    std::cout << "Fine 8) \n";
                } else {
                    throw std::runtime_error(
                            "Solver failed with status " + status);
                }
            }
        }
        catch (CasadiException) {
//            std::cout << opti->stats()["return_status"] << std::endl;
//            std::cout << opti.stats()["return_status"] << std::endl;
        }
//        solution = opti->value(ca_variables);
        solution = opti.value(ca_variables);
    }
}

DM OptiWrapper::value(int idx0, int idx1) {
    return solution(Slice(idx0, idx1));
}

void OptiWrapper::set_value() {
//    opti->set_value(ca_parameter_vector, parameter_vector);
    opti.set_value(ca_parameter_vector, parameter_vector);
}

void OptiWrapper::update_parameter(int idx0, int idx1, const DM &value) {
    parameter_vector(Slice(idx0, idx1)) = value;
}



