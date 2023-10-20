//
// Created by finn on 5/26/23.
//

#ifndef CASADI_TEST_CONNECTIVITY_CONSTRAINT_H
#define CASADI_TEST_CONNECTIVITY_CONSTRAINT_H
#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "bicycle_dynamics.h"
#include "configs.h"
#include "sqp_module.h"

using namespace casadi;
namespace E = Eigen;

/**
 * Connectivity Constraint for one timestep. The SQP should contain a vector
 * of these for each time step. Furthermore, the corresponding timestep needs
 * to be supplied to make sure the correct state information is used.
 */
class ConnectivityConstraint: public SQPModule {
private:
    SQPParams* params;
    int time_id;
    int num_lagrange;
    double sig_sq;
    Function Jac;
    Function Hess;
    Function L;
//    SubIndex<MX, int> slack;
//    SubIndex<MX, int> fiedler_eig_param;
//    SubIndex<MX, Slice> grad_param;
    // Params
    std::array<int, 2> sl_fiedler;
    std::array<int, 2> sl_grad;
    // Variables
    std::array<int, 2> sl_slack;
//    E::Matrix<double, E::Dynamic, E::Dynamic> dvdL;
    E::Matrix<double, E::Dynamic, E::Dynamic> dvdState;
    double fiedler_eig = 0.0;
    E::Matrix<double, E::Dynamic, E::Dynamic> A;
    E::Tensor<double, 3> dAdxi;
    E::Matrix<double, E::Dynamic, E::Dynamic> D;
    E::Matrix<double, E::Dynamic, E::Dynamic> pseudo_inverse;
    E::Matrix<double, E::Dynamic, 1> grad;
    E::Matrix<double, E::Dynamic, E::Dynamic> fiedlerHess;
    E::Matrix<double, E::Dynamic, 1> states;
    double connectivityFunction(int id0, int id1);
    std::array<double, 4> connectivityFunctionDeriv(int id0, int id1);
    void connectivityFunctionHess(int id0, int id1, E::Matrix<double, 2,
            2>& ref);
    void fillLGrad(const E::Matrix<double, E::Dynamic, 1>& v2);
    void fillFiedlerHess(const E::Matrix<double, E::Dynamic, 1>& v2);
    void filldvdState(const E::Matrix<double, E::Dynamic, 1>& v2, const
    E::Matrix<double, E::Dynamic, E::Dynamic>& L, const double&
    eigv);
    double dvdL(const E::Matrix<double, E::Dynamic, 1>& v2, int i, int j);
    void buildCasadiGradients();
public:
    int addToParams(int start_idx) noexcept;
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u,
                       OptiWrapper *opti_obj, DM &states, MX &ca_param_vec,
                       DM &param_vec);
    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec);
    void updateMatrix(DM &states);
    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) {}
    ConnectivityConstraint(SQPParams* params, int
    num_lagrange, int time_id);
    void iterationUpdate(OptiWrapper *opti_obj){}
    void draw(sf::RenderTarget& target, sf::RenderStates states) const{}
    bool checkFeasible(Opti *opti_obj, DM &states);

    int addToVariables(int start_idx) noexcept override;
};

#endif //CASADI_TEST_CONNECTIVITY_CONSTRAINT_H
