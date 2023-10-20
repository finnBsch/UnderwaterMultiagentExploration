//
// Created by finn on 5/26/23.
//
#include "../include/connectivity_constraint.h"
using namespace std::chrono;
ConnectivityConstraint::ConnectivityConstraint(SQPParams* params, int
num_lagrange, int time_id) {
    this->num_lagrange = num_lagrange;
    this->params = params;
    this->time_id = time_id;
    this->sig_sq = -pow(params->connectivity.r_c, 2)/(2 * log
            (params->connectivity.fiedler_eps));
    A.resize(params->num_agents, params->num_agents);
    A.setZero();
    dAdxi.resize(params->num_agents, params->num_agents, 2);
    D.resize(params->num_agents, params->num_agents);
    D.setZero();
    grad.resize(params->num_agents * params->num_states, 1);
    grad.setZero();
    fiedlerHess.resize(params->num_agents * params->num_states + num_lagrange,
                       params->num_agents *
                                                               params->num_states +
                                                               num_lagrange);
    fiedlerHess.setZero();
//    dvdL.resize(params->num_agents, params->num_agents);
    dvdState.resize(params->num_agents * 2, params->num_agents);
    states.resize(params->num_agents*params->num_states);
//    dvdL.setZero();
//    buildCasadiGradients();
}

void
ConnectivityConstraint::updateMatrix(DM& states) {
    // TODO: Potential huge speedups possible..
    std::memcpy(this->states.data(), states(Slice(), time_id).ptr(), sizeof
    (double) *
                                                      this->states.size());
//    std::cout << this->states << "\n";
    for(int i = 0; i < params->num_agents; i += 1){
        for(int j = 0; j < params->num_agents; j += 1){
            A(i, j) = connectivityFunction(i, j);
            auto deriv = connectivityFunctionDeriv(i, j);
            dAdxi(i, j, 0) = deriv[0];
            dAdxi(i, j, 1) = deriv[1];
        }
    }
    for(int i = 0; i < params->num_agents; i += 1){
        D(i, i) = A.row(i).sum();
    }
    auto L = D - A;
    E::SelfAdjointEigenSolver<E::MatrixXd> es(L);
    auto eigs = es.eigenvalues();
    // TODO: What happens if the smallest occurs more than once? Need to
    //  catch that.
    fiedler_eig = eigs[1];
    auto fiedler_v = es.eigenvectors().col(1);
//    Slice all;  // Basic slicer
//    auto res2 = Hess(DMDict({{"i0", states_vec}}));
//    std::cout << res2["o0"] <<"\n";
//
//    E::MatrixXd eig_mat(params->num_agents, params->num_agents);
//    auto start = high_resolution_clock::now();
//    auto res = Jac(DMDict({{"i0", states_vec}}));
//    auto mat = res["o0"];
//    for(int i = 0; i < params->num_agents*params->num_states; i++){
//        for(int k = 0; k < params->num_agents; k++){
//            for(int l = 0; l < params->num_agents; l++){
//                eig_mat(l, k) = mat(l + k * params->num_agents, i).scalar();
//            }
//        }
////        std::cout << eig_mat << "\n";
////        std::cout << mat(all, i) <<"\n";
//        grad(i, 0) = fiedler_v.transpose() * eig_mat * fiedler_v;
//    }
//    auto stop = high_resolution_clock::now();
//    auto duration = duration_cast<microseconds>(stop - start);
//
//
//    std::cout << duration.count() << std::endl;
//    std::cout << grad << "\n\n";
//    start = high_resolution_clock::now();
    fillLGrad(fiedler_v);
//    filldvdState(fiedler_v, L, fiedler_eig);
//    fillFiedlerHess(fiedler_v);
//    stop = high_resolution_clock::now();
//    duration = duration_cast<microseconds>(stop - start);
//
//    std::cout << duration.count() << std::endl;
//    std::cout << fiedlerHess << std::endl;
//    std::cout << grad << std::endl;

}

void ConnectivityConstraint::buildCasadiGradients() {
    Slice all;  // Basic slicer
    // Gradients of L
    MX ca_states = MX::sym("ca_states", params->num_agents*params->num_states,
                           1);
    MX ca_A = MX::sym("ca_A", params->num_agents, params->num_agents);
    for(int i = 0; i < params->num_agents; i++){
        for(int j = 0; j < params->num_agents; j++) {
            auto x0 = ca_states(i * params->num_states, 0);
            auto y0 = ca_states(i * params->num_states + 1, 0);
            auto x1 = ca_states(j * params->num_states, 0);
            auto y1 = ca_states(j * params->num_states + 1, 0);
            auto sq_norm = pow(x0 - x1, 2) + pow(y0 - y1, 2);
            ca_A(i, j) = exp(-sq_norm/(2.0 * sig_sq));
        }
    }
    MX ca_D = MX::sym("ca_D", params->num_agents, params->num_agents);
    for(int i = 0; i < params->num_agents; i++){
        for(int j = 0; j < params->num_agents; j++){
            if(i != j){
                ca_D(i, j) = 0.0;
            }
            else{
                ca_D(i, i) = sum2(ca_A(i, all));
            }
        }
    }
    MX ca_L = ca_D - ca_A;
//    std::cout << ca_L.dim() << "\n";
    auto jac = simplify(jacobian(ca_L, ca_states));
    auto hess = jacobian(jac, ca_states);
    Jac = Function("jac_con", {ca_states}, {jac}, Dict({{"compiler", "shell"},
                                                        {"jit", params->jit},
                                                        {"jit_options"
                                                         ".compiler", "gcc"},
                                                        {"jit_options"
                                                         ".flags", "-O3"}}));
    Hess = Function("hess_con", {ca_states}, {hess}, Dict({{"compiler", "shell"},
                                                           {"jit", params->jit},
                                                           {"jit_options"
                                                            ".compiler", "gcc"},
                                                           {"jit_options"
                                                            ".flags", "-O3"}}));
    // Gradients of fiedler-vector w.r.t. the states

}

double
ConnectivityConstraint::connectivityFunction(int id0, int id1) {
    if (id0 == id1) {
        return 1.0;
    }
    double x0 = states(id0 * params->num_states);
    double y0 = states(id0 * params->num_states + 1);
    double x1 = states(id1 * params->num_states);
    double y1 = states(id1 * params->num_states + 1);
    double sq_norm = pow(x0 - x1, 2) + pow(y0 - y1, 2);
    return exp(-sq_norm / (2.0 * sig_sq));
}

/**
 * Calculated the derivative with respect to the first argument.
 * @param id0
 * @param id1
 * @return
 */
std::array<double, 4> ConnectivityConstraint::connectivityFunctionDeriv(int
                                                                        id0, int id1) {
    if(id0 == id1){
        return std::array<double, 4>{0.0, 0.0, 0.0, 0.0};
    }
    double x0 = states(id0 * params->num_states);
    double y0 = states(id0 * params->num_states + 1);
    double x1 = states(id1 * params->num_states);
    double y1 = states(id1 * params->num_states + 1);
    return std::array<double, 4>{
            A(id0, id1) * 0.5*(-2*x0 + 2*x1)/sig_sq,
            A(id0, id1) * 0.5 *(-2*y0 + 2*y1)/sig_sq,
            A(id0, id1) * 0.5*(2*x0 - 2*x1)/sig_sq,
            A(id0, id1) * 0.5*(2*y0 - 2*y1)/sig_sq
    };
}


void ConnectivityConstraint::fillLGrad(const E::Matrix<double, E::Dynamic,
        1>& v2) {
    grad.setZero();
    for(int i = 0; i < params->num_agents; i++) {
        for(int j = 0; j < params->num_agents; j++) {
            grad(i * params->num_states, 0) += dAdxi(i, j, 0) * pow(v2(i, 0) -
                    v2(j,
                                                                          0),2);
            grad(i * params->num_states + 1, 0) += dAdxi(i, j, 1) * pow(v2(i,
                                                                           0)
                    - v2
                    (j, 0), 2);
        }
    }
}


void
ConnectivityConstraint::fillFiedlerHess(const E::Matrix<double, E::Dynamic, 1> &v2) {
    fiedlerHess.setZero();
    Eigen::Matrix<double, 2, 2> hessian;
    for(int i = 0; i < params->num_agents; i++){
        for(int j = 0; j < params->num_agents; j++){
            // v2^T * d^2 L/(dxi dxj) * v2
            if (i == j){
                // This calculation is obviously not required, just here for
                // clarity
                fiedlerHess(i * params->num_states, j*params->num_states) +=
                        0.0;
                fiedlerHess(i * params->num_states + 1, j*params->num_states)
                +=
                        0.0;
                fiedlerHess(i * params->num_states, j*params->num_states + 1)
                +=
                        0.0;
                fiedlerHess(i * params->num_states + 1, j*params->num_states
                + 1)
                +=
                        0.0;
            }
            else{
                double diff = pow(v2(i, 0) - v2(j, 0), 2);
                connectivityFunctionHess(i, j, hessian);
                fiedlerHess(i + params->num_states, j * params->num_states) +=
                        hessian(0, 0) * diff;
                fiedlerHess(i + params->num_states + 1, j *
                params->num_states) +=
                        hessian(1, 0) * diff;
                fiedlerHess(i + params->num_states, j * params->num_states +
                1) +=
                        hessian(0, 1) * diff;
                fiedlerHess(i + params->num_states + 1, j *
                params->num_states +
                1) +=
                        hessian(1, 1) * diff;
            }
            // dv2/dxj^T * dL/dxi * v2 = (v2^T * dL/dxi * dv2/dxj)^T, but
            // it's scalar, I think, so it should be the same

            // dv2/dxj^T is of shape (params->num_agents, 2). The second
            // dimension
            // corresponds to the two coordinates of the agent j. In the same
            // manner, dL/dxi has an extra dimension for the two coordinates
            // for agent i. This gives the 4 entries of the local hessian.

            for(int k = 0; k < params->num_agents; k++){
                // build the scalar products
                if (k != i){
                    // d/dxixj
                    fiedlerHess(i * params->num_states, j *
                    params->num_states) +=
                            dAdxi(i, k, 0) * (v2(k, 0) - v2(i, 0)) * dvdState
                                    (j*2 , k) * 2;
                    // d/dxiyj
                    fiedlerHess(i * params->num_states, j *
                    params->num_states +
                    1) +=
                            dAdxi(i, k, 0) * (v2(k, 0) - v2(i, 0)) * dvdState
                                    (j*2 + 1 , k) * 2;
                    // d/dyixj
                    fiedlerHess(i * params->num_states + 1, j *
                    params->num_states) +=
                            dAdxi(i, k, 1) * (v2(k, 0) - v2(i, 0)) * dvdState
                                    (j*2 , k) * 2;
                    // d/dyiyj
                    fiedlerHess(i * params->num_states + 1, j *
                    params->num_states +
                    1) +=
                            dAdxi(i, k, 1) * (v2(k, 0) - v2(i, 0)) * dvdState
                                    (j*2 + 1 , k) * 2;
                }
                else{
                    // catch the special case. not dAik/dxi * (vk - vi), but
                    // -dAi:/dxi * v (scalar product)
                    // + (sum_over_m ( dAim/dxi) ) * v_i.
                    double A_sum_xi = 0;  // dAim/dxi sum
                    double A_sum_yi = 0;  // dAim/dxi sum
                    double sum_xi = 0;
                    double sum_yi = 0;
                    for(int m = 0; m < params->num_agents; m++){
                        sum_xi -= dAdxi(i, m, 0) * v2(m, 0);
                        sum_yi -= dAdxi(i, m, 1) * v2(m, 0);
                        A_sum_xi += dAdxi(i, m, 0);
                        A_sum_yi += dAdxi(i, m, 1);
                    }
                    double val_xi = sum_xi + A_sum_xi * v2(i, 0);
                    double val_yi = sum_yi + A_sum_yi * v2(i, 0);
                    // d/dxixj
                    fiedlerHess(i * params->num_states, j *
                    params->num_states) +=
                            val_xi * dvdState(j * 2, k) * 2;
                    // d/dxiyj
                    fiedlerHess(i * params->num_states, j *
                    params->num_states +
                    1) +=
                            val_xi * dvdState(j * 2 + 1, k) * 2;
                    // d/dyixj
                    fiedlerHess(i * params->num_states + 1, j *
                    params->num_states) +=
                            val_yi * dvdState(j * 2, k) * 2;
                    // d/dyiyj
                    fiedlerHess(i * params->num_states + 1, j *
                    params->num_states +
                    1) +=
                            val_yi * dvdState(j * 2 + 1, k)  * 2;
                }
            }
            // v2^T * dL/dxi * dv2/dxj is included in the "*2" above.
        }
    }
    fiedlerHess *= - states(params->num_agents*params->num_states, 0); //
    // multiply
    // with
    // negative lagrange multiplier.
    // Fill the columns for dLambda/dxdlagrange
    fiedlerHess(E::seq(0, params->num_agents * params->num_states - 1),
                params->num_agents*params->num_states)
            = - grad;
    fiedlerHess(params->num_agents*params->num_states, E::seq(0,
                                                         params->num_agents *
    params->num_states -
    1))
    = -
            grad
            .transpose();
}

void ConnectivityConstraint::connectivityFunctionHess(int id0, int id1,
                                                      E::Matrix<double, 2, 2>&
                                                      ref) {
    // x0, y0, x1, y1
    if(id0 == id1){
        ref.setZero();
    }
    else{
        double Aij = A(id0, id1);
        double x0 = states(id0 * params->num_states);
        double y0 = states(id0 * params->num_states + 1);
        double x1 = states(id1 * params->num_states);
        double y1 = states(id1 * params->num_states + 1);
        ref(0, 0) = Aij * (sig_sq - pow(x0,2) + 2.0*x0*x1
                           - pow(x1,2))/pow(sig_sq,2);
        ref(0, 1) = Aij * (-(x0 - x1)*(y0 - y1)/pow(sig_sq,2));
        ref(1, 0) = Aij * (-(x0 - x1)*(y0 - y1)/pow(sig_sq,2));
        ref(1, 1) = Aij * ((sig_sq - pow(y0,2) + 2.0*y0*y1
                            - pow(y1,2))/pow(sig_sq,2));
    }
}

void
ConnectivityConstraint::filldvdState(const E::Matrix<double, E::Dynamic, 1>& v2, const
E::Matrix<double, E::Dynamic, E::Dynamic>& L, const double&
eigv) {
    // dv/dL (i, j) is the i'th entry of dv differentiated by L(i, j).
    // all the other entries of dv/dL(i,j) are zero, hence one value is enough
    // to be stored. This is important to remember when assembling for
    // further calculation!
    E::MatrixXd identity = E::MatrixXd::Identity(params->num_agents,
                                                 params->num_agents);
    pseudo_inverse = (eigv * identity - L).completeOrthogonalDecomposition()
            .pseudoInverse();
    // now we build the dv/dxi = dv/dL * dL/dxi
    dvdState.setZero();
    for(int i = 0; i < params->num_agents; i++){
        for(int j = 0; j < params->num_agents; j++){
            auto deriv = connectivityFunctionDeriv(i, j);  // dAij/dxi
            // diagonal
            dAdxi(i, j, 0) = dAdxi(i, j, 0);
            dAdxi(i, j, 1) = dAdxi(i, j, 1);
            dvdState(i * 2, j) += dAdxi(i, j, 0) * dvdL(v2,  j, j);
            dvdState(i * 2 + 1, j) += dAdxi(i, j, 1) * dvdL(v2,  j, j);
            // term at cross
            dvdState(i * 2, i) += dAdxi(i, j, 0) * dvdL(v2,  i, i);
            dvdState(i * 2 + 1, i) += dAdxi(i, j, 1) * dvdL(v2,  i, i);
            // column
            dvdState(i * 2, j) += -dAdxi(i, j, 0) * dvdL(v2,  j, i);
            dvdState(i * 2 + 1, j) += - dAdxi(i, j, 1) * dvdL(v2,  j, i);
            // row
            dvdState(i * 2, i) += -dAdxi(i, j, 0) * dvdL(v2,  i, j);
            dvdState(i * 2 + 1, i) += -dAdxi(i, j, 1) * dvdL(v2,  i, j);
        }
        // Technically, a correction term would be required because the row
        // and column two times sum up (-dAii/dxi). However since dAii/dxi =
        // 0, we don't need to do that
    }
}

double ConnectivityConstraint::dvdL(const E::Matrix<double, E::Dynamic, 1> &v2,
                                    int i, int j) {
    return pseudo_inverse(i, j) * v2(j, 0);
}

void
ConnectivityConstraint::linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) {
    updateMatrix(states);
    // TODO: Speedup: Make DM sparse, copy only the correct values.
//    DM Hess_param_vals = DM::zeros(params->num_states*params->num_agents +
//            num_lagrange,
//                                   params->num_states*params->num_agents +
//                                   num_lagrange);
    DM grad_param_vals = DM::zeros(params->num_states*params->num_agents, 1);
//    std::memcpy(Hess_param_vals.ptr(), fiedlerHess.data(), sizeof(double)
//                                                           *fiedlerHess.size());
    std::memcpy(grad_param_vals.ptr(), grad.data(), sizeof(double)
                                                    *grad.size());
    param_vec(Slice(sl_grad[0], sl_grad[1])) = grad_param_vals;
    param_vec(Slice(sl_fiedler[0], sl_fiedler[1])) = fiedler_eig;

}



void
ConnectivityConstraint::initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u,
                                      MX &delta_z, MX &delta_u,
                                      OptiWrapper *opti_obj, DM &states,
                                      MX &ca_param_vec, DM &param_vec) {
    auto slack = opti_obj->ca_variables(Slice(sl_slack[0], sl_slack[1]));

    auto grad_param = ca_param_vec(Slice(sl_grad[0], sl_grad[1]));


    auto fiedler_eig_param = ca_param_vec(Slice(sl_fiedler[0], sl_fiedler[1]));

    MX c = -fiedler_eig_param + params->connectivity.fiedler_eps - mtimes(grad_param.T()
            , delta_z(Slice(),time_id));
    opti_obj->subject_to(c <= pow(slack, 2));
    cost_fun += pow(slack, 2) * params->connectivity.slack_weight;
    Slice all;
//    cost_fun -= 0.5*mtimes(transpose(delta_z(all, time_id)),mtimes(
//                                  Hess_param,
//                           delta_z(all, time_id)));
}

bool ConnectivityConstraint::checkFeasible(Opti *opti_obj, DM &states) {
    return false;
}

int ConnectivityConstraint::addToParams(int start_idx) noexcept {
    param_start_idx = start_idx;
    sl_fiedler[0] = start_idx;
    sl_fiedler[1] = start_idx + 1;
    sl_grad[0] = start_idx + 1;
    sl_grad[1] = start_idx + 1 + params->num_agents * params->num_states;
    num_params = 1 // fiedler_eig_param
            + params->num_agents * params->num_states; // grad
    return num_params;
}

int ConnectivityConstraint::addToVariables(int start_idx) noexcept {
    sl_slack[0] = start_idx;
    sl_slack[1] = start_idx + 1;

    return 1;
}




