//
// Created by finn on 6/20/23.
//

#ifndef MASTERMAIN_OPTI_WRAPPER_H
#define MASTERMAIN_OPTI_WRAPPER_H

#include <casadi/casadi.hpp>
#include "configs.h"

using namespace casadi;

class OptiWrapper {
private:
    Opti opti;
    Function solve_compiled;
    SQPParams* params;

    Function costHess;

public:
    OptiWrapper(SQPParams* params);
    void initializeParams(int num_params);
    void initializeVariables(int num_variables);
    void initializeSolver();

    void subject_to(const MX& expr);
    void subject_to(const std::vector<MX>& g);

    void solve();

    DM value(int idx0, int idx1);

    native_DM value(const MX& x, const std::vector<MX>& values=std::vector<MX>()) const;
    native_DM value(const DM& x, const std::vector<MX>& values=std::vector<MX>()) const;
    native_DM value(const SX& x, const std::vector<MX>& values=std::vector<MX>()) const;

    void set_value();

    void update_parameter(int idx0, int idx1, const DM& value);

    MX f;
    MX ca_parameter_vector;  // All parameters
    MX ca_variables;  // All variables

    MX ca_states; // All states as parameters!
    MX ca_u;  // All control inputs as parameters!

    MX delta_z;  // All delta states as variables!
    MX delta_u;  // All delta control inputs as variables!

    // Cannot split these because slicing/reshaping is copying contents for DM!
    DM parameter_vector;
    DM variables;

    DM solution;


};

#endif //MASTERMAIN_OPTI_WRAPPER_H
