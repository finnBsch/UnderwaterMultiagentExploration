//
// Created by finn on 6/16/23.
//

#ifndef MASTERMAIN_DYNAMICS_H
#define MASTERMAIN_DYNAMICS_H

#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include "configs.h"
#include "sqp_module.h"

class Dynamics : public SQPModule {
protected:
    int id;
public:
};

#endif //MASTERMAIN_DYNAMICS_H
