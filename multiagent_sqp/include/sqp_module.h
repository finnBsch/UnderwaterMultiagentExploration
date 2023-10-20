//
// Created by finn on 5/26/23.
//

#ifndef CASADI_TEST_SQP_CONSTRAINT_H
#define CASADI_TEST_SQP_CONSTRAINT_H
#include <casadi/casadi.hpp>
#include <SFML/Graphics.hpp>
#include "opti_wrapper.h"

using namespace casadi;

class SQPModule : public sf::Drawable{
protected:
    int param_start_idx;
    int num_params;
public:
    virtual void
    initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                  MX &delta_u,
                  OptiWrapper *opti_obj, DM &states, MX &ca_param_vec, DM &param_vec) = 0;
    virtual int addToParams(int start_idx) noexcept = 0;
    virtual int addToVariables(int start_idx) noexcept = 0;
    virtual void iterationUpdate(OptiWrapper *opti_obj) = 0;
    virtual void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) = 0;
    virtual void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) = 0;
    virtual void draw(sf::RenderTarget& target, sf::RenderStates states)
    const = 0;
};


#endif //CASADI_TEST_SQP_CONSTRAINT_H
