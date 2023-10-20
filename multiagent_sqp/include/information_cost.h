//
// Created by finn on 5/29/23.
//

#ifndef SQP_INFORMATION_COST_H
#define SQP_INFORMATION_COST_H

#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include "info_field.h"
#include "sqp_module.h"
#include "configs.h"

using namespace casadi;
namespace E = Eigen;

class InformationCost : public SQPModule{
private:
    int id;
    std::array<int, 2> sl_gradient_x;
    std::array<int, 2> sl_gradient_y;
    std::array<int, 2> sl_value;
    SQPParams* params;
    InfoField* field;
public:
    int addToVariables(int start_idx) noexcept override;

public:
    int addToParams(int start_idx) noexcept;
    void initializeSQP(MX &cost_fun, MX &ca_states, MX &ca_u, MX &delta_z,
                       MX &delta_u,
                       OptiWrapper *opti_obj, DM &states, MX &ca_param_vec,
                       DM &param_vec) override;
    void iterationUpdate(OptiWrapper *opti_obj) override{}
    void linearize(OptiWrapper *opti_obj, DM &states, DM &param_vec) override;
    void warmStart(OptiWrapper *opti_obj, DM &states, DM &param_vec) override{}
    InformationCost(SQPParams* params, int id, InfoField* field);
    bool checkFeasible(Opti *opti_obj, DM &states){
        return true;
    }
    void draw(sf::RenderTarget& target, sf::RenderStates states) const{}
};
#endif //SQP_INFORMATION_COST_H
