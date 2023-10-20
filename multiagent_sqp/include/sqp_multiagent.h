//
// Created by finn on 5/22/23.
//

#ifndef CASADI_TEST_SQP_BUILDER_H
#define CASADI_TEST_SQP_BUILDER_H
#include <casadi/casadi.hpp>
#include <eigen3/Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "bicycle_dynamics.h"
#include "connectivity_constraint.h"
#include "configs.h"
#include "collision_constraint.h"
#include "information_cost.h"
#include "path_tracking.h"
#include "com_path_tracking.h"
#include "corridor_constraint.h"
#include "nearest_obstacle_avoidance.h"
#include <SFML/Graphics.hpp>
#include "distribution_cost.h"
#include "opti_wrapper.h"


using namespace casadi;
namespace E = Eigen;


class SQPMultiAgent : public sf::Drawable{
private:

    InfoField* info_field;
    SplinePath* path;
    SQPParams params;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_states;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigen_u;

    std::array<int, 2> sl_states;
    std::array<int, 2> sl_u;

    DM states_vals;
    DM u_vals;

    OptiWrapper opti;

    std::vector<SQPModule*> sqp_modules;
    std::vector<ConnectivityConstraint*> connectivity_constraints;
    std::vector<CollisionConstraint*> collision_constraints;
    std::vector<CorridorConstraint*> corridor_constraints;
    std::vector<Dynamics*> agents;
    COMPathTracking* path_tracking;
    std::vector<sf::VertexArray> plans;

    // Debug
    Function f_Hess;

public:
    double getTheta() const;
    double getMinTheta() const;
    void solve(const Eigen::Matrix<double, Eigen::Dynamic, 1>& init_states, bool reset);
    SQPMultiAgent(SQPParams params_, SplinePath *path,
                  InfoField *info_field,
                  Eigen::Matrix<double, Eigen::Dynamic, 1>& initial_state,
                  std::vector<SQPObstacle*>* obstacles);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& getSol();
    void draw(sf::RenderTarget& target, sf::RenderStates states) const;

};
#endif //CASADI_TEST_SQP_BUILDER_H
