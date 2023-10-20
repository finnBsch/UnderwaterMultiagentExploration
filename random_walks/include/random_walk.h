//
// Created by finn on 7/25/23.
//

#ifndef MASTERMAIN_RANDOM_WALK_H
#define MASTERMAIN_RANDOM_WALK_H

#include <vector>
#include <eigen3/Eigen/Dense>
#include <random>

struct RandomWalkObstacle {
    double x0;
    double y0;
    double x1;
    double y1;

};

class RandomWalk {
protected:
    Eigen::Matrix<bool, Eigen::Dynamic, 1> needed_correction;
    std::mt19937 generator;
    Eigen::Matrix<double, Eigen::Dynamic, 1> states;
    double dt;
    int num_agents = 0;
    std::vector<RandomWalkObstacle*> obstacles;
    double collision_radius = 0.2;
    bool collisionFree(int agent_id, double x0, double y0, double x1, double y1, double& xn, double& yn);
    std::vector<std::array<double, 2>> getBlockedDirections(int agent_id, double x0, double y0, double x1, double y1);
    bool correctAgentPairs();
    bool correctWallPositions(Eigen::MatrixXd &old_states);
public:
    virtual void step() = 0;
    RandomWalk(Eigen::Matrix<double, Eigen::Dynamic, 1>& states, double dt,
               std::vector<RandomWalkObstacle*>& obstacles, int num_agents);
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& getStates();
};

#endif //MASTERMAIN_RANDOM_WALK_H
