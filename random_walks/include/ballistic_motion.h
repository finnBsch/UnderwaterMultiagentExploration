//
// Created by finn on 7/27/23.
//

#ifndef MASTERMAIN_BALLISTIC_MOTION_H
#define MASTERMAIN_BALLISTIC_MOTION_H
#include "random_walk.h"

class BallisticMotion : public RandomWalk {
private:
    std::uniform_real_distribution<double> distribution;
    std::vector<int> num_rotations;
    double yaw_rate = M_PI/4;
public:
    BallisticMotion(Eigen::Matrix<double, Eigen::Dynamic, 1>& states, double dt,
                    std::vector<RandomWalkObstacle*>& obstacles, int num_agents);
    void step();
};

#endif //MASTERMAIN_BALLISTIC_MOTION_H
