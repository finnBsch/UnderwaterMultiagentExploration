//
// Created by finn on 7/25/23.
//

#ifndef MASTERMAIN_BROWNIAN_RANDOM_WALK_H
#define MASTERMAIN_BROWNIAN_RANDOM_WALK_H
#include "random_walk.h"


class BrownianRandomWalk : public RandomWalk {
private:
//    int num_steps = 2;
    std::vector<int> num_steps;
    std::uniform_real_distribution<double> distribution;
public:
    BrownianRandomWalk(Eigen::Matrix<double, Eigen::Dynamic, 1>& states, double dt,
                       std::vector<RandomWalkObstacle*>& obstacles, int num_agents);
    void step();
};


#endif //MASTERMAIN_BROWNIAN_RANDOM_WALK_H
