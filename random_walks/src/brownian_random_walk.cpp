//
// Created by finn on 7/25/23.
//
/**
 * Simple brownian random walk implementation.
 * TODO: At least some dynamics
 */
#include <iostream>
#include "../include/brownian_random_walk.h"

BrownianRandomWalk::BrownianRandomWalk(
        Eigen::Matrix<double, Eigen::Dynamic, 1> &states, double dt,
        std::vector<RandomWalkObstacle *> &obstacles , int num_agents) :
                                       RandomWalk(states, dt, obstacles, num_agents), distribution(0.0, 1.0){
    for (int i = 0; i < num_agents; i++) {
        num_steps.push_back(0);
    }

}
inline double getAbsoluteDiff2Angles(const double x, const double y)
{
    // c can be PI (for radians) or 180.0 (for degrees);
    return M_PI - fabs(fmod(fabs(x - y), 2*M_PI) - M_PI);
}

inline double constrainAngle(double x){
    x = fmod(x + 180,360);
    if (x < 0)
        x += 360;
    return x - 180;
}


void BrownianRandomWalk::step() {
    Eigen::MatrixXd old_states = states;
    for(int i = 0; i < num_agents; i++){
        double v = 1;
        if( num_steps[i] == 0){
            // TODO correct random generator
            // Generate new angle at random and write into state
            double d_angle = distribution(generator) *  M_PI/6 * 2 - M_PI/6;
            num_steps[i] = distribution(generator) * 10 + 1;
            states(3 * i + 2) += d_angle;

        }

        // Update the positions
        double x = states(3 * i);
        double y = states(3 * i + 1);
        double angle = states(3 * i + 2);
        double dx = cos(angle) * dt * v;
        double dy = sin(angle) * dt * v;
        double x1 = x + dx;
        double y1 = y + dy;
        states(3 * i) = x1;
        states(3 * i + 1) = y1;
        num_steps[i] -= 1;
    }
    bool collision_free = false;
    while(!collision_free) {
        collision_free = true;
        if(correctWallPositions(old_states)){
            collision_free = false;
        }
        if(correctAgentPairs()){
            collision_free = false;
        }
    }
}
