//
// Created by finn on 7/27/23.
//
#include "ballistic_motion.h"

BallisticMotion::BallisticMotion(
        Eigen::Matrix<double, Eigen::Dynamic, 1> &states, double dt,
        std::vector<RandomWalkObstacle *> &obstacles, int num_agents)
        : RandomWalk(states, dt, obstacles, num_agents), distribution(0.0, 1.0) {
    for (int i = 0; i < num_agents; i++) {
        num_rotations.push_back(0);
    }
}

void BallisticMotion::step() {
    needed_correction = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Zero(num_agents);
    Eigen::MatrixXd old_states = states;
    for(int i = 0; i < num_agents; i++){
        double v = 1;
        if( num_rotations[i] > 0){
            // TODO correct random generator
            // Generate new angle at random and write into state
            states(3 * i + 2) += yaw_rate*dt;
            num_rotations[i] --;

        }
        else if (num_rotations[i] < 0){
            // TODO correct random generator
            // Generate new angle at random and write into state
            states(3 * i + 2) -= yaw_rate*dt;
            num_rotations[i] ++;
        }
        else {
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
        }
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
    for (int i = 0; i < num_agents; i++) {
        if (needed_correction(i)) {
            double num_steps_max = M_PI / (dt * yaw_rate);
            num_rotations[i] = (int) ((distribution(generator) * 2 - 1)
                                      * num_steps_max);
        }
    }
}
