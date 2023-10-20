//
// Created by finn on 6/23/23.
//
#include <iostream>
#include "obstacle_interface.h"

ObstacleInterface::ObstacleInterface(std::vector<RRTObstacle *>*
        obstacles_rrt) {
    for(auto & obs : *obstacles_rrt){
        obstacles.push_back(new SQPObstacle(obs->x0, obs->y0, obs->x1, obs->y1));
        obstacles_rw.push_back(new RandomWalkObstacle{obs->x0, obs->y0,
                                                      obs->x1, obs->y1});
    }
}

std::vector<SQPObstacle *> *ObstacleInterface::getObstacles() {
    return &obstacles;
}

std::vector<RandomWalkObstacle *> *ObstacleInterface::getObstaclesRW() {
    return &obstacles_rw;
}

ObstacleInterface::~ObstacleInterface() {
    for(auto & obs : obstacles){
//        std::cout << "x" << std::endl;
//        std::cout << obs->x0 << std::endl;
        delete obs;
    }
    for(auto & obs : obstacles_rw){
        delete obs;
    }
}
