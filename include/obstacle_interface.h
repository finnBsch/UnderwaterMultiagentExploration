//
// Created by finn on 6/23/23.
//

#ifndef MASTERMAIN_OBSTACLE_INTERFACE_H
#define MASTERMAIN_OBSTACLE_INTERFACE_H
#include "rrt.h"
#include "path.h"
#include "random_walk.h"

class ObstacleInterface {
private:
    std::vector<SQPObstacle*> obstacles;
    std::vector<RandomWalkObstacle*> obstacles_rw;

public:
    std::vector<SQPObstacle*>* getObstacles();
    std::vector<RandomWalkObstacle*>* getObstaclesRW();
    ObstacleInterface(std::vector<RRTObstacle *> * obstacles_rrt);
    ~ObstacleInterface();
};


#endif //MASTERMAIN_OBSTACLE_INTERFACE_H
