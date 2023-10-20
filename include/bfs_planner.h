//
// Created by finn on 9/2/23.
//

#ifndef MASTERMAIN_BFS_PLANNER_H
#define MASTERMAIN_BFS_PLANNER_H
#include "obstacle_interface.h"
#include "random_walk.h"
class BFSPlanner{
private:
    std::vector<RRTObstacle *> * obstacles;
    int N_X;
    int N_Y;
    int x;
    int y;
    double delta_x;
    double delta_y;
public:
    BFSPlanner(std::vector<RRTObstacle *> * obstacles, int N_X, int N_Y, double x0, double y0, double sx, double sy);
    void step(int goal_x, int goal_y);
    bool collisionFree(std::pair<int, int> p1, std::pair<int, int> p2);
    double getX();
    double getY();

};

#endif //MASTERMAIN_BFS_PLANNER_H
