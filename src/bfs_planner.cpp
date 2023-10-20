//
// Created by finn on 9/2/23.
//
#include <queue>
#include <iostream>
#include "bfs_planner.h"


BFSPlanner::BFSPlanner(std::vector<RRTObstacle *> *obstacles, int N_X, int N_Y, double x0, double y0, double sx, double sy) {
    this->obstacles = obstacles;
    this->N_X = N_X;
    this->N_Y = N_Y;
    this->x = std::floor(x0/sx *(N_X-1));
    this->y = std::floor(y0/sy *(N_Y-1));
    this->delta_x = sx/(N_X - 1);
    this->delta_y = sy/(N_Y - 1);

}
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);

        if (hash1 != hash2) {
            return hash1 ^ hash2;
        }

        // If hash1 == hash2, their XOR is zero.
        return hash1;
    }
};
/** \brief Get the next action to take using breadth first search
 *
 */
void BFSPlanner::step(int goal_x, int goal_y) {
    if(x == goal_x && y == goal_y){
//        // make random move, sample number between 0 and 3
//        int action = rand() % 4;
//        if(action == 0){
//            x += 1;
//        }
//        else if(action == 1){
//            x -= 1;
//        }
//        else if(action == 2){
//            y += 1;
//        }
//        else if(action == 3){
//            y -= 1;
//        }
        return;
    }
    std::set<std::pair<int, int>> visited;
    std::queue<std::pair<int, int>> q;
    std::map<std::pair<int, int>, std::pair<int, int>> parent;
    std::pair<int, int> start = std::make_pair(x, y);
    q.push(start);
    visited.insert(start);
//    std::cout << "[BFS] Finding to " << goal_x << ", " << goal_y << std::endl;
    while(!q.empty()){
        std::pair<int, int> curr = q.front();
        q.pop();
        if(curr.first == goal_x && curr.second == goal_y){
//            std::cout << "x" << std::endl;
            break;
        }
        // check the four neighboring nodes, NOT all 8
        std::pair<int, int> neighbors[4] = {std::make_pair(curr.first + 1, curr.second),
                                            std::make_pair(curr.first - 1, curr.second),
                                            std::make_pair(curr.first, curr.second + 1),
                                            std::make_pair(curr.first, curr.second - 1)};
        for(auto next : neighbors){
            if(!visited.contains(next) && next.first >= 0 && next.first < N_X && next.second >= 0 && next.second < N_Y){
                if(collisionFree(curr, next)){
                    q.push(next);
                    visited.insert(next);
                    parent[next] = curr;
                }
            }
        }
    }
    std::pair<int, int> curr = std::make_pair(goal_x, goal_y);
    while(parent[curr] != start){
        curr = parent[curr];
    }
    x = curr.first;
    y = curr.second;
}

bool BFSPlanner::collisionFree(std::pair<int, int> p1, std::pair<int, int> p2) {
    if(p2.first < 0 || p2.first >= N_X || p2.second < 0 || p2.second >= N_Y){
        return false;
    }
    double x_eps1 = 0;
    double y_eps1 = 0;
    double x_eps2 = 0;
    double y_eps2 = 0;
    if(p2.first < p1.first){
        x_eps2 = 0.01;
        x_eps1 = -0.01;
    }
    else if(p2.first > p1.first){
        x_eps2 = -0.01;
        x_eps1 = 0.01;
    }
    if(p2.second < p1.second){
        y_eps2 = 0.01;
        y_eps1 = -0.01;
    }
    else if(p2.second > p1.second){
        y_eps2 = -0.01;
        y_eps1 = 0.01;
    }
    if(p2.second == 0){
        y_eps2 = 0.01;
    }
    if(p2.first == 0){
        x_eps2 = 0.01;
    }
    if(p2.second == N_Y - 1){
        y_eps2 = -0.01;
    }
    if(p2.first == N_X - 1){
        x_eps2 = -0.01;
    }
    // same for p1
    if(p1.first == 0){
        x_eps1 = 0.01;
    }
    if(p1.second == 0){
        y_eps1 = 0.01;
    }
    if(p1.first == N_X - 1){
        x_eps1 = -0.01;
    }
    if(p1.second == N_Y - 1){
        y_eps1 = -0.01;
    }

    for(auto obstacle : *obstacles){
        if(doIntersect(p1.first * delta_x + x_eps1, p1.second * delta_y + y_eps1, p2.first * delta_x + x_eps2, p2.second * delta_y + y_eps2, obstacle->x0, obstacle->y0, obstacle->x1, obstacle->y1)){
            return false;
        }
    }
    return true;
}

double BFSPlanner::getX() {
    return x*delta_x;
}

double BFSPlanner::getY() {
    return y*delta_y;
}

