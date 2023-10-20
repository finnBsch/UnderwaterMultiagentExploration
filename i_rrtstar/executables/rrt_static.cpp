#include <iostream>
#include "rrt.h"
#include "viz.h"
#include <cstdlib>

#include <chrono>
using namespace std::chrono;
// simple test script
int main() {
    RRTParams params;
    std::vector<RRTObstacle*> obstacles;
    // Small rectangle in the center
    float xc = 5.0;
    float yc = 2.0;
    float w = 1.0;
    float h = 2.5;
    obstacles.push_back(new RRTObstacle(xc-w/2, yc-h/2, xc-w/2, yc + h/2));
    obstacles.push_back(new RRTObstacle(xc-w/2, yc + h/2, xc + w/2, yc + h/2));
    obstacles.push_back(new RRTObstacle(xc + w/2, yc + h/2, xc + w/2, yc - h/2));
    obstacles.push_back(new RRTObstacle(xc + w/2, yc - h/2, xc - w/2, yc - h/2));


//    obstacles.push_back(new SQPObstacle(6.0f, 5.0f, 11.0f, 5.0f));
    RRT rrt = RRT(&params, &obstacles, new SineTestField, 2.5, 2.0);
    auto start = high_resolution_clock::now();
    auto start2 = high_resolution_clock::now();
    rrt.buildTree();
    rrt.updateTree();
    rrt.updateFromRoot();
    rrt.updateFromRoot();
    rrt.updateFromRoot();
    rrt.updateFromRoot();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << std::endl;
    start = high_resolution_clock::now();
    Viz v(&params, &rrt);
    float pos = 0;
    while(!rrt.rewireDone()){}
    rrt.exportTree("/home/finn/CLionProjects/mastersthesis_main/Experiments"
                   "/tree.csv");
    v.draw();
    auto goal_node = rrt.getPathTo(25.0f, 20.0f);
    while(true){
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        if(duration.count() > 10){
            rrt.buildPathWithCorridor(goal_node);
            v.draw();
        }
    }
}
