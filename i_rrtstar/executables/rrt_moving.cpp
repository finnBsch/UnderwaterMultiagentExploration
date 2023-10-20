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
    obstacles.push_back(new RRTObstacle(5.0f, -0.1f, 5.0f, 4.0f));
    obstacles.push_back(new RRTObstacle(6.0f, 5.0f, 11.0f, 5.0f));
    RRT rrt = RRT(&params, &obstacles, new SineTestField() ,0, 0);
    auto start = high_resolution_clock::now();
    auto start2 = high_resolution_clock::now();
    rrt.buildTree();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << std::endl;
    start = high_resolution_clock::now();
    Viz v(&params, &rrt);
    float pos = 0;
//    while(!rrt.rewireDone()){}
    v.draw();
    auto goal_node = rrt.getBestPath();
    while(true){
        stop = high_resolution_clock::now();
        duration = duration_cast<milliseconds>(stop - start);
        if(duration.count() > 10){
            start = stop;
            if(rrt.rewireDone()){
                pos += 0.1;
                if(pos >=10.0f){
                    std::cout << "bye";
                    stop = high_resolution_clock::now();
                    duration = duration_cast<milliseconds>(stop - start2);
                    std::cout << "Took : " << duration.count() << "ms" <<
                              std::endl;
                    return 0;
                }
                rrt.moveRoot(pos, pos, 0);
            }
            rrt.buildPathWithCorridor(goal_node);
            v.draw();
        }
    }
}
