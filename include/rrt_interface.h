//
// Created by finn on 6/22/23.
//

#ifndef MASTERMAIN_RRT_INTERFACE_H
#define MASTERMAIN_RRT_INTERFACE_H
#include "rrt.h"
#include "path.h"



class RRTSplinePathInterface : public SplinePath {
private:
    RRT* rrt;
    int current_root_id = 0;
public:
    explicit RRTSplinePathInterface(RRT* rrt, std::vector<SQPObstacle*>* obs);
    void update(int offset_replan, int min_id);
    const std::vector<double>* getThetas() const;
    const std::vector<double>* getXs() const;
    const std::vector<double>* getYs() const;
    void getPos(int id, double& x, double& y, double& angle) const;
    double getAngle(int id_to_root);
    int getRootId() const{
        return current_root_id;
    }

protected:
//    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};

class RRTPathInterface : public Path {
private:
    int num_segs = 0;
    RRT* rrt;
public:
    explicit RRTPathInterface(RRT* rrt);
    void draw(sf::RenderTarget &target, sf::RenderStates states) const;
    void update(int start_id, int min_id);
};

#endif //MASTERMAIN_RRT_INTERFACE_H
