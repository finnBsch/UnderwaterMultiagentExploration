//
// Created by finn on 6/22/23.
//

#include <iostream>
#include "rrt_interface.h"
#include "spline.h"


RRTPathInterface::RRTPathInterface(RRT *rrt) {
    this->rrt = rrt;
    update(0.0, 0);
}

void RRTPathInterface::update(int start_id, int min_id) {
    auto segs = rrt->getSegments();
    int id0 = start_id;
    min_id = std::max(min_id, 0);
//    segments.erase(segments.begin(), segments.begin() + min_id);
//    start_id = start_id - min_id;
    int id_max = segments.size() - 1;
    double theta0;

    if(start_id == 0){
        theta0 = 0.0;
    }
    else {
        theta0 = segments[id0].theta0;
    }
    int counter = 0;
//    std::cout <<"ROOT " <<  segs->back().x1 << std::endl;
    if(!segments.empty()) {
        if (id0 <= id_max) {
//            std::cout << segments[id0].x1 << std::endl;
        } else {
//            std::cout << segments.back().x1 << std::endl;
        }
    }
    for(int i = segs->size() - 1; i >= 0; i--){
        double x1 = segs->operator[](i).x0;
        double y1 = segs->operator[](i).y0;
        double x0 = segs->operator[](i).x1;
        double y0 = segs->operator[](i).y1;

        double angle = atan2(y1 - y0, x1 - x0);
        double length_seg = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
        if(length_seg < 0.05){
            continue;
        }
        double theta1 = theta0 + length_seg;
        double fl0 = segs->operator[](i).lower_bound;
        double fr0 = segs->operator[](i).upper_bound;
        double dalpha_l = 0;
        double dalpha_r = 0;
        if(counter + id0 <= id_max){
            segments[counter + id0] = PathSegment(x0, y0, x1, y1, angle, theta0,
                                                  theta1, fl0, fr0, dalpha_l,
                                                  dalpha_r);
        }
        else {
            segments.emplace_back(x0, y0, x1, y1, angle, theta0, theta1, fl0,
                                  fr0, dalpha_l, dalpha_r);
        }
        counter += 1;
        theta0 = theta1;
    }
    num_segs = counter + id0 - 1;
    segments.erase(segments.begin() + num_segs, segments.end());

//    for(int i = 0; i < segments.size() - 1; i++){
//        std::cout << segments[i].x1 << "\t";
//        if(segments[i + 1].x0 != segments[i].x1
//        || segments[i + 1].y0 != segments[i].y1
//        || segments[i + 1].theta0 != segments[i].theta1
//        || segments[i + 1].theta1 < segments[i + 1].theta0){
//            std::cout << "ERROR" << std::endl;
//        }
//    }
//    std::cout << std::endl;

}

void RRTPathInterface::draw(sf::RenderTarget &target,
                            sf::RenderStates states) const {
    Path::draw(target, states);
}

//void RRTSplinePathInterface::draw(sf::RenderTarget &target,
//                                  sf::RenderStates states) const {
//    for(auto& obs : *obstacles){
//        target.draw(*obs, states);
//    }
//
//}

RRTSplinePathInterface::RRTSplinePathInterface(RRT *rrt,
                                               std::vector<SQPObstacle*>*
                                                       obs): SplinePath(obs) {
    this->rrt = rrt;
    size_x = rrt->getSizeX();
    size_y = rrt->getSizeY();

    update(0.0, 0);
}


void RRTSplinePathInterface::update(int offset_replan, int min_id) {
    auto segs = rrt->getSegments();
    int id0 = offset_replan + current_root_id;

    int id_max = theta.size() - 1;
    double num_segs;
    double theta0;
    double x0;
    double y0;
    if(id0 == 0){
        theta0 = 0.0;
        x0 = segs->back().x1;
        y0 = segs->back().y1;
        if(theta.size() == 0){
            theta.emplace_back(theta0);
            x.emplace_back(x0);
            y.emplace_back(y0);
        }
        else{
            theta[0] = theta0;
            x[0] = x0;
            y[0] = y0;
        }
    }
    else {
        theta0 = theta[id0];
        x0 = x[id0];
        y0 = y[id0];
        if(x0 != segs->back().x1){
            std::cout << "ERROR" << std::endl;
        }
    }
    current_root_id = id0;
    int counter = 1;
    for(int i = segs->size() - 1; i >= 0; i--){
        double x1 = segs->operator[](i).x0;
        double y1 = segs->operator[](i).y0;
        double seg_length = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
        if(seg_length > 0.3 or true) {
            double theta1 = theta0 + seg_length;
            if (counter + id0 <= id_max) {
                theta[counter + id0] = theta1;
                x[counter + id0] = x1;
                y[counter + id0] = y1;
            } else {
                theta.emplace_back(theta1);
                x.emplace_back(x1);
                y.emplace_back(y1);
            }
            counter += 1;
            theta0 = theta1;
            x0 = x1;
            y0 = y1;
        }
    }
    num_segs = counter + id0;
    theta.erase(theta.begin() + num_segs, theta.end());
    x.erase(x.begin() + num_segs, x.end());
    y.erase(y.begin() + num_segs, y.end());

    sx.set_points(theta, x);
    sy.set_points(theta, y);
    for(int i = 1; i < x.size(); i++){
//        std::cout << theta[i]  - theta[i - 1] << ", ";
        if(x[i] == x[i - 1] && y[i] == y[i - 1]){
            std::cout << "ERROR 1" << std::endl;
        }
    }
//    std::cout << std::endl;
    postUpdate();
}

const std::vector<double> *RRTSplinePathInterface::getThetas() const {
    return &theta;
}

const std::vector<double> *RRTSplinePathInterface::getXs() const {
    return &x;
}

const std::vector<double> *RRTSplinePathInterface::getYs() const {
    return &y;
}

void RRTSplinePathInterface::getPos(int id, double &x, double &y,
                                    double &angle) const {
    x = this->x[id];
    y = this->y[id];
    angle = atan2(this->y[id + 1] - this->y[id], this->x[id + 1] - this->x[id]);

}

double RRTSplinePathInterface::getAngle(int id_to_root) {
    return atan2(sy.deriv(1, theta[id_to_root+current_root_id]), sx.deriv(1, id_to_root+current_root_id));
}
