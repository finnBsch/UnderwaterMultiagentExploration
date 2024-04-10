//
// Created by finn on 5/31/23.
//

#ifndef SQP_PATH_H
#define SQP_PATH_H
#include <vector>
#include <cmath>
#include <array>
#include <SFML/Graphics.hpp>
#include "spline.h"





class SQPObstacle: public sf::Drawable{
protected:
    sf::VertexArray line;
public:
    float x0;
    float y0;
    float x1;
    float y1;
    float length;
    float angle;
    bool left_block = true;
    SQPObstacle(float x0, float y0, float x1, float y1);
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    int dir_t = 0;
    // DEBUG
    bool isactive = true;

};

struct PathSegment : public sf::Drawable {
    double x0;
    double y0;
    double x1;
    double y1;
    double fl0 = -.5;
    double fr0 = 0.5;
    double dalpha_l = 0.0;
    double dalpha_r = 0.0;
    double angle;
    double theta0;
    double theta1;
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    sf::VertexArray line;
    sf::CircleShape c0;
    sf::CircleShape c1;
    PathSegment(double x0, double y0, double x1, double y1, double angle,
                double theta0, double theta1, double fl0, double fr0,
                double dalpha_l, double dalpha_r);
    PathSegment(const PathSegment& t);
};

class Path : public sf::Drawable {
protected:
    int num_segments;
    double max_theta;
    std::vector<PathSegment> segments;
public:
    Path();
    std::vector<PathSegment>& getSegments();
    void buildPath(std::vector<std::array<double, 2>>& points, bool reversed);
    void buildCirclePath();
    // TODO: Make faster
    std::array<double, 8> linSegment(double theta);
    double getTheta(double x, double y);
    double getTheta(double x, double y, double minTheta, double maxTheta);
    double getMaxTheta() const;
    void draw(sf::RenderTarget &target, sf::RenderStates states) const;
};

class SplinePath : public sf::Drawable{
protected:
    double size_x;
    double size_y;
    std::vector<SQPObstacle*>* obstacles;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> theta;
    tk::spline sx;
    tk::spline sy;

public:
    void buildTestPath();
    void postUpdate();
    SplinePath(std::vector<SQPObstacle*>* obstacles);
    std::array<double, 8> linSegment(double theta) const;
    double getMaxTheta() const;
    void draw(sf::RenderTarget &target, sf::RenderStates states) const;


};
#endif //SQP_PATH_H
