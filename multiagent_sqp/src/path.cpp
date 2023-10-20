//
// Created by finn on 5/31/23.
//
#include "path.h"
#include <iostream>
#include <SFML/Graphics/RenderTarget.hpp>

#include "rrt_util.h"
inline double getAbsoluteDiff2Angles(const double x, const double y)
{
    // c can be PI (for radians) or 180.0 (for degrees);
    return M_PI - fabs(fmod(fabs(x - y), 2*M_PI) - M_PI);
}

inline double constrainAngle(double x){
    x = fmod(x + 180,360);
    if (x < 0)
        x += 360;
    return x - 180;
}


SQPObstacle::SQPObstacle(float x0, float y0, float x1, float y1): x0(x0), y0(y0),
                                                                  x1(x1), y1(y1), line(sf::Lines, 2) {
    if(this->x0 < this->x1){
        float temp = this->x0;
        this->x0 = this->x1;
        this->x1 = temp;
        temp = this->y0;
        this->y0 = this->y1;
        this->y1 = temp;
    }
    else if(this->x0 == this->x1 && this->y0 > this->y1){
        float temp = this->y0;
        this->y0 = this->y1;
        this->y1 = temp;
    }
    // TODO MAke angle work correctly
    length = sqrtf(pow(x1 - x0, 2) + pow(y1 - y0, 2));
    angle = atan2f(y1 - y0, x1 - x0);
    line[0].position = sf::Vector2f(x0, y0);
    line[0].color = sf::Color(120, 120, 120);
    line[1].position = sf::Vector2f(x1, y1);
    line[1].color = sf::Color(120, 120, 120);
}

void SQPObstacle::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    sf::VertexArray line_(sf::Lines, 2);

    line_[0].position = sf::Vector2f(x0 - 0.1, y0 - 0.1);
    line_[1].position = sf::Vector2f(x1 - 0.1, y1 - 0.1);
    sf::CircleShape blockside(0.1);
    double normal_x = -(y1 - y0) / length;
    double normal_y = (x1 - x0) / length;
    blockside.setOrigin(0.1, 0.1);
    if(dir_t == 1) {
        line_[0].color = sf::Color(10, 10, 0);
        line_[1].color = sf::Color(255, 255, 0);
    }
    else if(dir_t == -1){
        line_[0].color = sf::Color(255, 255, 0);
        line_[1].color = sf::Color(10, 10, 0);
    }
    else{
        line_[0].color = sf::Color(120, 120, 120);
        line_[1].color = sf::Color(120, 120, 120);
    }
    if (dir_t != 0) {
        if (left_block) {
            blockside.setPosition((x1 + x0) / 2 + normal_x * 0.4,
                                  (y1 + y0) / 2 +
                                  normal_y * 0.4);
        } else {
            blockside.setPosition((x1 + x0) / 2 - normal_x * 0.4,
                                  (y1 + y0) / 2 -
                                  normal_y * 0.4);
        }
        target.draw(blockside, states);
    }
    target.draw(line_, states);
}

Path::Path() {
//    buildCirclePath();
}

std::array<double, 8> Path::linSegment(double theta) {
    if(theta < 0){
        return std::array<double, 8>{segments[0].x0, segments[0].y0,
                                     segments[0].theta0, segments[0].angle,
                                     segments[0].fl0, segments[0].fr0, segments[0].dalpha_l,
                                     segments[0].dalpha_r};
    }
    for(auto& seg : segments){
        if(seg.theta0 <= theta && seg.theta1 > theta){
//            std::cout << "x0 " << seg.x0 << " y0 " << seg.y0 << " theta0 " << seg.theta0 << " angle " << seg.angle << " fl0 " << seg.fl0 << " fr0 " << seg.fr0 << " dalpha_l " << seg.dalpha_l << " dalpha_r " << seg.dalpha_r << std::endl;
            return std::array<double, 8>{seg.x0, seg.y0, seg.theta0, seg
                    .angle, seg.fl0, seg.fr0, seg.dalpha_l, seg.dalpha_r};
        }
    }
    auto seg = segments.back();
//    std::cout << "WW " << std::endl;
    return std::array<double, 8>{seg.x0, seg.y0, seg.theta0, seg
            .angle, seg.fl0, seg.fr0, seg.dalpha_l, seg.dalpha_r};
}

//double Path::getTheta(double x, double y) {
//    double min_dist = pow(segments[0].x0 - x, 2) + pow(segments[0].y0 - y, 2);
//    double d_check;
//    int min_id = -1;
//    for(int i = 0; i < segments.size(); i++){
//        auto & seg = segments[i];
//        double d =pow(seg.x1 - x, 2) + pow(seg.y1 - y, 2);
//        if(d < min_dist){
//            min_dist = d;
//            min_id = i;
//        }
//    }
//    if(min_id == -1){
//        min_id = 0;
//        d_check = min_dist;
//    }
//    else if(min_id == segments.size() - 1){
//        d_check = pow(segments[min_id].x0 - x, 2) + pow
//                (segments[min_id].y0 -y, 2);
//    }
//    else{
//        // check next and prev.
//        double d0 = pow(segments[min_id].x0 - x, 2) + pow(segments[min_id].y0 -
//                                                                 y, 2);
//        double d1 = pow(segments[min_id + 1].x1 - x, 2) + pow(segments[min_id
//                                                                       + 1].y1 - y, 2);
//        if(d1 < d0){
//            min_id = min_id + 1;
//            d_check = min_dist;
//        }
//        else{
//            d_check = d0;
//        }
//    }
//    double d_angle = atan2(y - segments[min_id].y0, x - segments[min_id].x0)
//            - segments[min_id].angle;
//    return fmin(segments[min_id].theta0 + sqrt(d_check) * cos(d_angle),
//               segments.back().theta1);
//}

double Path::getTheta(double x, double y) {
    double min_dist = pow(segments[0].x0 - x, 2) + pow(segments[0].y0 - y, 2);
    double theta0 = segments[0].theta0;
    for(auto & seg : segments){
        double d = pow(seg.x0 - seg.x1, 2) + pow(seg.y0 - seg.y1, 2);
        double t = ((x - seg.x0)*(seg.x1 - seg.x0) + (y - seg.y0)*(seg.y1 -
                                                                   seg.y0))/d;
        if(t >= 0 && t <= 1){
            double xp = seg.x0 + t * (seg.x1 - seg.x0);
            double yp = seg.y0 + t * (seg.y1 - seg.y0);
            double dist = pow(x - xp, 2) + pow(y - yp, 2);
            if(dist < min_dist){
                min_dist = dist;
                theta0 = seg.theta0 + t * (seg.theta1 - seg.theta0);
            }
        }
    }
    return theta0;
}

void Path::buildCirclePath() {
    segments.clear();
    num_segments = 10;
    double e_c_max = 0.5;
    double r = 6.7;
    double cx = 0;
    double cy = r;
    double max_angle = 4.0/3.0 * M_PI;// * 3.0/4.0;
    double d_alpha = max_angle/((double)num_segments - 1.0);
    double x0 = 0;
    double y0 = 0;
    double theta0 = 0;
    for(int i = 0; i < num_segments; i++){
        double m = i%2;
        double x1 = cx + sin(d_alpha * (i + 1))*r;
        double y1 = cy - cos(d_alpha * (i + 1))*r;
        double length_seg = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
        double ang_corridor = atan2(m * 2.0, length_seg) * 0;
        double ang = atan2(y1 - y0, x1 - x0);
        double theta1 =  theta0 + length_seg;
        segments.emplace_back(x0, y0, x1, y1, ang, theta0,
                              theta1, -(e_c_max + m*2.0), + e_c_max + m*2.0, -ang_corridor, ang_corridor);
//        std::cout << "Point at x: " << x1 << ", y: " << y1 << ", angle: " <<
//                  ang * 180 / M_PI << ", theta0: " << theta0 << ", theta1: "
//                                                                "" << theta1<<
//                  "\n";
        theta0 = theta1;
        x0 = x1;
        y0 = y1;

    }
}

std::vector<PathSegment> &Path::getSegments() {
    return segments;
}

void
Path::buildPath(std::vector<std::array<double, 2>> &points, bool reversed) {

    segments.clear();
    if(reversed){
        double x0 = points.back()[0];
        double y0 = points.back()[1];
        double theta0 = 0;
        for(int i = points.size() - 2; i >= 0; i--){
            double x1 = points[i][0];
            double y1 = points[i][1];
            double angle = atan2(y1 - y0, x1 - x0);
            double length_seg = sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
            double theta1 = theta0 + length_seg;
            double fl0 = -0.5;
            double fr0 = 0.5;
            double dalpha_l = 0;
            double dalpha_r = 0;
            segments.emplace_back(x0, y0, x1, y1, angle, theta0, theta1, fl0,
                                  fr0, dalpha_l, dalpha_r);
            theta0 = theta1;
            x0 = x1;
            y0 = y1;
        }
    }
    else{

    }
}

void Path::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for(auto & seg : segments){
        target.draw(seg, states);
    }
}

double Path::getMaxTheta() const {
    return segments.back().theta1;
}

double Path::getTheta(double x, double y, double minTheta, double maxTheta) {
    double min_dist = pow(segments[0].x0 - x, 2) + pow(segments[0].y0 - y, 2);
    double theta0 = minTheta;
    for(auto & seg : segments){
        if(seg.theta0 < minTheta){
            continue;
        }
        double d = pow(seg.x0 - seg.x1, 2) + pow(seg.y0 - seg.y1, 2);
        double t = ((x - seg.x0)*(seg.x1 - seg.x0) + (y - seg.y0)*(seg.y1 -
                                                                   seg.y0))/d;
        if(t >= 0 && t <= 1){
            double xp = seg.x0 + t * (seg.x1 - seg.x0);
            double yp = seg.y0 + t * (seg.y1 - seg.y0);
            double dist = pow(x - xp, 2) + pow(y - yp, 2);
            if(dist < min_dist){
                min_dist = dist;
                theta0 = seg.theta0 + t * (seg.theta1 - seg.theta0);
            }
        }
    }
    return theta0;
}

PathSegment::PathSegment(double x0, double y0, double x1, double y1, double angle,
                         double theta0, double theta1, double fl0, double fr0,
                         double dalpha_l, double dalpha_r): line(sf::Lines,
                                                                 2),
                                                            c0(0.05),
                                                            c1(0.05){
    this->x0 = x0;
    this->y0 = y0;
    this->x1 = x1;
    this->y1 = y1;
    this->angle = angle;
    this->theta0 = theta0;
    this->theta1 = theta1;
    this->fl0 = fl0;
    this->fr0 = fr0;
    this->dalpha_l = dalpha_l;
    this->dalpha_r = dalpha_r;
    line[0].position = sf::Vector2f(x0, y0);
    line[1].position = sf::Vector2f(x1, y1);
    c0.setOrigin(0.05, 0.05);
    c1.setOrigin(0.05, 0.05);
    c0.setPosition(x0, y0);
    c1.setPosition(x1, y1);
}

void
PathSegment::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    target.draw(line, states);
    target.draw(c0, states);
    target.draw(c1, states);
}

PathSegment::PathSegment(const PathSegment &t) {
    this->x0 = t.x0;
    this->y0 = t.y0;
    this->x1 = t.x1;
    this->y1 = t.y1;
    this->angle = t.angle;
    this->theta0 = t.theta0;
    this->theta1 = t.theta1;
    this->fl0 = t.fl0;
    this->fr0 = t.fr0;
    this->dalpha_l = t.dalpha_l;
    this->dalpha_r = t.dalpha_r;
    this->line = t.line;
}

SplinePath::SplinePath(std::vector<SQPObstacle*>* obstacles) {
    this->obstacles = obstacles;
}

std::array<double, 8> SplinePath::linSegment(double theta) const {
//    if(theta < this->theta[0]){
//        theta = this->theta[0] + 0.1;
//    }
//    else if(theta > this->theta.back()){
    theta = std::max(0.01, theta);
    // TODO Fix these obstacle angles, add artificial segment length
    double x0 = sx(theta);
    double y0 = sy(theta);
    double theta0 = theta;
    double lower_bound = -std::numeric_limits<double>::max();
    double upper_bound = std::numeric_limits<double>::max();
    double delta_y = sy.deriv(1, theta);
    double delta_x = sx.deriv(1, theta);
    double angle_center = atan2(delta_y, delta_x);
    double d_angle_upper = 0.0;
    double d_angle_lower = 0.0;
    double virtual_length = 1.0;
    for(int d_a = -1; d_a < 2; d_a++) {
        double theta_ = theta + virtual_length/2 * d_a;
        double x0_ = sx(theta_);
        double y0_ = sy(theta_);
        double angle = atan2(sy.deriv(1, theta_), sx.deriv(1, theta_));
        double x_0 = x0_ - virtual_length/2 * cos(angle);
        double y_0 = y0_ - virtual_length/2 * sin(angle);


        for (auto &obs: *obstacles) {
            // Transform Obstacle into path segment coordinate system
            double x0_obs =
                    cos(angle) * (obs->x0 - x_0) + sin(angle) * (obs->y0 - y_0);
            double y0_obs = -sin(angle) * (obs->x0 - x_0) +
                            cos(angle) * (obs->y0 - y_0);
            double x1_obs =
                    cos(angle) * (obs->x1 - x_0) + sin(angle) * (obs->y1 - y_0);
            double y1_obs = -sin(angle) * (obs->x1 - x_0) +
                            cos(angle) * (obs->y1 - y_0);

            // Correct Sorting
            if (x0_obs > x1_obs) {
                std::swap(x0_obs, x1_obs);
                std::swap(y0_obs, y1_obs);
            }
            // Vertical Obstacle
            if (x0_obs == x1_obs) {
                if (x0_obs > 0.0 && x0_obs < virtual_length) {
                    double min_y = std::min(y0_obs, y1_obs);
                    if (min_y > 0) {
                        upper_bound = std::min(upper_bound, min_y);
                    }
                    double max_y = std::max(y0_obs, y1_obs);
                    if (max_y < 0) {
                        lower_bound = std::max(lower_bound, max_y);
                    }
                }
            }
                // Check if within bounds
            else if (x0_obs < virtual_length && x1_obs > 0) {
                double elevation = (y1_obs - y0_obs) / (x1_obs - x0_obs);
                double min_x = std::max(0.0, x0_obs);
                double max_x = std::min(virtual_length, x1_obs);
                double y_min_x = y0_obs + elevation * (min_x - x0_obs);
                double y_max_x = y0_obs + elevation * (max_x - x0_obs);
                double min_y = std::min(y_min_x, y_max_x);
                double max_y = std::max(y_min_x, y_max_x);

                if(min_y >= 0){
                    upper_bound = std::min(upper_bound, min_y);
                }
                if (max_y <= 0){
                    lower_bound = std::max(lower_bound, max_y);
                }
            }
        }
    }
    lower_bound += 0.5;
    upper_bound -= 0.5;
    std::swap(lower_bound, upper_bound);
    lower_bound *= -1;
    upper_bound *= -1;
//    std::cout << "x0 " << x0 << " y0 " << y0 << " theta0 " << theta0 << " angle " << angle_center << " lower_bound " << lower_bound << " upper_bound " << upper_bound << " sin d_angle_lower " << sin(d_angle_lower) << " sin d_angle_upper " << sin(d_angle_upper) << std::endl;
// seg.x0, seg.y0, seg.theta0, seg
//                    .angle, seg.fl0, seg.fr0, seg.dalpha_l, seg.dalpha_r
    return std::array<double, 8>({x0, y0, theta0, angle_center, lower_bound,
                                  upper_bound, d_angle_lower, d_angle_upper});
}



//std::array<double, 8> SplinePath::linSegment(double theta) const {
////    if(theta < this->theta[0]){
////        theta = this->theta[0] + 0.1;
////    }
////    else if(theta > this->theta.back()){
//    theta = std::max(0.01, theta);
//    // TODO Fix these obstacle angles, add artificial segment length
//    double x0 = sx(theta);
//    double y0 = sy(theta);
//    double theta0 = theta;
//    double lower_bound = -std::numeric_limits<double>::max();
//    double upper_bound = std::numeric_limits<double>::max();
//    double delta_y = sy.deriv(1, theta);
//    double delta_x = sx.deriv(1, theta);
//    double angle_center = atan2(delta_y, delta_x);
//    double d_angle_upper = 0.0;
//    double d_angle_lower = 0.0;
//    for(int d_a = -1; d_a < 2; d_a++) {
//        double angle = angle_center + d_a * 4.0 * M_PI/180.0;
//        double dx = cos(angle + M_PI / 2);
//        double dy = sin(angle + M_PI / 2);
//        double x1_normal = x0 + dx;
//        double y1_normal = y0 + dy;
//        double x1 = x0 + cos(angle);
//        double y1 = y0 + sin(angle);
//
//        for (auto &obs: *obstacles) {
//            double den = (obs->x0 * y0 - obs->x0 * y1_normal - obs->x1 * y0 +
//                          obs->x1 * y1_normal - obs->y0 * x0 +
//                          obs->y0 * x1_normal + obs->y1 * x0 -
//                          obs->y1 * x1_normal);
//            double t = (obs->x0 * y0 - obs->x0 * y1_normal - obs->y0 * x0 +
//                        obs->y0 * x1_normal + x0 * y1_normal - x1_normal * y0) /
//                       den;
//            if (t >= 0 && t <= 1) {
//                double angle_obs = obs->angle;
//                double xp = obs->x0 + t * (obs->x1 - obs->x0);
//                double yp = obs->y0 + t * (obs->y1 - obs->y0);
//                double dist = sqrt(pow(x0 - xp, 2) + pow(y0 - yp, 2));
//                bool left = isLeft(x0, y0, x1, y1, xp, yp);
//                if (!left && dist < upper_bound) {
//                    upper_bound = dist;
//                    double delta_y_obs = obs->y1 - obs->y0;
//                    double delta_x_obs = obs->x1 - obs->x0;
//                    d_angle_upper = 0;//atan2((delta_x * delta_y_obs - delta_y *delta_x_obs), (delta_x * delta_x_obs + delta_y * delta_y_obs));
//                } else if (left && -dist > lower_bound) {
//                    lower_bound = -dist;
//                    double delta_y_obs = obs->y1 - obs->y0;
//                    double delta_x_obs = obs->x1 - obs->x0;
//                    d_angle_lower = 0; //atan2((delta_x * delta_y_obs - delta_y *delta_x_obs), (delta_x * delta_x_obs + delta_y * delta_y_obs)) - M_PI;
//                }
//            }
//        }
//    }
//    lower_bound += 0.5;
//    upper_bound -= 0.5;
////    std::cout << "x0 " << x0 << " y0 " << y0 << " theta0 " << theta0 << " angle " << angle << " lower_bound " << lower_bound << " upper_bound " << upper_bound << " sin d_angle_lower " << sin(d_angle_lower) << " sin d_angle_upper " << sin(d_angle_upper) << std::endl;
//    // seg.x0, seg.y0, seg.theta0, seg
//    //                    .angle, seg.fl0, seg.fr0, seg.dalpha_l, seg.dalpha_r
//    return std::array<double, 8>({x0, y0, theta0, angle_center, lower_bound,
//                                  upper_bound, d_angle_lower, d_angle_upper});
//}



double SplinePath::getMaxTheta() const {
    return theta.back();
}

void SplinePath::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    sf::VertexArray line(sf::LinesStrip, theta.size()*10);
    double theta_max = theta.back();
    for(int i = 0; i < theta.size()*10; i++){
        double t = theta_max * i / (double)(theta.size()*10);
        line[i].position = sf::Vector2f(sx(t), sy(t));
    }
    target.draw(line);
//    for(int i = 0; i < theta.size(); i++){
//        sf::CircleShape circle(0.1 + 0.1 * (double)(i%2));
//        circle.setOrigin(0.1 + 0.1 * (double)(i%2), 0.1 + 0.1 * (double)(i%2));
//        circle.setPosition(x[i], y[i]);
//        target.draw(circle);
//    }
    for(auto & obs : *obstacles){
        target.draw(*obs, states);
    }
}

void SplinePath::postUpdate() {
    for(auto& obs: *obstacles){
        if(obs->x0 == obs->x1 && (obs->x0 == 0 || obs->x0 == size_x)){
            continue;
        }
        if(obs->y0 == obs->y1 && (obs->y0 == 0 || obs->y0 == size_y)){
            continue;
        }
        double min_dist_0 = std::numeric_limits<double>::max();
        double min_dist_1 = std::numeric_limits<double>::max();
        double min_id_0;
        double min_id_1;
        for(int i = 0; i < theta.size(); i++){
            double dist = pow(x[i] - obs->x0, 2) + pow(y[i] - obs->y0, 2);
            if(dist < min_dist_0){
                min_dist_0 = dist;
                min_id_0 = i;
            }
            dist = pow(x[i] - obs->x1, 2) + pow(y[i] - obs->y1, 2);
            if(dist < min_dist_1){
                min_dist_1 = dist;
                min_id_1 = i;
            }
        }
        double dx_path;
        double dy_path;
        bool is_left;
        if(min_dist_0 < min_dist_1){
            is_left = isLeft(x[min_id_0], y[min_id_0], x[min_id_0+1], y[min_id_0+1], obs->x0, obs->y0);
            dx_path = x[min_id_0+1] - x[min_id_0];
            dy_path = y[min_id_0+1] - y[min_id_0];
        }
        else{
            is_left = isLeft(x[min_id_1], y[min_id_1], x[min_id_1+1], y[min_id_1+1], obs->x1, obs->y1);
            dx_path = x[min_id_1+1] - x[min_id_1];
            dy_path = y[min_id_1+1] - y[min_id_1];
        }
        double dx_obs = obs->x1 - obs->x0;
        double dy_obs = obs->y1 - obs->y0;
        double cos_term = (dx_path * dx_obs + dy_path * dy_obs) / (sqrt(pow(dx_path, 2) + pow(dy_path, 2)) * sqrt(pow(dx_obs, 2) + pow(dy_obs, 2)));
        double angle = std::abs(acos(cos_term));
        double cross = dx_path * dy_obs - dy_path * dx_obs;
        int direction = 0;
        bool negative = false;
        if(cross < 0){
            angle = -angle;
            negative = true;
        }
        if(is_left && angle < 0) {
            direction = 1;
        }
        else if(is_left && angle > 0) {
            direction = -1;
        }
        else if(!is_left && angle < 0) {
            direction = -1;
        }
        else if(!is_left && angle > 0){
            direction = 1;
        }
        if(angle < 0 and is_left){
            angle = angle + M_PI;
        }
        else if(angle > 0 and !is_left){
            angle = angle - M_PI;
        }
        angle = std::abs(angle);
        bool problematic = angle < 2 * M_PI/4;
        if(problematic){
            obs->dir_t = -direction;
            if (negative) {
                obs->left_block = false;
            }
            else{
                obs->left_block = true;
            }
        }
        else{
            obs->dir_t = 0;
        }
    }

}
