//
// Created by finn on 7/25/23.
//
#include <iostream>

#include "random_walk.h"

bool isLeft2(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax)*(cy - ay) - (by - ay)*(cx - ax) >= 0;
}

bool onSegment2(float px, float py, float qx, float qy, float rx, float ry){
    if (qx <= std::max(px, rx) && qx >= std::min(px, rx) &&
        qy <= std::max(py, ry) && qy >= std::min(py, ry)) {
        return true;
    }
    return false;
}

int orientation2(float px, float py, float qx, float qy, float rx, float ry) {
    float val = (qy - py) * (rx - qx) -
                (qx - px) * (ry - qy);

    if (val == 0) return 0;  // collinear

    return (val > 0)? 1: 2; // clock or counterclock wise
}

bool doIntersect2(float px1, float py1, float qx1, float qy1,
                  float px2, float py2, float qx2, float qy2){
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation2(px1, py1, qx1, qy1, px2, py2);
    int o2 = orientation2(px1, py1, qx1, qy1, qx2, qy2);
    int o3 = orientation2(px2, py2, qx2, qy2, px1, py1);
    int o4 = orientation2(px2, py2, qx2, qy2, qx1, qy1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment2(px1, py1, px2, py2, qx1, qy1)) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment2(px1, py1, qx2, qy2, qx1, qy1)) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment2(px2, py2, px1, py1, qx2, qy2)) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment2(px2, py2, qx1, qy1, qx2, qy2)) return true;

    return false; // Doesn't fall in any of the above cases
}

RandomWalk::RandomWalk(Eigen::Matrix<double, Eigen::Dynamic, 1>& states,
                       double dt,
                       std::vector<RandomWalkObstacle *> &obstacles, int num_agents):
        generator(std::random_device{}()){
    this->states = states;
    this->dt = dt;
    this->obstacles = obstacles;
    this->num_agents = num_agents;
    needed_correction = Eigen::Matrix<bool, Eigen::Dynamic, 1>::Zero(num_agents);
}

const Eigen::Matrix<double, Eigen::Dynamic, 1> &RandomWalk::getStates() {
    return states;
}

/**
 * Checks if a step is collision free.
 * @param agent_id
 * @param x0
 * @param y0
 * @param x1
 * @param y1
 * @return
 */
bool RandomWalk::collisionFree(int agent_id, double x0, double y0, double x1, double y1, double& xn, double& yn) {
    // First pairwise with other agents
    for(int i = 0; i < num_agents; i++){
        if(i!=agent_id){
            double x2 = states(3 * i);
            double y2 = states(3 * i + 1);
            double dx = x1 - x0;
            double dy = y1 - y0;
            double d = sqrt(dx * dx + dy * dy);
            if(d < collision_radius){
                xn = dx;
                yn = dy;
                return false;
            }
        }
    }
    // Then with line obstacles. Needs 2 checks: Circle/Linesegment for x1, y1
    //              and line/line for x0->x1, y0->y1 with obstacle
    for(int i = 0; i < obstacles.size(); i++){
        double x2 = obstacles[i]->x0;
        double y2 = obstacles[i]->y0;
        double x3 = obstacles[i]->x1;
        double y3 = obstacles[i]->y1;
        // get obstacle dir
        double dx_o = x3 - x2;
        double dy_o = y3 - y2;
        // Check if circle intersects line segments
        double dx = x2 - x1;
        double dy = y2 - y1;
        double d_sq = dx_o * dx_o + dy_o * dy_o; // Lenght squared of the obstacle
        double t = std::max(0.0, std::min(1.0, (dx_o * dx + dy_o * dy)/d_sq));
        double px = x2 + t * dx_o;
        double py = y2 + t * dy_o;
        if (sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1)) < collision_radius){
            xn = x3 - x2;
            yn = y3 - y2;
            return false;
        }
        // Check if line intersects line segments
        if (doIntersect2(x0, y0, x1, y1, x2, y2, x3, y3)){
            xn = x3 - x2;
            yn = y3 - y2;
            return false;
        }
    }
    return true;
}

std::vector<std::array<double, 2>>
RandomWalk::getBlockedDirections(int agent_id, double x0, double y0, double x1,
                                 double y1) {
    std::vector<std::array<double, 2>> blocked_directions;
    // First pairwise with other agents
    for(int i = 0; i < num_agents; i++){
        if(i!=agent_id){
            double x2 = states(3 * i);
            double y2 = states(3 * i + 1);
            double dx = x2 - x1;
            double dy = y2 - y1;
            double d = sqrt(dx * dx + dy * dy);
            if(d < collision_radius){
                std::array<double, 2> dir = {dx, dy};
                blocked_directions.push_back(dir);
            }
        }
    }
    // Then with line obstacles. Needs 2 checks: Circle/Linesegment for x1, y1
    //              and line/line for x0->x1, y0->y1 with obstacle
    for(int i = 0; i < obstacles.size(); i++){
        double x2 = obstacles[i]->x0;
        double y2 = obstacles[i]->y0;
        double x3 = obstacles[i]->x1;
        double y3 = obstacles[i]->y1;
        // get obstacle dir
        double dx_o = x3 - x2;
        double dy_o = y3 - y2;
        // Check if circle intersects line segments
        double dx = x2 - x1;
        double dy = y2 - y1;
        double d_sq = dx_o * dx_o + dy_o * dy_o; // Lenght squared of the obstacle
        double t = std::max(0.0, std::min(1.0, (dx_o * dx + dy_o * dy)/d_sq));
        double px = x2 + t * dx_o;
        double py = y2 + t * dy_o;
        if (sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1)) < collision_radius){
            double xn = x3 - x2;
            double yn = y3 - y2;
            std::array<double, 2> dir = {xn, yn};
            blocked_directions.push_back(dir);
        }
        // Check if line intersects line segments
        if (doIntersect2(x0, y0, x1, y1, x2, y2, x3, y3)){
            double xn = x3 - x2;
            double yn = y3 - y2;
            std::array<double, 2> dir = {xn, yn};
            blocked_directions.push_back(dir);
        }
    }
    return blocked_directions;
}

bool RandomWalk::correctAgentPairs() {
    bool collision_free = false;
    bool required_correction = false;
    // Get all possible agent pairs
    std::vector<std::array<int, 2>> agent_pairs;
    for(int i = 0; i < num_agents; i++){
        for(int j = i+1; j < num_agents; j++){
            std::array<int, 2> pair = {i, j};
            agent_pairs.push_back(pair);
        }
    }
    while(!collision_free){
        collision_free = true;
        for (auto& pair: agent_pairs){
            int i = pair[0];
            int j = pair[1];
            double x1 = states(3 * j);
            double y1 = states(3 * j + 1);
            double x2 = states(3 * i);
            double y2 = states(3 * i + 1);
            double dx = x2 - x1;
            double dy = y2 - y1;
            double d = sqrt(dx * dx + dy * dy);
            if(d < collision_radius*2) {
                needed_correction(i) = true;
                needed_correction(j) = true;
                double difference = collision_radius*2 - d + 0.01;
                // Correct the position such that the distance is collision_radius
                x1 -= 0.5*difference * dx / d;
                y1 -= 0.5*difference * dy / d;
                x2 += 0.5*difference * dx / d;
                y2 += 0.5*difference * dy / d;
                collision_free = false;
                required_correction = true;
            }
            states(3 * j) = x1;
            states(3 * j + 1) = y1;
            states(3 * i) = x2;
            states(3 * i + 1) = y2;
        }
    }
    return required_correction;
}

bool RandomWalk::correctWallPositions(Eigen::MatrixXd& old_states) {
    bool required_correction = false;
    for (int a = 0; a < num_agents; a++) {
        double x1 = states(3 * a);
        double y1 = states(3 * a + 1);
        double x0 = old_states(3 * a);
        double y0 = old_states(3 * a + 1);
        bool collision_free = false;
        while (!collision_free) {
            collision_free = true;
            // Now check with line obstacles
            for (int i = 0; i < obstacles.size(); i++) {
                double x2 = obstacles[i]->x0;
                double y2 = obstacles[i]->y0;
                double x3 = obstacles[i]->x1;
                double y3 = obstacles[i]->y1;
                // get obstacle dir
                double dx_o = x3 - x2;
                double dy_o = y3 - y2;
                // Check if circle intersects line segments
                double dx = x1 - x2;
                double dy = y1 - y2;
                double d_sq =
                        dx_o * dx_o + dy_o * dy_o; // Length squared of the
                double dx_o_n = dx_o / sqrt(d_sq); // Normalized obstacle dir
                double dy_o_n = dy_o / sqrt(d_sq);
                // obstacle
                double t = std::max(0.0,
                                    std::min(1.0,
                                             (dx_o * dx + dy_o * dy) / d_sq));
                double px = x2 + t * dx_o;
                double py = y2 + t * dy_o;
                double d = sqrt((px - x1) * (px - x1) + (py - y1) * (py - y1));
                if (d < collision_radius) {
                    needed_correction(a) = true;
                    collision_free = false;
                    required_correction = true;
                    // Then move the circle to the correct position. Need to
                    // determine if it was left or right before in order to move
                    // it to the correct side.
                    bool was_left = isLeft2(x2, y2, x3, y3, x0, y0);
                    bool is_left = isLeft2(x2, y2, x3, y3, x1, y1);
                    double difference;
                    if (was_left == is_left) {
                        difference = collision_radius - d + 0.01;
                    } else {
                        difference = d + collision_radius + 0.01;
                    }
                    if (was_left) {
                        x1 += difference * -dy_o_n;
                        y1 += difference * dx_o_n;
                    } else {
                        x1 += difference * dy_o_n;
                        y1 += difference * -dx_o_n;
                    }
                }
            }
        }
        states(3 * a) = x1;
        states(3 * a + 1) = y1;
    }
    return required_correction;
}
