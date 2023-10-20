//
// Created by finn on 6/25/23.
//
#include "rrt_util.h"


bool isLeft(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax)*(cy - ay) - (by - ay)*(cx - ax) >= 0;
}
