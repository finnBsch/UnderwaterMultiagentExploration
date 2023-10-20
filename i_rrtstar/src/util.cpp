//
// Created by finn on 5/12/23.
//
#include <random>
#include "util.h"

float sampleScalar(float &max) {
    static thread_local std::mt19937 generator(std::random_device{}());
//    static thread_local std::mt19937 generator(123);
    std::uniform_real_distribution<float> distribution(0.0f,max);
    return distribution(generator);
}

float euclideanDistance(const float &x0, const float &y0, const float &x1,
                        const float &y1) {
    return sqrtf(powf(x1 - x0, 2) + powf(y1 - y0, 2));
}

float squaredDistance(const float &x0, const float &y0, const float &x1,
                      const float &y1) {
    return powf(x1 - x0, 2) + powf(y1 - y0, 2);
}

void
linearInterp(const float &x0, const float &y0, const float &x1, const float &y1,
             const float &dist, float& x_interp, float & y_interp) {
    float fac = dist/euclideanDistance(x0, y0, x1, y1);
    x_interp = (x1 - x0) * fac + x0;
    y_interp = (y1 - y0) * fac + y0;
}


bool doIntersect(float px1, float py1, float qx1, float qy1,
                 float px2, float py2, float qx2, float qy2){
    // Find the four orientations needed for general and
    // special cases
    int o1 = orientation(px1, py1, qx1, qy1, px2, py2);
    int o2 = orientation(px1, py1, qx1, qy1, qx2, qy2);
    int o3 = orientation(px2, py2, qx2, qy2, px1, py1);
    int o4 = orientation(px2, py2, qx2, qy2, qx1, qy1);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    // Special Cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(px1, py1, px2, py2, qx1, qy1)) return true;

    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(px1, py1, qx2, qy2, qx1, qy1)) return true;

    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(px2, py2, px1, py1, qx2, qy2)) return true;

    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(px2, py2, qx1, qy1, qx2, qy2)) return true;

    return false; // Doesn't fall in any of the above cases
}


bool onSegment(float px, float py, float qx, float qy, float rx, float ry){
    if (qx <= std::max(px, rx) && qx >= std::min(px, rx) &&
        qy <= std::max(py, ry) && qy >= std::min(py, ry)) {
        return true;
    }
    return false;
}


int orientation(float px, float py, float qx, float qy, float rx, float ry) {
    float val = (qy - py) * (rx - qx) -
                (qx - px) * (ry - qy);

    if (val == 0) return 0;  // collinear

    return (val > 0)? 1: 2; // clock or counterclock wise
}