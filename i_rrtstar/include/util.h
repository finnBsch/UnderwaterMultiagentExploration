//
// Created by finn on 6/6/23.
//

#ifndef I_RRTSTAR_UTIL_H
#define I_RRTSTAR_UTIL_H
#include <cmath>
float sampleScalar(float & max);
float euclideanDistance(const float& x0, const float& y0, const float& x1,
                        const float& y1);
float squaredDistance(const float& x0, const float& y0, const float& x1,
                      const float& y1);
void linearInterp(const float& x0, const float& y0, const float& x1,
                  const float& y1, const float& dist, float& x_interp,
                  float& y_interp);


bool doIntersect(float px1, float py1, float qx1, float qy1,
                 float px2, float py2, float qx2, float qy2);

bool onSegment(float px, float py, float qx, float qy, float rx, float ry);

int orientation(float px, float py, float qx, float qy, float rx, float ry);

inline float getAbsoluteDiff2Angles(const float x, const float y)
{
    // c can be PI (for radians) or 180.0 (for degrees);
    return M_PI - fabs(fmod(fabsf(x - y), 2*M_PI) - M_PI);
}

#endif //I_RRTSTAR_UTIL_H
