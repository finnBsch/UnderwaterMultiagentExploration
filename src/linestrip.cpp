//
// Created by finn on 10/11/23.
//
#include <cassert>
#include <cmath>
#include "linestrip.h"

/**
 * Construct a linestrip given a set of points and a thickness
 * @param x_coords
 * @param y_coords
 * @param thickness
 */
LineStrip::LineStrip(std::vector<float> x_coords, std::vector<float> y_coords,
                     float thickness) {
    this->x_coords = x_coords;
    this->y_coords = y_coords;
    this->thickness = thickness;
    triangle_strip.setPrimitiveType(sf::TriangleStrip);
}

/**
 * Construct a linestrip given a thickness
 * @param thickness
 */
LineStrip::LineStrip(float thickness) {
    this->thickness = thickness;
    triangle_strip.setPrimitiveType(sf::TriangleStrip);
}

/**
 * Overwrite all points in the linestrip
 * @param x_coords
 * @param y_coords
 */
void
LineStrip::setPoints(std::vector<float> x_coords, std::vector<float> y_coords) {
    this->x_coords = x_coords;
    this->y_coords = y_coords;
    triangle_strip.resize(2 * x_coords.size());
    updateTriangleStrip();
}

/**
 * Appends to the points of the linestrip. Hereby, x_coords is expected to be
 * longer than the previous x_coords such that just the new points are appended
 * @param x_coords
 * @param y_coords
 */
void LineStrip::updatePoints(std::vector<float> x_coords,
                             std::vector<float> y_coords) {
    assert(x_coords.size() >= this->x_coords.size());
    assert(y_coords.size() >= this->y_coords.size());
    int id0 = this->x_coords.size() - 1;
    this->x_coords = x_coords;
    this->y_coords = y_coords;
    triangle_strip.resize(2 * x_coords.size());
    updateTriangleStrip(id0);
}

/**
 * Updates the triangle strip from the given id0
 * @param id0
 */
void LineStrip::updateTriangleStrip(int id0) {
    int num_pts = x_coords.size();
    for (int i = 0; i < num_pts; i++){
        if (i >= id0) {
            float angle;
            if (i == 0) {
                angle = atan2f(y_coords[i + 1] - y_coords[i],
                               x_coords[i + 1] - x_coords[i]);
            } else if (i == x_coords.size() - 1) {
                angle = atan2f(y_coords[i] - y_coords[i - 1],
                               x_coords[i] - x_coords[i - 1]);
            } else {
                // TODO Make sure that this is correct
                angle = atan2f(y_coords[i + 1] - y_coords[i - 1],
                               x_coords[i + 1] - x_coords[i - 1]);
            }
            float x1 = x_coords[i] + thickness * cosf(angle + M_PIf / 2.0f);
            float y1 = y_coords[i] + thickness * sinf(angle + M_PIf / 2.0f);
            float x2 = x_coords[i] + thickness * cosf(angle - M_PIf / 2.0f);
            float y2 = y_coords[i] + thickness * sinf(angle - M_PIf / 2.0f);
            triangle_strip[2 * i].position = sf::Vector2f(x1, y1);
            triangle_strip[2 * i + 1].position = sf::Vector2f(x2, y2);
        }
        triangle_strip[2 * i].color = sf::Color(255, 255, 255, (int)round((double)i/(double)num_pts * 255));
        triangle_strip[2 * i + 1].color = sf::Color(255, 255, 255, (int)round((double)i/(double)num_pts * 255));

    }
}

void LineStrip::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    target.draw(triangle_strip, states);
}

void DashedLineStrip::updateTriangles(int id0) {
    // lengths can adopt the length of x_coords
    int num_pts = x_coords.size();
    lengths.resize(num_pts);
    for (int i = 0; i < num_pts; i++){
        if(i>=id0){
            if (i == 0){
                lengths[i] = 0.0;
            }
            else {
                lengths[i] = sqrtf(powf(x_coords[i] - x_coords[i - 1], 2.0f) +
                                  powf(y_coords[i] - y_coords[i - 1], 2.0f)) + lengths[i - 1];
            }
        }
    }
    int num_dash_segments = (int) ceilf(2.0f * lengths[num_pts - 1] / (dash_length + dash_gap));
    triangles.resize(6 * num_dash_segments);
    for (int i = 0; i < num_dash_segments; i++){
        float length0 = (float)i * (dash_length + dash_gap);
        float length1 = (float)i * (dash_length + dash_gap) + dash_length;
        float x0;
        float y0;
        float x1;
        float y1;
        get_coords(length0, x0, y0);
        get_coords(length1, x1, y1);
        float angle = atan2f(y1 - y0, x1 - x0);
        float x0_0 = x0 + thickness * cosf(angle + M_PIf / 2.0f);
        float y0_0 = y0 + thickness * sinf(angle + M_PIf / 2.0f);
        float x0_1 = x0 + thickness * cosf(angle - M_PIf / 2.0f);
        float y0_1 = y0 + thickness * sinf(angle - M_PIf / 2.0f);
        float x1_0 = x1 + thickness * cosf(angle + M_PIf / 2.0f);
        float y1_0 = y1 + thickness * sinf(angle + M_PIf / 2.0f);
        float x1_1 = x1 + thickness * cosf(angle - M_PIf / 2.0f);
        float y1_1 = y1 + thickness * sinf(angle - M_PIf / 2.0f);
        triangles[6 * i].position = sf::Vector2f(x0_0, y0_0);
        triangles[6 * i + 1].position = sf::Vector2f(x0_1, y0_1);
        triangles[6 * i + 2].position = sf::Vector2f(x1_0, y1_0);
        triangles[6 * i + 3].position = sf::Vector2f(x0_1, y0_1);
        triangles[6 * i + 4].position = sf::Vector2f(x1_0, y1_0);
        triangles[6 * i + 5].position = sf::Vector2f(x1_1, y1_1);
        triangles[6 * i].color = sf::Color(255, 255, 255, 255);
        triangles[6 * i + 1].color = sf::Color(255, 255, 255, 255);
        triangles[6 * i + 2].color = sf::Color(255, 255, 255, 255);
        triangles[6 * i + 3].color = sf::Color(255, 255, 255, 255);
        triangles[6 * i + 4].color = sf::Color(255, 255, 255, 255);
        triangles[6 * i + 5].color = sf::Color(255, 255, 255, 255);
    }
}

DashedLineStrip::DashedLineStrip(std::vector<float> x_coords,
                                 std::vector<float> y_coords, float thickness,
                                 float dash_length):
                                 triangles(sf::Triangles){
    this->thickness = thickness;
    this->dash_length = dash_length;
    setPoints(x_coords, y_coords);
}

DashedLineStrip::DashedLineStrip(float thickness, float dash_length):
                                triangles(sf::Triangles){
    this->thickness = thickness;
    this->dash_length = dash_length;
    this->dash_gap = dash_length;
}

DashedLineStrip::DashedLineStrip(float thickness, float dash_length,
                                 float dash_gap):triangles(sf::Triangles) {
    this->thickness = thickness;
    this->dash_length = dash_length;
    this->dash_gap = dash_gap;
}

void DashedLineStrip::setPoints(std::vector<float> x_coords,
                                std::vector<float> y_coords) {
    this->x_coords = x_coords;
    this->y_coords = y_coords;
    updateTriangles();
}

void DashedLineStrip::updatePoints(std::vector<float> x_coords,
                                   std::vector<float> y_coords) {
    this->x_coords = x_coords;
    this->y_coords = y_coords;
    int id0 = this->x_coords.size() - 1;
    updateTriangles(id0);
}

void
DashedLineStrip::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    target.draw(triangles, states);
}

void DashedLineStrip::get_coords(float length, float &x, float &y) {
    if (length <= 0){
        x = x_coords[0];
        y = y_coords[0];
        return;
    }
    if (length >= lengths[lengths.size() - 1]){
        x = x_coords[x_coords.size() - 1];
        y = y_coords[y_coords.size() - 1];
        return;
    }
    int id0 = 0;
    while (lengths[id0] < length){
        id0++;
    }
    // Found the first entry with length > length, now need to interpolate
    float length0 = lengths[id0 - 1];
    float length1 = lengths[id0];
    float x0 = x_coords[id0 - 1];
    float x1 = x_coords[id0];
    float y0 = y_coords[id0 - 1];
    float y1 = y_coords[id0];
    float frac = (length - length0) / (length1 - length0);
    x = x0 + frac * (x1 - x0);
    y = y0 + frac * (y1 - y0);
}
