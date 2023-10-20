//
// Created by finn on 10/11/23.
//

#ifndef MASTERMAIN_LINESTRIP_H
#define MASTERMAIN_LINESTRIP_H
#include <SFML/Graphics.hpp>

class DashedLineStrip: public sf::Drawable, public sf::Transformable {
private:
    sf::VertexArray triangles; // The vertex array to draw the triangles. Cannot use triangle strip here for obvious reasons
    float thickness = 0.01;
    float dash_length = 0.1;
    float dash_gap = 0.1;
    std::vector<float> x_coords;
    std::vector<float> y_coords;
    std::vector<float> lengths;

    void updateTriangles(int id0 = 0);

    void get_coords(float length, float&x, float& y);
public:
    DashedLineStrip(std::vector<float> x_coords, std::vector<float> y_coords, float thickness, float dash_length);
    DashedLineStrip(float thickness, float dash_length);
    DashedLineStrip(float thickness, float dash_length, float dash_gap);

    void setPoints(std::vector<float> x_coords, std::vector<float> y_coords);
    void updatePoints(std::vector<float> x_coords, std::vector<float> y_coords);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};


class LineStrip: public sf::Drawable, public sf::Transformable {
private:
    sf::VertexArray triangle_strip; // The vertex array to draw as a triangle strip
    float thickness = 0.01;
    std::vector<float> x_coords;
    std::vector<float> y_coords;

    void updateTriangleStrip(int id0 = 0);
public:
    LineStrip(std::vector<float> x_coords, std::vector<float> y_coords, float thickness);
    LineStrip(float thickness);

    void setPoints(std::vector<float> x_coords, std::vector<float> y_coords);
    void updatePoints(std::vector<float> x_coords, std::vector<float> y_coords);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};

#endif //MASTERMAIN_LINESTRIP_H
