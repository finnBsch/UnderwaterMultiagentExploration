//
// Created by finn on 6/6/23.
//

#ifndef I_RRTSTAR_VIZ_H
#define I_RRTSTAR_VIZ_H
#include <SFML/Graphics.hpp>
#include "rrt.h"


class Viz {
private:
    sf::ContextSettings settings;
    RRTParams* params;
    RRT* rrt;
    sf::View view;
    sf::RenderWindow window;
public:
    Viz(RRTParams* params, RRT* rrt);
    bool draw();
};

#endif //I_RRTSTAR_VIZ_H
