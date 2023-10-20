//
// Created by finn on 6/6/23.
//
#include "viz.h"
#include <iostream>

Viz::Viz(RRTParams *params, RRT *rrt):window(), params(params), rrt(rrt),
                                   view(sf::Vector2f((float)params->size_x/2, (float)params->size_y/2),
                                        sf::Vector2f((float)params->size_x*1.2f, - (float)params->size_y*1.2f)),
                                   settings(0, 0, 8){
    window.create(sf::VideoMode (1000, 500), "RRT", sf::Style::None,
                  settings);
    window.setView(view);
    window.setPosition(sf::Vector2i(0, 100));
}

bool Viz::draw() {
    window.clear(sf::Color::Black);
    sf::Event ev;
    bool a = false;
    while (window.pollEvent(ev)) {
        if (ev.type == sf::Event::MouseButtonPressed) {
            a = true;
        }
    }
    window.draw(*rrt);
    window.display();
    return a;
}

