//
// Created by finn on 5/26/23.
//
#include "bicycle_dynamics.h"
#include <casadi/casadi.hpp>
#include <SFML/Graphics.hpp>
#include <eigen3/Eigen/Dense>
#include "../include/sqp_multiagent.h"
#include <chrono>
#include "info_field.h"

using namespace std::chrono;

using namespace casadi;
namespace E = Eigen;

int main() {
    TestField info_field;
    SplinePath path;
    int N = 200;
    SQPParams params;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> trajectory;
    trajectory.resize(params.num_states * params.num_agents, N);
    Eigen::Matrix<double, Eigen::Dynamic, 1> current_state;
    current_state.resize(params.num_states * params.num_agents, 1);
    current_state.setZero();
    for(int i = 0; i < params.num_agents; i++){
        current_state(i * params.num_states, 0) = i * 0.45;
    }
    float sx = 15;
    float sy = 15;
    std::vector<SQPObstacle*> obstacles;
    SQPMultiAgent opt(params, &path, &info_field,
                      current_state, &obstacles);
    opt.solve(current_state, true);
    sf::ContextSettings settings( 0, 0, 8);
    auto window = sf::RenderWindow(sf::VideoMode (1000, 1000), "sqp_test",
                                   sf::Style::None, settings);
    auto view = sf::View(sf::Vector2f(0, sy/2), sf::Vector2f(sx, -sy));
    window.setView(view);
    std::vector<sf::CircleShape> circs;
    for(int i = 0; i < params.num_agents; i++){
        circs.emplace_back(0.2);
        circs[i].setOrigin(0.2, 0.2);
    }
    window.setPosition(sf::Vector2i(0, 0));
    auto sol = opt.getSol();
    auto start = high_resolution_clock::now();
    float t;
    int counter = 0;
    trajectory.col(counter) = current_state;
    while(true){
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        t = (float)duration.count()/1000.0f;
        if(t > params.dt){
            start = stop;
            current_state(Eigen::seq(0, params.num_states * params.num_agents - 1)) = sol(Eigen::seq(0, params.num_states * params.num_agents - 1), 1);
            counter += 1;
            if(counter > N-1){
                std::cout << "[";
                for(int i = 0; i < trajectory.rows(); i++){
                    std::cout << "[";
                    std::cout << trajectory(i, 0);
                    for(int j = 1; j < trajectory.cols(); j++){
                        std::cout << ", " << trajectory(i, j);
                    }
                    if(i < trajectory.rows() - 1) {
                        std::cout << "],\n";
                    }

                }
                std::cout << "]]\n";
                return 0;
            }
            else {
                trajectory.col(counter) = current_state;
            }
            opt.solve(current_state, false);
            sol = opt.getSol();
        }
        window.clear(sf::Color::Black);
        sf::Event ev;
        bool a = false;
        while (window.pollEvent(ev)) {
            if (ev.type == sf::Event::MouseButtonPressed) {
                a = true;
            }
        }
        int id = floor(t/0.1);
        if(id>sol.cols() - 2){
            return 0;
        }
        float dt0 = id * params.dt;
        float dt1 = (id + 1) * params.dt;
        float dt = t - dt0;
        for(int i = 0; i < params.num_agents; i++){
            circs[i].setPosition((sol(i * params.num_states, id + 1) - sol(i*params.num_states, id))*dt/
            (dt1-dt0) + sol(i*params.num_states, id), (sol(i * params.num_states +1, id + 1) - sol(i*params.num_states +1,
                                                                   id))*dt/
                                      (dt1-dt0) + sol(i*params.num_states+ 1, id) );
            window.draw(circs[i]);
            window.draw(opt);
        }
        window.display();
    }
}