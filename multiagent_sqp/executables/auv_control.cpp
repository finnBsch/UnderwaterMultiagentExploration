//
// Created by finn on 5/18/23.
//
// TODO: Test code-gen, extend to multi, test SQP..

#include "bicycle_dynamics.h"
#include <casadi/casadi.hpp>
#include <SFML/Graphics.hpp>
#include <eigen3/Eigen/Dense>
#include "../include/sqp_multiagent.h"
#include <chrono>
using namespace std::chrono;

using namespace casadi;
namespace E = Eigen;

int main(){
    float sx = 20;
    float sy = 20;
    Slice all;  // Basic slicer
    BicycleDynamics agent;
    int N = 10;
    double dt = 0.1;
    Opti opti = Opti();

    MX x = opti.variable(4, N + 1);
    MX x2 = opti.variable(4, N + 1);
    MX x_init = opti.parameter(4);
    MX x2_init = opti.parameter(4);
    MX u = opti.variable(2, N);
    MX u2 = opti.variable(2, N);
    opti.minimize(sum2(pow(x(0, all) - (4.9 + 1.2), 2) + pow(x(1, all) - 12.9, 2))
                  + sum2(pow(x2(0, all) - 4.9, 2) + pow(x2(1, all) - 12.9, 2)));  // Minimize

    for(int k = 0; k < N; k++){
        MX x_next = x(all, k) + dt * agent.getDeriv(x(all, k), u(all, k));
        MX x_next2 = x2(all, k) + dt * agent.getDeriv(x2(all, k), u2(all, k));
        opti.subject_to(x(all, k + 1) == x_next);
        opti.subject_to(x2(all, k + 1) == x_next2);
    }
    opti.subject_to(-agent.a_params.speed_max <= x(3,  all)
                    <= agent.a_params.speed_max);
    opti.subject_to(-agent.a_params.yaw_rate_max <= u(1,all) <= agent.a_params
            .yaw_rate_max);
    opti.subject_to(-agent.a_params.accel_max <= u(0, all) <= agent.a_params
            .accel_max);
    opti.subject_to(1 <= pow(x(0, all) - x2(0, all), 2) + pow(x(1, all) - x2(1, all), 2) <= 4);
    opti.subject_to(x(0, 0) == x_init(0));
    opti.subject_to(x(1, 0) == x_init(1));
    opti.subject_to(x(2, 0) == x_init(2));
    opti.subject_to(x(3, 0) == x_init(3));

    opti.subject_to(-agent.a_params.speed_max <= x2(3,  all)
                    <= agent.a_params.speed_max);
    opti.subject_to(-agent.a_params.yaw_rate_max <= u2(1,all) <= agent.a_params
            .yaw_rate_max);
    opti.subject_to(-agent.a_params.accel_max <= u2(0, all) <= agent.a_params
            .accel_max);
    opti.subject_to(x2(0, 0) == x2_init(0));
    opti.subject_to(x2(1, 0) == x2_init(1));
    opti.subject_to(x2(2, 0) == x2_init(2));
    opti.subject_to(x2(3, 0) == x2_init(3));
//    opti.subject_to(x(0, N) == 4.9);
//    opti.subject_to(x(1, N) == 12.9);
//    opti.subject_to(x(3, N) == 0.0);

    opti.set_initial(x, 0.0);
    opti.set_initial(x2, 0.0);
    opti.solver("ipopt");//, Dict({{"print_time", false}}), Dict({{"print_level", 0}}));
    e_state state;
    e_state state2;
    state.setZero();
    state2.setZero();
    state2(0) = 1.1;
    state2(2) = 0.1;
    DM st = DM::zeros(4, 1);
    DM st2 = DM::zeros(4, 1);
    std::memcpy(st.ptr(), state.data(), sizeof(double)*4*1);
    std::memcpy(st2.ptr(), state2.data(), sizeof(double)*4*1);
    opti.set_value(x_init, st);
    opti.set_value(x2_init, st2);
    OptiSol sol = opti.solve();
    std::cout << sol.value(x) << "\n";
    auto sol_x1 = sol.value(x);
    auto sol_x2 = sol.value(x2);
    auto sol_u = sol.value(u);
    auto sol_u2 = sol.value(u2);
    sf::ContextSettings settings( 0, 0, 8);
    auto window = sf::RenderWindow(sf::VideoMode (1000, 1000), "blah",
                                   sf::Style::None, settings);
    auto view = sf::View(sf::Vector2f(sx/2, sy/2), sf::Vector2f(sx, -sy));
    window.setView(view);
    auto start = high_resolution_clock::now();
    auto circ = sf::CircleShape(0.2);
    auto circ2 = sf::CircleShape(0.2);
    float t = 0;
    window.setPosition(sf::Vector2i(0, 0));
    e_u u_rt;
    e_u u_rt2;
    float t_last = 0.0f;
    t = 0.0f;
    float t_last_control = 0.0f;
    int counter = 0;
    while(true){
        window.clear(sf::Color::Black);
        sf::Event ev;
        bool a = false;
        while (window.pollEvent(ev)) {
            if (ev.type == sf::Event::MouseButtonPressed) {
                a = true;
            }
        }
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
//        t = (float)duration.count()/2000.0f;
        if(t - t_last > 0.01) {
            int id = floor((t - t_last_control) / dt);
            u_rt(0, 0) = sol_u(0, 0)->data()[0];
            u_rt(1, 0) = sol_u(1, 0)->data()[0];

            u_rt2(0, 0) = sol_u2(0, 0)->data()[0];
            u_rt2(1, 0) = sol_u2(1, 0)->data()[0];
            state += agent.getDeriv(state, u_rt) *dt;
            std::cout << state << " " << sol_x1 << std::endl;
            state2 += agent.getDeriv(state2, u_rt2) *dt;
            std::cout << state2 << " " << sol_x2 << std::endl;
            t_last = t;
            circ.setPosition(sf::Vector2f(state(0, 0), state(1, 0)));
            circ2.setPosition(sf::Vector2f(state2(0, 0), state2(1, 0)));
            std::cout << "x: " << state(0,0) << ", y: " << state(1, 0) << " " << state(3, 0) << "\n";
            window.draw(circ);
            window.draw(circ2);
            window.display();
        }
        std::cout << counter << std::endl;
        t_last_control = t;
        std::memcpy(st.ptr(), state.data(), sizeof(double)*4*1);
        std::memcpy(st2.ptr(), state2.data(), sizeof(double)*4*1);
        opti.set_value(x_init, st);
        opti.set_value(x2_init, st2);
        sol = opti.solve();
        sol_x1 = sol.value(x);
        sol_x2 = sol.value(x2);
        sol_u = sol.value(u);
        sol_u2 = sol.value(u2);

        t += dt;
        counter += 1;
    }

}