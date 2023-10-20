//
// Created by finn on 5/18/23.
//
// TODO: Test code-gen, extend to multi, test SQP..

#include "bicycle_dynamics.h"
#include <casadi/casadi.hpp>
#include <SFML/Graphics.hpp>

#include <chrono>
using namespace std::chrono;

using namespace casadi;

class NLPBuilder {
private:
    BicycleDynamics agent1;
    BicycleDynamics agent2;
    MX f;
    MX g;
    MX lbg;
    MX ubg;
    MX u1;
    MX x2;
    MX u2;
    MX x1_init;
    MX x2_init;
    template<typename T> void addLUBC(double lb, double ub, T expr);
    template<typename T> void addEQ(T expr);
    template<typename T> void addLBC(double lb, T expr);
    template<typename T> void addUBC(double ub, T expr);
public:
    MX x1;
    Function nlpsolver;
    NLPBuilder(int N, double dt);

};

NLPBuilder::NLPBuilder(int N, double dt) {
    Slice all;  // Basic slicer
    x1 = MX::sym("x1", 4, N + 1);
    u1 = MX::sym("u1", 2, N);
    x2 = MX::sym("x2", 4, N + 1);
    u2 = MX::sym("u2", 2, N);
    x1_init = MX::sym("x1_init", 4, 1);
    x2_init = MX::sym("x2_init", 4, 1);
    auto a = x1(0, all);
    f = sum2(pow(x1(0, all) - (4.9 + 1.2), 2) + pow(x1(1, all) - 12.9, 2))
        + sum2(pow(x2(0, all) - 4.9, 2) + pow(x2(1, all) - 12.9, 2));
    addEQ(x1(0, 0) - x1_init(0));
    addEQ(x1(1, 0) - x1_init(1));
    addEQ(x1(2, 0) - x1_init(2));
    addEQ(x1(3, 0) - x1_init(3));
    addEQ(x2(0, 0) - x2_init(0));
    addEQ(x2(1, 0) - x2_init(1));
    addEQ(x2(2, 0) - x2_init(2));
    addEQ(x2(3, 0) - x2_init(3));

    for(int k = 0; k < N; k++){
        MX x_next = x1(all, k) + dt * agent1.getDeriv(x1(all, k), u1(all, k));
        MX x_next2 = x2(all, k) + dt * agent2.getDeriv(x2(all, k), u2(all, k));
        addEQ(x_next - x1(all, k + 1));
        addEQ(x_next2 - x2(all, k + 1));
    }
    addLUBC(1, 4, pow(x1(0, all) - x2(0, all), 2)
                            + pow(x1(1, all) - x2(1, all), 2));
    addLUBC(-agent1.a_params.speed_max, agent1.a_params.speed_max, x1(3,  all));
    addLUBC(-agent2.a_params.speed_max, agent2.a_params.speed_max, x2(3,  all));
    addLUBC(-agent1.a_params.yaw_rate_max, agent1.a_params.yaw_rate_max, u1(1,all));
    addLUBC(-agent2.a_params.yaw_rate_max, agent2.a_params.yaw_rate_max, u2(1,all));
    addLUBC(-agent1.a_params.accel_max, agent1.a_params.accel_max,  u1(0, all));
    addLUBC(-agent2.a_params.accel_max, agent2.a_params.accel_max,  u2(0, all));

    MXDict nlp;
    nlp["x"] = vertcat( reshape(x1, (N+1)*4, 1), reshape(x2, (N+1)*4, 1),
                        reshape(u1, N*2, 1), reshape(u2, N*2, 1));
    nlp["p"] = vertcat(x1_init, x2_init);
    nlp["f"] = f;
    nlp["g"] = g;
    nlpsolver = nlpsol("my_nlp_sol", "ipopt", nlp);//, Dict({{"ipopt.print_level", 0}}));
}

template<typename T>
void NLPBuilder::addLUBC(double lb, double ub, T  expr) {
    g = vertcat(g, expr);
    lbg = vertcat(lbg, lb);
    ubg = vertcat(ubg, ub);
}

template<typename T>
void NLPBuilder::addEQ(T expr) {
    g = vertcat(g, expr);
    lbg = vertcat(lbg, 0);
    ubg = vertcat(ubg, 0);

}

template<typename T>
void NLPBuilder::addLBC(double lb, T expr) {
    g = vertcat(g, expr);
    lbg = vertcat(lbg, lb);
    ubg = vertcat(ubg, inf);
}

template<typename T>
void NLPBuilder::addUBC(double ub, T expr) {
    g = vertcat(g, expr);
    lbg = vertcat(lbg, -inf);
    ubg = vertcat(ubg, ub);
}

/*
 * Simple test using the opti helpers from casadi.
 */
int main(){
    float sx = 20;
    float sy = 20;
    Slice all;  // Basic slicer
    BicycleDynamics agent;
    int N = 10;
    double dt = 0.1;
    NLPBuilder nlp(N, dt);
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
    auto sol = nlp.nlpsolver(DMDict{{"x0", 0}, {"p", vertcat(st, st2)}, {"lbg", 0}, {"ubg", 0}});
//    OptiSol sol = opti.solve();

//    auto sol_x1 = sol.value(x);
//    auto sol_x2 = sol.value(x2);
//    auto sol_u = sol.value(u);
//    auto sol_u2 = sol.value(u2);
//    sf::ContextSettings settings( 0, 0, 8);
//    auto window = sf::RenderWindow(sf::VideoMode (1000, 1000), "blah",
//                                   sf::Style::None, settings);
//    auto view = sf::View(sf::Vector2f(sx/2, sy/2), sf::Vector2f(sx, -sy));
//    window.setView(view);
//    auto start = high_resolution_clock::now();
//    auto circ = sf::CircleShape(0.2);
//    auto circ2 = sf::CircleShape(0.2);
//    float t = 0;
//    window.setPosition(sf::Vector2i(0, 0));
//    e_u u_rt;
//    e_u u_rt2;
//    float t_last = 0.0f;
//    t = 0.0f;
//    float t_last_control = 0.0f;
//    int counter = 0;
//    while(true){
//        window.clear(sf::Color::Black);
//        sf::Event ev;
//        bool a = false;
//        while (window.pollEvent(ev)) {
//            if (ev.type == sf::Event::MouseButtonPressed) {
//                a = true;
//            }
//        }
//        auto stop = high_resolution_clock::now();
//        auto duration = duration_cast<milliseconds>(stop - start);
////        t = (float)duration.count()/2000.0f;
//        if(t - t_last > 0.01) {
//            int id = floor((t - t_last_control) / dt);
//            u_rt(0, 0) = sol_u(0, 0)->data()[0];
//            u_rt(1, 0) = sol_u(1, 0)->data()[0];
//
//            u_rt2(0, 0) = sol_u2(0, 0)->data()[0];
//            u_rt2(1, 0) = sol_u2(1, 0)->data()[0];
//            state += agent.getDeriv(state, u_rt) *dt;
//            state2 += agent.getDeriv(state2, u_rt2) *dt;
//            t_last = t;
//            circ.setPosition(sf::Vector2f(state(0, 0), state(1, 0)));
//            circ2.setPosition(sf::Vector2f(state2(0, 0), state2(1, 0)));
//            std::cout << "x: " << state(0,0) << ", y: " << state(1, 0) << " " << state(3, 0) << "\n";
//            window.draw(circ);
//            window.draw(circ2);
//            window.display();
//        }
//            std::cout << counter << std::endl;
//            t_last_control = t;
//            std::memcpy(st.ptr(), state.data(), sizeof(double)*4*1);
//            std::memcpy(st2.ptr(), state2.data(), sizeof(double)*4*1);
//            opti.set_value(x_init, st);
//            opti.set_value(x2_init, st2);
//            sol = opti.solve();
//            sol_x1 = sol.value(x);
//            sol_x2 = sol.value(x2);
//            sol_u = sol.value(u);
//            sol_u2 = sol.value(u2);
//
//        t += dt;
//        counter += 1;
//    }

}