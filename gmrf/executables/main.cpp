//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//
#include "viz.h"
#include "gmrf/gmrf.h"
#include <chrono>
#include <iostream>
#include <random>
#include <scenario.h>
using namespace std::chrono;


double tanh_custom(double x, double up, double down, double scale){
    return tanhf((x-up)/scale) + 1 - (tanhf((x-down)/scale) +1);
}

double myMap(double x, double y, double t, Gmrf_Params & params) {
//    return (sinf(x ) * sinf(y)) + (powf(x-10.0/2,1)/10+powf(y-10.0/2, 1)/10);
//    return sinf(powf(x, 2)/2)* sinf(y/2);
//    return sinf(x)* sinf(y);
//    if(x < 2.63*2 && y < 7.5){
//        return tanh_custom(x, 10.0 / 4, 10.0 * 3 / 4, 1.0) * tanh_custom(y, 10.0/4, 10.0*3/4 , 1.0) + (7.5 - y);
//    }
//    return tanh_custom(x, 10.0 / 4, 10.0 * 3 / 4, 1.0) * tanh_custom(y, 10.0/4, 10.0*3/4 , 1.0);
//    if(x < 10.0/(N_X - 1) * 5 && y < 7.5){
//        return (2.0)- 1.0;
//    }
//    return 2.0*(sinf(2.0*M_PI/20*t)+1.0)/2.0-1.0;
    if(x < 10.0/(params.N_X - 1) * 5 && y < 7.5){
        return 0.5 + 2.0;
    }
    return 0.5;
}

Eigen::Matrix<double, Dynamic, Dynamic> true_field;
Eigen::Matrix<double, Dynamic, Dynamic> true_flow_x;
Eigen::Matrix<double, Dynamic, Dynamic> true_flow_y;
int n_meas = 20;
int main(int argc, char *argv[]){
    auto map_obj = ConfigFileScenario("test_map_2");


    double t;
    Gmrf_Params params;
    Gmrf_Params params2;
    Viz_Params viz_p;
    GMRF gmrf(params);
    true_field.resize(MAP_X, MAP_Y);
    true_flow_x.resize(params.N_X, params.N_Y);
    true_flow_y.resize(params.N_X, params.N_Y);
    RealtimeViz a(&gmrf, &true_field, &true_flow_x, &true_flow_y, viz_p, params);
    gmrf.addToViz(a, 0, 4, 0, 5);
    auto start = high_resolution_clock::now();
    gmrf.addWall();
    for(int i = 1; i < n_meas; i++){
        for(int j = 1; j < n_meas; j++) {
            auto temp = high_resolution_clock::now();
            t = (double)(duration_cast<milliseconds>(temp - start).count()) / 1000.0;
            double x = (double)i * params.size_x / ((double)n_meas + 0.01);
            double y = (double)j * params.size_y / ((double)n_meas + 0.01);
            gmrf.addMeasurement(x, y, map_obj.getFieldValue(x, y), t, field);
            double flow_x_meas;
            double flow_y_meas;
            map_obj.getFlowValue(x, y, flow_x_meas, flow_y_meas);
            gmrf.addMeasurement(x, y, flow_x_meas, t, flow_x);
            gmrf.addMeasurement(x, y, flow_y_meas, t, flow_y);
            gmrf.updateEstimates(t);
            for (int i_ = 0; i_ < MAP_X; i_++) {
                for (int j_ = 0; j_ < MAP_Y; j_++) {
                    true_field(i_, j_) = map_obj.getFieldValue(
                            params.size_x / MAP_X * (double) i_,
                            params.size_y / MAP_Y * (double) j_);
                }
            }
            for (int i_ = 0; i_ < params.N_X; i_++) {
                for (int j_ = 0; j_ < params.N_Y; j_++) {
                    map_obj.getFlowValue(
                            params.size_x / (double) params.N_X * (double) i_,
                            params.size_y / (double) params.N_Y * (double) j_,
                            true_flow_x.coeffRef(i_, j_),
                            true_flow_y.coeffRef(i_, j_));
                }
            }
            a.draw();
        }
    }
    std::cout << t << std::endl;
    while(true){
        a.draw();
    }
//    while(true){
//        auto temp = high_resolution_clock::now();
//        t = (double)(duration_cast<milliseconds>(temp - start).count()) / 1000.0;
//        gmrf.updateEstimates(t);
////        a.draw();
//    }
}