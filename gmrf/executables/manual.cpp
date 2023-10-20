//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//
#include "gmrf_viz.h"
#include "gmrf/gmrf.h"
#include <chrono>
#include <iostream>
#include <random>
#include <scenario.h>


using namespace std::chrono;


Eigen::Matrix<double, Dynamic, Dynamic> true_field;
Eigen::Matrix<double, Dynamic, Dynamic> true_flow_x;
Eigen::Matrix<double, Dynamic, Dynamic> true_flow_y;

int main(int argc, char *argv[]){
    auto map_obj = ConfigFileScenario("CaseB");
    std::string log_base_path = "/home/finn/CLionProjects/mastersthesis_main"
                           "/Experiments/Paper"
                           "/CaseB/";

    true_field.resize(MAP_X, MAP_Y);
    double t = 0;
    GmrfParams params;
    params.size_x = map_obj.get_sx();
    params.size_y = map_obj.get_sy();
    GmrfParams params2;
    params2.size_x = map_obj.get_sx();
    params2.size_y = map_obj.get_sy();
    Viz_Params viz_p;
    viz_p.offset_y = 100;
    Viz_Params viz_p2;
    viz_p2.offset_y = 100;
    viz_p2.offset_x = viz_p.res_x * 1.1;
    params2.flow_field = false;
    GMRF gmrf(params, &map_obj);
    GMRF gmrf2(params2, &map_obj);
    true_field.resize(MAP_X, MAP_Y);
    true_flow_x.resize(params.N_X, params.N_Y);
    true_flow_y.resize(params.N_X, params.N_Y);
    RealtimeViz a(&gmrf, &true_field, &true_flow_x, &true_flow_y, viz_p, params);
    RealtimeViz a2(&gmrf2, &true_field, &true_flow_x, &true_flow_y, viz_p2, params2);
    ObjectsViz obj_viz = ObjectsViz(map_obj.getObjects());
    a.addVizObject(&obj_viz, map_types::std_map);
    gmrf.addToViz(a, 0, 4, 0, 5);
    gmrf2.addToViz(a2, 0, 4, 0, 5);
    auto start = high_resolution_clock::now();
    int i = 0;
//    std::vector<double> x_pos {5, 5, 5, 5, 4, 3, 2, 1, 1, 1, 1};
//    std::vector<double> y_pos {1, 3, 5, 8, 9, 9, 8, 7, 6, 3, 1};
//    std::vector<double> x_pos_2 {5, 6.5, 8, 9, 9, 9};
//    std::vector<double> y_pos_2 {1, 1, 2, 2, 4, 7, 9};
//    for(int j = 0; j < x_pos.size(); j++){
//        double x = x_pos[j]/10 * params.size_x;
//        double y = y_pos[j]/10 * params.size_y;
//        gmrf.addMeasurement(x, y, map_obj.getFieldValue(x, y), t, field);
//        gmrf.updateEstimates(t);
//    }
//    for(int j = 0; j < x_pos_2.size(); j++){
//        double x = x_pos_2[j]/10 * params.size_x;
//        double y = y_pos_2[j]/10 * params.size_y;
//        gmrf.addMeasurement(x, y, map_obj.getFieldValue(x, y), t, field);
//        gmrf.updateEstimates(t);
//    }
    while(true){
        auto temp = high_resolution_clock::now();
        t = (double)(duration_cast<milliseconds>(temp - start).count()) / 1000.0f;
        if(a.meas_ev.newMeasurement) {
            i += 1;
            a.meas_ev.newMeasurement = false;
            double x = a.meas_ev.measurement_location[0];
            double y = a.meas_ev.measurement_location[1];
            auto iter_start = high_resolution_clock::now();
            gmrf.addMeasurement(x, y, map_obj.getFieldValue(x, y), t, field);
            gmrf2.addMeasurement(x, y, map_obj.getFieldValue(x, y), t, field);
            double flow_x_meas;
            double flow_y_meas;
            map_obj.getFlowValue(x, y, flow_x_meas, flow_y_meas);
            gmrf.addMeasurement(x, y, flow_x_meas, t, flow_x);
            gmrf2.addMeasurement(x, y, flow_x_meas, t, flow_x);
            gmrf.addMeasurement(x, y, flow_y_meas, t, flow_y);
            gmrf2.addMeasurement(x, y, flow_y_meas, t, flow_y);
            gmrf.updateEstimates(t);
//            gmrf2.updateEstimates(t);
            double x_; double y_; bool success;
            gmrf.sourceHypothesisNeighbor(x_, y_, success);
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - iter_start);
             std::cout << "(x : " << x << ", y: " << y << ") - " << map_obj.getFieldValue(x, y) << std::endl;
            std::cout << "Duration with " << i << " measurements: " << duration.count() << ". Current time passed is "
                      << t << "s." << std::endl;
            if (i>5){
            }
        }
        else{
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

        }

//        a.set_data(2, &true_field, 0.0f, 2.0f);
        a.draw();
        a2.draw();
    }
}