//
// Created by finn on 7/13/23.
//
#include "metrics.h"

FieldLikelihoodMetric::FieldLikelihoodMetric(GMRF *gmrf,
                                             ConfigFileScenario *scenario,
                                             double resolution) {
    name = "FieldLikelihoodMetric";
    this->gmrf = gmrf;
    this->scenario = scenario;
    this->resolution = resolution;
    size_x = scenario->get_sx();
    size_y = scenario->get_sy();
    n_x = (int) floor(size_x / resolution);
    n_y = (int) floor(size_y / resolution);
    this->info = "";
    this->info += std::to_string(n_x);
    this->info += ",";
    this->info += std::to_string(n_y);
}

double FieldLikelihoodMetric::updateMetric() {
    metric = 0;
    int counter = 0;
    for(int i = 0; i < n_x; i++){
        for(int j = 0; j < n_y; j++){

            double x = i * resolution;
            double y = j * resolution;
            if(gmrf->isReachable(x, y)) {
                counter += 1;
                double field_value = scenario->getFieldValue(x, y);
                double field_likelihood = gmrf->getLikelihood(x, y,
                                                              field_value);
                metric += field_likelihood;
            }
        }
    }
    return metric/counter;
}

std::string &Metrics::genDataStr() {
    data_str = std::to_string(updateMetric());
    return data_str;
}

std::string &Metrics::getInfo() {
    return info;
}

FieldLogLikelihoodMetric::FieldLogLikelihoodMetric(GMRF *gmrf,
                                                   ConfigFileScenario *scenario,
                                                   double resolution) {
    name = "FieldLogLikelihoodMetric";
    this->gmrf = gmrf;
    this->scenario = scenario;
    this->resolution = resolution;
    size_x = scenario->get_sx();
    size_y = scenario->get_sy();
    n_x = (int) floor(size_x / resolution);
    n_y = (int) floor(size_y / resolution);
    this->info = "";
    this->info += std::to_string(n_x);
    this->info += ",";
    this->info += std::to_string(n_y);
}

double FieldLogLikelihoodMetric::updateMetric() {
    metric = 0;
    int counter = 0;
    for(int i = 0; i < n_x; i++){
        for(int j = 0; j < n_y; j++){
            double x = i * resolution;
            double y = j * resolution;
            if(gmrf->isReachable(x, y)) {
                counter += 1;
                double field_value = scenario->getFieldValue(x, y);
                double field_likelihood = gmrf->getLogLikelihood(x, y,
                                                                 field_value);
                metric += field_likelihood;
            }
        }
    }
    return metric;
}

RMSEMetric::RMSEMetric(GMRF *gmrf, ConfigFileScenario *scenario, std::vector<RRTObstacle*> obstacles,
                       double resolution) {
    this->obstacles = obstacles;
    name = "RMSEMetric";
    this->gmrf = gmrf;
    this->scenario = scenario;
    this->resolution = resolution;
    size_x = scenario->get_sx();
    size_y = scenario->get_sy();
    n_x = (int) floor(size_x / resolution);
    n_y = (int) floor(size_y / resolution);
    this->info = "";
    this->info += std::to_string(n_x);
    this->info += ",";
    this->info += std::to_string(n_y);

}

double RMSEMetric::updateMetric() {
    metric = 0;
    int counter = 0;
//    double max = 0;
//    double max_x;
//    double max_y;
//


//    std::cout << "[";
    for(int i = 0; i < n_x; i++){
//        std::cout << "[";
        for(int j = 0; j < n_y; j++){
            double x = i * resolution;
            double y = j * resolution;
            // Check if point is closer then 0.1 to an obstacle
            bool close_to_obstacle = false;
            for (auto &obstacle : obstacles) {
                double dist = obstacle->getDistance(x, y);
                if (dist < 0.1) {
                    close_to_obstacle = true;
                    break;
                }
            }
            if(gmrf->isReachable(x, y) && !close_to_obstacle or true) {
                counter += 1;
                double field_value = scenario->getFieldValue(x, y);
                double estimate = gmrf->getEstimate(x, y);
                metric += pow(field_value - estimate, 2);
//                std::cout << pow(field_value - estimate, 2) << ", ";
//                double error_v = pow(field_value - estimate, 2);
//                if (error_v > max) {
//                    max = error_v;
//                    max_x = x;
//                    max_y = y;
//                }
            }
            else{
//                std::cout << "0, ";
            }
        }
//        std::cout << "],";
    }
//    std::cout << "]" << std::endl;
//    std::cout << "Max error is: " << max << " at: " << max_x << ", " << max_y << std::endl;
    return sqrt(metric/counter);
}

SourceDistMetric::SourceDistMetric(GMRF *gmrf, ConfigFileScenario *scenario) {
    name = "SourceDistMetric";
    this->gmrf = gmrf;
    this->scenario = scenario;
    this->info = "";
}

double SourceDistMetric::updateMetric() {
    if(gmrf->source_found) {
        return scenario->getDistanceToSource(gmrf->source_x, gmrf->source_y);
    }
    else{
        // Return max distance in the field
        return sqrt(pow(scenario->get_sx(), 2) + pow(scenario->get_sy(), 2));
    }
}
