//
// Created by finn on 7/12/23.
//
#include "logger.h"
#include <stdexcept>


LoggableEigen::LoggableEigen(std::string name, const Eigen::MatrixXd *data) {
    this->name = name;
    this->data = data;
    this->info = "Eigen\n";
    this->info += std::to_string(data->rows());
    this->info += ",";
    this->info += std::to_string(data->cols());
}

std::string &LoggableEigen::genDataStr() {
    data_str = "";
    for(int i = 0; i < data->size(); i++){
        data_str += std::to_string((*data)(i)) + ",";
    }
    return data_str;
}

std::string &LoggableEigen::getInfo() {
    return info;
}

SimLogger::SimLogger(std::string base_path) {
    this->base_path = base_path;
    std::string path_to_file = base_path  + "metrics.csv";
    if(exists(path_to_file)){
        throw std::runtime_error("Logging into existing file! Clear first!");
    }
    file.open(path_to_file, std::ios::out);
}

void SimLogger::addLoggable(LoggableType *logable) {
    logables.push_back(logable);
    std::string path_to_file = base_path + logable->getName() + ".csv";
    if(exists(path_to_file)){
        throw std::runtime_error("Logging into existing file! Clear first!");
    }
    logable->file.open(path_to_file, std::ios::out);
    logable->file << logable->getInfo() << std::endl;
}

void SimLogger::log(double t) {
    if(!initialized){
        initialized = true;
        file << metrics_header << std::endl;
    }
    for(auto logable: logables){
        logable->file << t << "," << logable->genDataStr() << std::endl;
    }
    file << t;
    for(auto metric: metrics){
        file << "," << metric->genDataStr();
    }
    file << std::endl;
}

void SimLogger::finalize() {
    for (auto logable: logables) {
        logable->file.close();
    }
    file.close();
}

void SimLogger::addMetric(LoggableType *metric) {
    metrics.push_back(metric);
    metrics_header += "," + metric->getName();


}

SimLogger::~SimLogger() {
    finalize();
//    for (auto logable: logables) {
//        delete logable;
//    }
//    for (auto metric: metrics) {
//        delete metric;
//    }
}

LoggableEigenVector::LoggableEigenVector(std::string name,
                                         const Eigen::Matrix<double,
                                                 Eigen::Dynamic, 1> *data) {
    this->name = name;
    this->data = data;
    this->info = "EigenVector\n";
    this->info += std::to_string(data->rows());
    this->info += ",";
    this->info += std::to_string(data->cols());
}

std::string &LoggableEigenVector::genDataStr() {
    data_str = "";
    for(int i = 0; i < data->size(); i++){
        data_str += std::to_string((*data)(i)) + ",";
    }
    return data_str;
}

std::string &LoggableEigenVector::getInfo() {
    return info;
}

LoggableTimer::LoggableTimer() {
    name = "timings";
}

void LoggableTimer::addModuleToLog(std::string module_name) {
    times.emplace(module_name, 0);
    started.emplace(module_name, std::chrono::high_resolution_clock::now());
}

void LoggableTimer::startTimer(std::string module_name) {
    started.at(module_name) = std::chrono::high_resolution_clock::now();

}

void LoggableTimer::stopTimer(std::string module_name) {
    times.at(module_name) =
            std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - started.at(module_name)
    ).count();
}

std::string &LoggableTimer::genDataStr() {
    data_str = "";
    for(auto &time: times){
        data_str += std::to_string(time.second) + ",";
    }
    // delete last comma
    data_str.pop_back();
    return data_str;
}

std::string &LoggableTimer::getInfo() {
    info = "t";
    for(auto &time: times){
        info += "," + time.first;
    }
    return info;
}

std::string &LoggableDynamicEigen::getInfo() {
    return info;
}

LoggableDynamicEigen::LoggableDynamicEigen(std::string name,
                                           const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *data,
                                           int offset_row, int offset_col) {
    this->offset_row = offset_row;
    this->offset_col = offset_col;
    this->name = name;
    this->data = data;
    this->info = "Eigen\n";
    this->info += "d";
    this->info += ",";
    this->info += "d";

}

std::string &LoggableDynamicEigen::genDataStr() {
    if(offset_row > 0 || offset_col > 0) {
        Eigen::MatrixXd mat = (*data)(Eigen::seq(offset_row, data->rows() - 1),
                                      Eigen::seq(offset_col, data->cols() - 1));
        data_str = std::to_string(mat.rows()) + "," + std::to_string(mat.cols()) + ",";
        for (int i = 0; i < mat.size(); i++) {
            data_str += std::to_string((mat)(i)) + ",";
        }
        return data_str;
    }
    else{
        data_str = std::to_string(data->rows()) + "," + std::to_string(data->cols()) + ",";
        for (int i = 0; i < data->size(); i++) {
            data_str += std::to_string((*data)(i)) + ",";
        }
        return data_str;
    }
}

std::string &LoggablePoint::getInfo() {
    info = "x,y";
    return info;
}

LoggablePoint::LoggablePoint(std::string name, const double *x, const double *y) {
    this->name = name;
    this->x = x;
    this->y = y;
}

std::string &LoggablePoint::genDataStr() {
    data_str = std::to_string(*x) + "," + std::to_string(*y);
    return data_str;
}
