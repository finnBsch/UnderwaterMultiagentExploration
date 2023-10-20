//
// Created by finn on 7/12/23.
//

#ifndef MASTERMAIN_LOGGER_H
#define MASTERMAIN_LOGGER_H

#include "string"
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <chrono>
#include "rrt.h"

inline bool exists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}


class LoggableType {
protected:
    std::string name;
    std::string data_str;
    std::string info;
    std::string column_names;
public:
    ~LoggableType(){
        file.close();
    }
    std::fstream file;
    std::string& getName(){
        return name;
    };
    virtual std::string& getInfo() = 0;
    std::string& getColumnNames(){
        return column_names;
    }
    virtual std::string& genDataStr() = 0;
};

class LoggablePoint: public LoggableType {
private:
    const double* x;
    const double* y;
public:
    std::string &getInfo() override;

    LoggablePoint(std::string name, const double* x, const double* y);
    std::string& genDataStr();
};

class LoggableEigen: public LoggableType {
private:
    const Eigen::MatrixXd* data;
public:
    std::string &getInfo() override;

    LoggableEigen(std::string name, const Eigen::MatrixXd* data);
    std::string& genDataStr();
};


class LoggableDynamicEigen: public LoggableType {
private:
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* data;
    int offset_row = 0;
    int offset_col = 0;
public:
    std::string &getInfo() override;

    LoggableDynamicEigen(std::string name, const Eigen::Matrix<double,
            Eigen::Dynamic, Eigen::Dynamic>* data, int offset_row = 0, int offset_col = 0);
    std::string& genDataStr();
};



class LoggableEigenVector: public LoggableType {
private:
    const Eigen::Matrix<double, Eigen::Dynamic, 1> * data;
public:
    std::string &getInfo() override;

    LoggableEigenVector(std::string name, const Eigen::Matrix<double,
            Eigen::Dynamic, 1> * data);
    std::string& genDataStr();
};

class LoggableTimer: public LoggableType {
private:
    std::unordered_map<std::string, double> times;
    std::unordered_map<std::string,
    std::chrono::high_resolution_clock::time_point> started;
public:
    LoggableTimer();

    std::string &getInfo() override;

    void addModuleToLog(std::string module_name);
    void startTimer(std::string module_name);
    void stopTimer(std::string module_name);
    std::string& genDataStr();
};

class SimLogger {
private:
    bool initialized = false;
    std::string base_path;
    std::string metrics_header = "t";
    std::fstream file; // for metrics
    std::vector<LoggableType*> logables;
    std::vector<LoggableType*> metrics;
public:
    SimLogger(std::string base_path);
    void addLoggable(LoggableType* logable);
    void addMetric(LoggableType* metric);
    void log(double t);
    void finalize();
    ~SimLogger();
};

#endif //MASTERMAIN_LOGGER_H
