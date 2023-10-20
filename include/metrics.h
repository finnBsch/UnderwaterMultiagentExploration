//
// Created by finn on 7/13/23.
//

#ifndef MASTERMAIN_METRICS_H
#define MASTERMAIN_METRICS_H
#include "logger.h"
#include "gmrf/gmrf.h"

class Metrics : public LoggableType {
protected:
    double metric = 0;
public:
    ~Metrics() {};
    virtual double updateMetric() = 0;

    std::string &getInfo() override;

    std::string &genDataStr() override;
};

class FieldLikelihoodMetric : public Metrics {
private:
    GMRF* gmrf;
    double resolution; // Every 10 cm.
    double size_x;
    double size_y;
    int n_x;
    int n_y;
    ConfigFileScenario* scenario;

public:
    FieldLikelihoodMetric(GMRF* gmrf, ConfigFileScenario* scenario, double
        resolution = 0.1);
    double updateMetric() override;
};

class FieldLogLikelihoodMetric : public Metrics {
private:
    GMRF* gmrf;
    double resolution; // Every 10 cm.
    double size_x;
    double size_y;
    int n_x;
    int n_y;
    ConfigFileScenario* scenario;
public:
    FieldLogLikelihoodMetric(GMRF* gmrf, ConfigFileScenario* scenario, double
        resolution = 0.1);
    double updateMetric() override;
};

class RMSEMetric : public Metrics {
private:
    GMRF* gmrf;
    double resolution; // Every 10 cm.
    double size_x;
    double size_y;
    int n_x;
    int n_y;
    ConfigFileScenario* scenario;
    std::vector<RRTObstacle*> obstacles;
public:
    RMSEMetric(GMRF* gmrf, ConfigFileScenario* scenario, std::vector<RRTObstacle*> obstacles, double resolution = 0.1);
    double updateMetric() override;
};


class SourceDistMetric : public Metrics {
private:
    GMRF* gmrf;
    ConfigFileScenario* scenario;
public:
    SourceDistMetric(GMRF* gmrf, ConfigFileScenario* scenario);
    double updateMetric() override;
};
#endif //MASTERMAIN_METRICS_H
