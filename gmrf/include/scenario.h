//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//

#ifndef GMRF_SCENARIO_H
#define GMRF_SCENARIO_H

#include <array>
#include <string>
#include <eigen3/Eigen/Dense>
#include <unordered_map>
#include <SFML/Graphics.hpp>
#include <random>

using namespace Eigen;
class Objects;

class SourceObject {
private:
public:
    virtual ~SourceObject() = 0;
    virtual double getDistance(double x, double y) const = 0;
};

class LineSource: public SourceObject {
private:
    double x0;
    double y0;
    double x1;
    double y1;
public:
    LineSource(double x0, double y0, double x1, double y1);
    double getDistance(double x, double y) const;
};


class Scenario {
protected:
    double sx;
    double sy;
    std::vector<SourceObject*> source_obj;
public:
    virtual void getFieldValue(double x, double y, double& v, double stddev=0) = 0;
    virtual double getFieldValue(double x, double y, double stddev=0) = 0;
    virtual void getFlowValue(double x, double y, double& flow_x, double& flow_y, double stddev=0) = 0;
    virtual  Objects* getObjects() = 0;
    double getDistanceToSource(double x, double y) const;
};

struct BasicField{
    int n_x;
    int n_y;
    double step_size_x;
    double step_size_y;
    bool uniform;
    bool interpolate;
    Eigen::Matrix<double, Eigen::Dynamic, 1> x_values;
    Eigen::Matrix<double, Eigen::Dynamic, 1> y_values;

};

struct Field: BasicField {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> field_value;
};

struct FlowField: BasicField {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> field_value_x;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> field_value_y;
};

enum class wallType{
    wall = 0,
    inlet = 1,
    outlet = 2,
    source = 3
};


class Objects{
private:
public:
    std::unordered_map<wallType, std::vector<sf::Vertex>> objects;
    Objects();
    void addObject(wallType _type, const Eigen::Matrix<double, 1, 4>& verts);
};

class ConfigFileScenario: public Scenario {
private:
    std::mt19937 gen;
    bool has_obstacles;
    FlowField flow_field;
    Field field;
    std::string name;
    Objects objects;
public:
    const double& get_sx() const;
    const double& get_sy() const;
    ConfigFileScenario(std::string name);
    void getFieldValue(double x, double y, double &v, double stddev=0);
    double getFieldValue(double x, double y, double stddev=0);
    void getFlowValue(double x, double y, double &flow_x, double& flow_y, double stddev=0);
    Objects* getObjects();
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* getFlowFieldX();
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* getFlowFieldY();
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* getField();
    const std::string& getName() const;
    ~ConfigFileScenario();


};


#endif //GMRF_SCENARIO_H
