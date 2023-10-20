//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//


#include "scenario.h"
#include <iostream>
#include <toml++/toml.h>
#include "res_loc.h"

template<typename T> T read_toml_value(toml::node_view<toml::node> node_view){
    auto optional_v = node_view.value<T>();
    if(optional_v){
        return optional_v.value();
    }
    else{
        std::cerr << "TOML field not populated." << std::endl;
        return T();
    }
}

template<typename T> T read_toml_value(toml::node& node){
    auto optional_v = node.value<T>();
    if(optional_v){
        return optional_v.value();
    }
    else{
        std::cerr << "TOML field not populated." << std::endl;
        return T();
    }
}

void read_toml_matrix(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& mat, toml::node_view<toml::node> node_view){
    int x = (int)node_view.as_array()->size();
    int y = (int)node_view.as_array()->operator[](0).as_array()->size();
    mat.resize(x, y);
    auto arr = node_view.as_array();
    for(int i = 0; i < x; i++){
        auto arr_x = arr->operator[](i).as_array();
        for(int j = 0; j < y; j++){
            mat(i, j) = read_toml_value<double>(arr_x->operator[](j));
        }
    }
}

void read_toml_matrix(Eigen::Matrix<double, Eigen::Dynamic, 1>& mat, toml::node_view<toml::node> node_view){
    int x = (int)node_view.as_array()->size();
    mat.resize(x, 1);
    auto arr = node_view.as_array();
    for(int i = 0; i < x; i++){
        mat(i, 0) = read_toml_value<double>(arr->operator[](i));
    }
}

void uniform_initialize(Eigen::Matrix<double, Eigen::Dynamic, 1>& mat, int num_pts, double size_dim){
    double step_size = size_dim/((double)num_pts - 1.0f);
    mat.resize(num_pts, 1);
    for(int i = 0; i < num_pts; i++){
        mat(i, 0) = (double) i * step_size;
    }
}

ConfigFileScenario::ConfigFileScenario(std::string name):
    gen(std::random_device()()){
    this->name = name;
    toml::table tbl;
    try {
        tbl = toml::parse_file(RES_LOC "/configs/maps/" + name + ".toml");
    }
    catch(const toml::parse_error& err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
    }
    has_obstacles = read_toml_value<bool>(tbl["general"]["has_obstacles"]);
    sx = read_toml_value<double>(tbl["general"]["size_x"]);
    sy = read_toml_value<double>(tbl["general"]["size_y"]);

    // Field
    field.uniform = read_toml_value<bool>(tbl["value_field"]["uniform_discretization"]);
    field.interpolate = read_toml_value<bool>(tbl["value_field"]["interpolate"]);
    read_toml_matrix(field.field_value, tbl["value_field"]["values"]);
    field.n_x = (int)field.field_value.rows();
    field.n_y = (int)field.field_value.cols();
    field.step_size_x = sx/((double)field.n_x - 1.0f);
    field.step_size_y = sy/((double)field.n_y - 1.0f);
    if(!field.uniform){
        // TODO ..
        read_toml_matrix(field.x_values, tbl["value_field"]["x_values"]);
        read_toml_matrix(field.y_values, tbl["value_field"]["y_values"]);
    }
    else{
        uniform_initialize(field.x_values, field.n_x, sx);
        uniform_initialize(field.y_values, field.n_y, sy);
    }
    // Flow Field
    flow_field.uniform = read_toml_value<bool>(tbl["flow_field"]["uniform_discretization"]);
    flow_field.interpolate = read_toml_value<bool>(tbl["flow_field"]["interpolate"]);
    read_toml_matrix(flow_field.field_value_x, tbl["flow_field"]["values_x"]);
    read_toml_matrix(flow_field.field_value_y, tbl["flow_field"]["values_y"]);
    flow_field.n_x = (int)flow_field.field_value_x.rows();
    flow_field.n_y = (int)flow_field.field_value_y.cols();
    flow_field.step_size_x = sx/((double)flow_field.n_x - 1.0f);
    flow_field.step_size_y = sy/((double)flow_field.n_y - 1.0f);
    Eigen::Matrix<double, Dynamic, Dynamic> walls;
    int num_walls = read_toml_value<int>(tbl["objects"]["walls"]["num_walls"]);
    if(num_walls > 0) {
        read_toml_matrix(walls, tbl["objects"]["walls"]["wall_vertices"]);
        for (int i = 0; i < walls.rows(); i++) {
            objects.addObject(wallType::wall, walls(i, all));
        }
    }
    Eigen::Matrix<double, Dynamic, Dynamic> inlet_;
    num_walls = read_toml_value<int>(tbl["objects"]["inlet"]["num_walls"]);
    if(num_walls > 0) {
        read_toml_matrix(inlet_, tbl["objects"]["inlet"]["wall_vertices"]);
        for (int i = 0; i < inlet_.rows(); i++) {
            objects.addObject(wallType::inlet, inlet_(i, all));
        }
    }
    Eigen::Matrix<double, Dynamic, Dynamic> outlet_;
    num_walls = read_toml_value<int>(tbl["objects"]["outlet"]["num_walls"]);
    if(num_walls > 0) {
        read_toml_matrix(outlet_, tbl["objects"]["outlet"]["wall_vertices"]);
        for (int i = 0; i < outlet_.rows(); i++) {
            objects.addObject(wallType::outlet, outlet_(i, all));
        }
    }
    Eigen::Matrix<double, Dynamic, Dynamic> source_;
    num_walls = read_toml_value<int>(tbl["objects"]["source"]["num_walls"]);
    if (num_walls > 0) {
        read_toml_matrix(source_, tbl["objects"]["source"]["wall_vertices"]);
        for (int i = 0; i < source_.rows(); i++) {
            objects.addObject(wallType::source, source_(i, all));
            source_obj.push_back(new LineSource(source_(i, 0), source_(i, 1), source_(i, 2), source_(i, 3)));
//            std::cout << "source obj: " << source_(i, 0) << " " << source_(i, 1) << " " << source_(i, 2) << " " << source_(i, 3) << std::endl;
        }
    }
    if(!flow_field.uniform){
        // TODO ..
    }
    else{
        uniform_initialize(flow_field.x_values, flow_field.n_x, sx);
        uniform_initialize(flow_field.y_values, flow_field.n_y, sy);
    }
}

void ConfigFileScenario::getFieldValue(double x, double y, double &v, double stddev) {
    if(!field.interpolate){
        int i = (int)roundf(x/field.step_size_x);
        int j = (int)roundf(y/field.step_size_y);
        v = field.field_value(i, j);
        if (stddev!=0){
            std::normal_distribution<> distribution(v, stddev);
            v = distribution(gen);
        }
    }
    else{
        // TODO ..
    }
}

void ConfigFileScenario::getFlowValue(double x, double y, double &flow_x, double &flow_y, double stddev) {
    if(!flow_field.interpolate){
        int i = (int)roundf(x/flow_field.step_size_x);
        int j = (int)roundf(y/flow_field.step_size_y);
        flow_x = flow_field.field_value_x(i, j);
        flow_y = flow_field.field_value_y(i, j);
        if (stddev!=0){
            std::normal_distribution<> distribution(flow_x, stddev);
            flow_x = distribution(gen);
            distribution = std::normal_distribution<>(flow_y, stddev);
            flow_y = distribution(gen);
        }
    }
    else{
        // TODO ..
    }
}

double ConfigFileScenario::getFieldValue(double x, double y, double stddev) {
    if (!field.interpolate) {
        int i = (int) roundf(x / field.step_size_x);
        int j = (int) roundf(y / field.step_size_y);
        if(stddev != 0){
            std::normal_distribution<> distribution(field.field_value(i, j), stddev);
            return distribution(gen);
        }
        else {
            return field.field_value(i, j);
        }
    } else {
        return 0;
    }
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *ConfigFileScenario::getFlowFieldX() {
    return &flow_field.field_value_x;
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *ConfigFileScenario::getFlowFieldY() {
    return &flow_field.field_value_y;
}

Objects *ConfigFileScenario::getObjects() {
    return &objects;
}

const double &ConfigFileScenario::get_sx() const {
    return sx;
}

const double &ConfigFileScenario::get_sy() const {
    return sy;
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *
ConfigFileScenario::getField() {
    return &field.field_value;
}

const std::string &ConfigFileScenario::getName() const {
    return name;
}

ConfigFileScenario::~ConfigFileScenario() {
    for (auto &i : source_obj) {
        delete i;
    }
}

SourceObject::~SourceObject() {
}

void Objects::addObject(wallType _type, const Eigen::Matrix<double, 1, 4>& verts) {
    objects[_type].emplace_back(sf::Vector2f(verts(0, 0), verts(0, 1)), sf::Color::Magenta);
    objects[_type].emplace_back(sf::Vector2f(verts(0, 2), verts(0, 3)), sf::Color::Magenta);
}

Objects::Objects() {
//    for(int i = 0; i < 4; i++){
//        objects[static_cast<wallType>(i)] = std::vector<sf::Vertex>();
//    }
}

LineSource::LineSource(double x0, double y0, double x1, double y1) {
    this->x0 = x0;
    this->y0 = y0;
    this->x1 = x1;
    this->y1 = y1;
}
/**
 * Get min distance to the line segment.
 * @param x
 * @param y
 * @return
 */
double LineSource::getDistance(double x, double y) const {
    double dx = x1 - x0;
    double dy = y1 - y0;
    double d = dx*dx + dy*dy;
    double u = ((x - x0) * dx + (y - y0) * dy) / d;
    if (u > 1) {
        u = 1;
    } else if (u < 0) {
        u = 0;
    }
    double x_ = x0 + u * dx;
    double y_ = y0 + u * dy;
    dx = x_ - x;
    dy = y_ - y;
    return sqrt(dx*dx + dy*dy);
}

double Scenario::getDistanceToSource(double x, double y) const {
    double dist = std::numeric_limits<double>::max();
    for (auto &source : source_obj) {
        double d = source->getDistance(x, y);
        if (d < dist) {
            dist = d;
        }
    }
    return dist;
}
