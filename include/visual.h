//
// Created by finn on 10/7/23.
//

#ifndef MASTERMAIN_VISUAL_H
#define MASTERMAIN_VISUAL_H
#include <SFML/Graphics.hpp>
#include "gmrf/gmrf.h"
#include <libInterpolate/Interpolate.hpp>
#include "colormaps.h"
#include "linestrip.h"


enum class ColorMapMode {
    estimates,
    variance,
    source_predictions,
    ground_truth
};

struct VizConf {
    bool renderVideo = false;
    double speed = 3;
//    int renderWidth = 1920;
//    int renderHeight = 1080;
    int renderWidth = 4096;
    int renderHeight = 2304;
    std::string base_path = "/home/finn/mastersthesis_main/Experiments/Video/1";

    bool draw_gmrf_grid = true;
    bool draw_connectivity = true;
    bool close_up = false;
    ColorMapMode colorMapMode = ColorMapMode::source_predictions;
};

class MyDrawable : public sf::Drawable, public sf::Transformable {
public:
    virtual void updateData() = 0;

};

class ScenarioViz : public MyDrawable {
private:
    Objects *objects;
    std::vector<sf::RectangleShape> obstacles;
public:
    ScenarioViz(Objects* objects);
    void updateData(){};
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};


class CrossShape : public sf::Drawable, public sf::Transformable {
private:
    double x;
    double y;
    double size;
    double line_thickness;

public:
    CrossShape(double x, double y, double size, double line_thickness);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};

/**
 * Draw the GMRF grid as a grid of crosses
 */
class GmrfGrid : public MyDrawable {
private:
    std::vector<CrossShape> crosses;
    const GmrfParams *gmrf_params;
public:
    GmrfGrid(const  GmrfParams *gmrf_params);

    void updateData() override;

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};

/**
 * Agent viz as triangle
 */
class Agent: public MyDrawable {
private:
    double x = 0;
    double y = 0;
    double phi = 0;
    std::vector<float> x_coords;
    std::vector<float> y_coords;
    sf::ConvexShape agent_shape;
    LineStrip line_strip;
public:
    Agent();
    void setState(double x, double y, double phi);
    void updateData() override;
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    double get_x() const;
    double get_y() const;
};


class Connectivity: public MyDrawable{
private:
    std::vector<Agent*>* agents;
    std::vector<bool> connected;
    std::vector<DashedLineStrip> connectivity_lines;
    double max_dist = 5.0;
public:
    Connectivity(std::vector<Agent*>* agents);

protected:
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

public:
    void updateData() override;
};


class ValueMap : public MyDrawable {
private:
    std::array<std::array<double, 3>, 256> colormap;
//    _2D::BilinearInterpolator<double> interp;
    _2D::NearestInterpolator<double> interp;
    Eigen::Matrix<sf::Uint8, Eigen::Dynamic, 1> pixel_values;
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *values_ptr;
    double sx;
    double sy;
    int subsampling_fac = 20;
    double interpolate_value(double x, double y) const;
public:
    ValueMap(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
             *values_ptr, double sx, double sy, std::string colormap_name);

    void updateData() override;

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};


class Visual {
private:
    bool is_closed_up = false;
    float sx;
    float sy;
    int frame_counter = 0;
    VizConf conf;
    sf::RenderWindow window;
    sf::View view;

    std::vector<MyDrawable*> drawables;
    std::vector<Agent*> agents;
    std::vector<ValueMap*> value_maps;
    void exportFrame();
public:
    bool paused = false;

    Visual(VizConf conf, GMRF *gmrf,
           const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *gt_ptr,
           const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *est_ptr,
           const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *source_ptr,
           const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *variance_ptr,
           GmrfParams gmrf_params);

    void addDrawable(MyDrawable* drawable);
    void addAgent(Agent* agent);
    void draw();
};

void renderVideo(std::string path);
#endif //MASTERMAIN_VISUAL_H
