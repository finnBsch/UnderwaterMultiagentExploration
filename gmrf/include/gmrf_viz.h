//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//

#ifndef GMRF_VIZ_H
#define GMRF_VIZ_H

#include <SFML/Graphics.hpp>
#include <eigen3/Eigen/Dense>
#include "settings.h"
#include "color_mixer.h"
#include "gmrf_configs.h"
#include "scenario.h"
#include "colormaps.h"

class GMRF;

struct MeasurementEvent {
    bool newMeasurement = false;
    std::array<double, 2> measurement_location;
};

enum class map_types {
    field_map = 0,
    std_map = 1,
    gt_map = 2,
    f_field_map = 3,
    f_std_map = 4,
    f_gt_map = 5
};

/**
 * Generic drawable object
 */
class DrawableObject : public sf::Drawable, public sf::Transformable {
public:
    virtual void
    draw(sf::RenderTarget &target, sf::RenderStates states) const = 0;

    virtual void prepareViz() = 0;
};


class InfoBox : public sf::Drawable, public sf::Transformable {
public:
    explicit InfoBox(int charactersize);

    void update(double x, double y, double v);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

private:
    int charSize;
    sf::Font font;
    sf::Text text;
    sf::RectangleShape box;
};

class GuiText : public DrawableObject {
public:
    void setTime(double t);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

    GuiText(int charactersize);

    void prepareViz() {

    };
private:
    int char_size;
    double t = 0.0f;
    sf::Font font;
    sf::Text text;
};


class ObjectsViz : public DrawableObject {
private:
    Objects *objects;
    std::vector<sf::RectangleShape> viz_objects;
public:
    void prepareViz() {};

    explicit ObjectsViz(Objects *objects);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};

class View : public DrawableObject {
protected:
    double sx;
    double sy;
public:
    View(double sx, double sy);

    sf::View view;
    std::vector<DrawableObject *> viz_objects;

    virtual void updateData() = 0;

    virtual double getValue(double x, double y) = 0;

};

class ColorMap : public View {
private:
    sf::Uint8 *pixels = new sf::Uint8[MAP_X * MAP_Y * 4];
//    std::array<int, 3> col0{102, 255, 71};
//    std::array<int, 3> col1{0, 102, 255};
    std::array<std::array<double, 3>, 256> colormap;
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *data_map;
public:
    ColorMap(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *data,
             double sx, double sy, std::string cmap);

    void prepareViz() override {};

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

    void updateData() override;

    double getValue(double x, double y);
};

class FlowMap : public View {
private:
    int n_x;
    int n_y;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> flow_x;
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> flow_y;
    std::vector<sf::RectangleShape> lines;
public:
    FlowMap(const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> *ptr_flow_x,
            const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> *ptr_flow_y,
            double sx, double sy, int n_x, int n_y);

    FlowMap(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *ptr_flow_x,
            const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *ptr_flow_y,
            double sx, double sy, int n_x, int n_y);

    void prepareViz() override {};

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;

    void updateData() override;

    double getValue(double x, double y) override;
};

class PrimitiveWrapper : public DrawableObject {
private:
    sf::Drawable *obj;
public:
    void prepareViz() override {};

    explicit PrimitiveWrapper(sf::Drawable *obj);

    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};


class RealtimeViz {
private:
    InfoBox info_box;
    double map_size_x;
    double map_size_y;
    Viz_Params viz_p;
    GmrfParams gmrf_p;
    sf::View gui_view;
    std::vector<DrawableObject *> gui_obs;
    sf::ContextSettings settings;
    sf::RenderWindow window;
    GuiText text_obj;
    std::unordered_map<map_types, View *> views;
    double t = 0.0f;
public:
    sf::View standard_view;

    sf::Vector2f mouse_pos;
    MeasurementEvent meas_ev;

    void addVizObject(DrawableObject *obj, map_types map_id);

    void addVizObject(sf::Drawable* obj, map_types map_id);

    RealtimeViz(GMRF *gmrf,
                const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *gt_ptr,
                const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *f_x_gt_ptr,
                const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *f_y_gt_ptr,
                Viz_Params viz_p_, GmrfParams gmrf_p);

    void setTime(double t);

    void draw();

    void update();

    ~RealtimeViz();
};

#endif //GMRF_VIZ_H
