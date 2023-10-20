//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//

#include <iostream>
#include "gmrf_viz.h"
#include <iomanip>
#include "gmrf/gmrf.h"
#include "res_loc.h"

void RealtimeViz::draw() {
    window.clear(sf::Color::Black);
    sf::Event ev;
    while (window.pollEvent(ev)) {
        if(ev.type == sf::Event::MouseButtonPressed){
            for(auto v: views){
                if(v.second->view.getViewport().contains(((double)ev.mouseButton.x) / viz_p.res_x, ((double)ev.mouseButton.y) / viz_p.res_y)){
                    window.setView(v.second->view);
                    sf::Vector2f worldPos = window.mapPixelToCoords(sf::Vector2i(ev.mouseButton.x, ev.mouseButton.y));
                    meas_ev.newMeasurement = true;
                    meas_ev.measurement_location[0] = worldPos.x;
                    meas_ev.measurement_location[1] = worldPos.y;
                    break;
                }
            }
        }
    }
    bool in_view = false;
    for(auto v: views){
        auto pos = sf::Mouse::getPosition(window);
        if(v.second->view.getViewport().contains((double)pos.x / viz_p.res_x, (double)pos.y / viz_p.res_y)){
            window.setView(v.second->view);
            mouse_pos = window.mapPixelToCoords(sf::Vector2i(pos.x, pos.y));
            info_box.setPosition(pos.x, pos.y);
            info_box.update(mouse_pos.x, mouse_pos.y,
                            v.second->getValue(mouse_pos.x, mouse_pos.y));
            in_view = true;
            break;
        }
    }

    for(auto v: views){
        window.setView(v.second->view);
        window.draw(*v.second);
        for(auto obj:v.second->viz_objects){
            window.draw(*obj);
        }
//        if(text_num == i){
//
//        }
    }
    window.setView(standard_view);
    if(in_view) {
        window.draw(info_box);
    }
    window.setView(gui_view);
    for(auto ob:gui_obs){
        window.draw(*ob);
    }

    window.display();
}

RealtimeViz::RealtimeViz(GMRF* gmrf, const Eigen::Matrix<double, Eigen::Dynamic,
                 Eigen::Dynamic>* gt_ptr, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* f_x_gt_ptr, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* f_y_gt_ptr , Viz_Params viz_p_, GmrfParams gmrf_p):
        window(),
        viz_p(viz_p_),
        gmrf_p(gmrf_p),
        settings(24, 0, 0, 4, 5, Default, false),
        info_box(14),
        text_obj(14)
{

    int num_maps_x = viz_p.draw_gt + viz_p.draw_uncertainty + viz_p.draw_map;
    double fac = gmrf_p.size_y/gmrf_p.size_x;

    map_size_x =  (1 - viz_p.margin - viz_p.margin*(double)num_maps_x - viz_p.infobox)/(double)num_maps_x;

    if(viz_p.res_y == 0){
        // s = x * (1 - margin - margin*num_maps - infobox)/num_maps
        // y = s/0.98
        viz_p.res_y = (double)viz_p.res_x * (2.0f * map_size_x*fac + viz_p.margin * 3);

        window.create(sf::VideoMode(viz_p.res_x, viz_p.res_y), "GMRF", sf::Style::None, settings);
    }
    double margin_y = viz_p.margin * (double)viz_p.res_x/(double)viz_p.res_y;
    map_size_y = (1.0f - 3.0f*margin_y)/2.0f;
    gui_view = sf::View(sf::Vector2f(viz_p.infobox*(double)viz_p.res_x/2.0f,  round(viz_p.res_y/2.0f)), sf::Vector2f(viz_p.infobox*(double)viz_p.res_x, viz_p.res_y - viz_p.margin*(double)viz_p.res_x * 2.0f));
    standard_view = window.getView();
    window.setPosition(sf::Vector2i(viz_p.offset_x, viz_p.offset_y));

    gui_obs.push_back(&text_obj);
    int num_maps = 0;
    gui_view.setViewport(sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * num_maps_x, margin_y, viz_p.infobox, (1 - 2.0f * margin_y)));
    if(viz_p.draw_map) {
        views.emplace(map_types::field_map, new ColorMap(&gmrf->estimates_map,
                                                     gmrf_p.size_x, gmrf_p
                                                     .size_y, "plasma"));
        views[map_types::field_map]->view.setViewport(
                sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * (double)num_maps, margin_y, map_size_x, map_size_y));
        views.emplace(map_types::f_field_map, new FlowMap(gmrf->getFlowXMap(), gmrf->getFlowYMap(), gmrf_p.size_x, gmrf_p.size_y, gmrf_p.N_X, gmrf_p.N_Y));
        views[map_types::f_field_map]->view.setViewport(sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * (double)num_maps, margin_y*2 + map_size_y, map_size_x, map_size_y));
        num_maps+=1;
    }
    if(viz_p.draw_uncertainty){
        views.emplace(map_types::std_map, new ColorMap(&gmrf->stds_map,
                                                       gmrf_p.size_x, gmrf_p
                                                       .size_y, "hot"));
        views[map_types::std_map]->view.setViewport(sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * (double)num_maps, margin_y, map_size_x, map_size_y));
        views.emplace(map_types::f_std_map, new ColorMap(&gmrf->f_stds_map, gmrf_p.size_x, gmrf_p.size_y, "hot"));
        views[map_types::f_std_map]->view.setViewport(sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * (double)num_maps, margin_y*2 + map_size_y, map_size_x, map_size_y));
        num_maps+=1;
    }
    if(viz_p.draw_gt){
        views.emplace(map_types::gt_map, new ColorMap(gt_ptr, gmrf_p.size_x,
                                                      gmrf_p.size_y, "plasma"));
        views[map_types::gt_map]->view.setViewport(sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * (double)num_maps, margin_y, map_size_x, map_size_y));
        views.emplace(map_types::f_gt_map, new FlowMap(f_x_gt_ptr, f_y_gt_ptr, gmrf_p.size_x, gmrf_p.size_y, f_x_gt_ptr->rows(), f_x_gt_ptr->cols()));
        views[map_types::f_gt_map]->view.setViewport(sf::FloatRect(viz_p.margin + (map_size_x + viz_p.margin) * (double)num_maps, margin_y*2 + map_size_y, map_size_x, map_size_y));
    }
    window.setView(views.begin()->second->view);
//    window.setFramerateLimit(30); // or 20 or 30

}

void RealtimeViz::addVizObject(DrawableObject *obj, map_types map_id) {
    if(!views.contains(map_id)){
        std::cout<<"WARNING: ColorMap ID doesn't exist!\n";
        return;
    }
    obj->prepareViz();
    views[map_id]->viz_objects.push_back(obj);
}


void RealtimeViz::setTime(double t) {
    this->t = t;
    text_obj.setTime(t);
}

void RealtimeViz::update() {
    for(auto v:views){
        v.second->updateData();
    }
}

void RealtimeViz::addVizObject(sf::Drawable *obj, map_types map_id) {
    if(!views.contains(map_id)){
        std::cout<<"WARNING: ColorMap ID doesn't exist!\n";
        return;
    }
    views[map_id]->viz_objects.push_back(new PrimitiveWrapper(obj));

}

RealtimeViz::~RealtimeViz() {
    for(auto v:views){
        delete v.second;
    }
}

View::View(double sx, double sy):sx(sx), sy(sy), view(sf::Vector2f((double)sx/2, (double)sy/2), sf::Vector2f((double)sx, -(double)sy)){

}
ColorMap::ColorMap(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* data, double sx, double sy,
                   std::string cmap): View(sx, sy){
    data_map = data;
    colormap = colormaps::getColorMap(cmap);

}

void ColorMap::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    sf::Texture tex;
    tex.create(data_map->rows(), data_map->cols());
    tex.update(pixels);
    sf::Sprite sprite;
    sprite.scale(sf::Vector2f(1.0f*sx/data_map->rows(), 1.0f*sy/data_map->cols()));
    sprite.setTexture(tex);
    target.draw(sprite, states);
}

double ColorMap::getValue(double x, double y) {
    int x_ = std::floor(data_map->rows()/sx * x);
    int y_ = std::floor(data_map->cols()/sy * y);
    return data_map->operator()(x_, y_);
}

void ColorMap::updateData() {
    int counter = 0;
    double lb = data_map->minCoeff();
    double ub = data_map->maxCoeff();
    for(int j = 0; j < data_map->cols(); j++){
        for(int i = 0; i < data_map->rows(); i++){
            double val = (data_map->operator()(i, j) - lb)/(ub - lb);
            if(ub == lb){
                val = 0;
            }
            if(std::isnan(val) || val < 0 || val > 1){
                val = 0;
            }
            int id_val = floor(val*255);
            pixels[counter] = std::max(std::min(255.0*colormap[id_val][0], 255.0), 0.0);
            pixels[counter+1] = std::max(std::min(255.0*colormap[id_val][1], 255.0), 0.0);
            pixels[counter+2] = std::max(std::min(255.0*colormap[id_val][2], 255.0), 0.0);
            pixels[counter+3] = 255;
            counter += 4;
        }
    }
}

InfoBox::InfoBox(int charactersize) {
    if(!font.loadFromFile(RES_LOC "/Roboto-Regular.ttf")){
        throw std::runtime_error("Font not found. Check path!");
    }
    text.setFont(font);
    charSize = charactersize;
    text.setCharacterSize(charactersize);
    text.setPosition(charactersize, -charactersize);
    text.setFillColor(sf::Color::White);
    box.setFillColor(sf::Color(0, 0, 0, 100));
    box.setOutlineColor(sf::Color::White);
    box.setOutlineThickness(0.8);
}

void InfoBox::update(double x, double y, double v) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << "(" << x << ", " << y << ")\n"  << v;
    text.setString(stream.str());
    sf::FloatRect text_bb = text.getLocalBounds();
    text.setOrigin(0.0f, text_bb.height);
    box.setOrigin(0.0f, text_bb.height);
    box.setPosition( charSize - text_bb.height*0.2, -charSize - text_bb.height*0.1);
    box.setSize(sf::Vector2f(text_bb.width + text_bb.height*0.4, text_bb.height*1.4));
}

void InfoBox::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    states.transform *= getTransform();
    target.draw(box, states);
    target.draw(text, states);
}

void GuiText::setTime(double t) {
    this->t = t;
    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << t << "s";
    text.setString(stream.str());
}

void GuiText::draw(sf::RenderTarget &target, sf::RenderStates states) const {
//    states.transform *= getTransform();
    target.draw(text, states);
}

GuiText::GuiText(int charactersize) {

    if(!font.loadFromFile(RES_LOC "/Roboto-Regular.ttf")){
        throw std::runtime_error("Font not found. Check path!");
    }
    text.setFont(font);
    char_size = charactersize;
    text.setCharacterSize(charactersize);
    text.setPosition(charactersize, charactersize);
    text.setFillColor(sf::Color::White);
    text.setOrigin(0.0f, 0.0f);

}

FlowMap::FlowMap(const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>* ptr_flow_x, const Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>>* ptr_flow_y, double sx, double sy, int n_x, int n_y): View(sx, sy),
                                          flow_x(ptr_flow_x->data(), n_x, n_y),
                                          flow_y(ptr_flow_y->data(), n_x, n_y)
{
    for(int i = 0; i < n_x; i++){
        for(int j = 0; j < n_y; j++){
            lines.push_back(sf::RectangleShape(sf::Vector2f(sx / (double)n_x * 0.45f, sy / (double)n_y * 0.2f)));
            lines.back().setPosition(sx / ((double)n_x - 1.0) * (double)i, sy / ((double)n_y - 1.0) * (double)j);

        }
    }
    this->n_x = n_x;
    this->n_y = n_y;

}

FlowMap::FlowMap(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *ptr_flow_x,
                 const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *ptr_flow_y, double sx, double sy, int n_x,
                 int n_y): View(sx, sy),
                           flow_x(ptr_flow_x->data(), n_x, n_y),
                           flow_y(ptr_flow_y->data(), n_x, n_y)
    {
    for(int i = 0; i < n_x; i++){
        for(int j = 0; j < n_y; j++){
            lines.push_back(sf::RectangleShape(sf::Vector2f(sx / (double)n_x * 0.6f, sy / (double)n_y * 0.2f)));
            lines.back().setPosition(sx / ((double)n_x - 1.0) * (double)i, sy / ((double)n_y - 1.0) * (double)j);
        }
    }
    this->n_x = n_x;
    this->n_y = n_y;

}

void FlowMap::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for(const auto& l: lines){
        target.draw(l, states);
    }
}

void FlowMap::updateData() {
    int counter = 0;
    auto magn = (flow_x.array().pow(2) + flow_y.array().pow(2)).sqrt();
    double max_magn = magn.maxCoeff();
//    auto angle = Eigen::atan(flow_y.array() / flow_x.array());
    for(int i = 0; i < flow_x.rows(); i++){
        for(int j = 0; j < flow_x.cols(); j++){
//            lines[counter].setFillColor(sf::Color(255, 255, 255, 255 * magn(i, j) / max_magn));
            lines[counter].setFillColor(sf::Color(255, 255, 255));
            lines[counter].setRotation(std::atan2(flow_y(i, j), flow_x(i, j)) * 180.0 / M_PI);
            counter += 1;
        }
    }
}

double FlowMap::getValue(double x, double y) {
    int x_id = (int)floor((double)n_x * x / sx);
    int y_id = (int)floor((double)n_y * y / sy);
    return(sqrt(pow(flow_x(x_id, y_id), 2) + pow(flow_y(x_id, y_id), 2)));
}

sf::RectangleShape rectFromPoints(sf::Vertex v1, sf::Vertex v2, double thickness){
    double angle = atan2(v2.position.y - v1.position.y, v2.position.x - v1.position.x);
    double length = sqrt(pow(v2.position.x - v1.position.x, 2) + pow(v2.position.y - v1.position.y, 2));
    sf::RectangleShape r(sf::Vector2f(length, thickness));
    double x0 = v1.position.x;
    double y0 = v1.position.y;
    r.rotate(angle * 180/M_PI);
    r.move(x0, y0);
    return r;

}

ObjectsViz::ObjectsViz(Objects *objects):objects(objects) {
    // Create all the Rectangles
    if(!objects->objects[wallType::wall].empty()){
        for(int i = 0; i < objects->objects[wallType::wall].size(); i+=2){
            viz_objects.push_back(
                    rectFromPoints(objects->objects[wallType::wall][i],
                                    objects->objects[wallType::wall][i + 1],
                                    0.05));
        }
    }
    if(!objects->objects[wallType::inlet].empty()){
        for(int i = 0; i < objects->objects[wallType::inlet].size(); i+=2){
            viz_objects.push_back(
                    rectFromPoints(objects->objects[wallType::inlet][i],
                                    objects->objects[wallType::inlet][i + 1],
                                    0.05));
        }
    }
    if(!objects->objects[wallType::outlet].empty()){
        for(int i = 0; i < objects->objects[wallType::outlet].size(); i+=2){
            viz_objects.push_back(
                    rectFromPoints(objects->objects[wallType::outlet][i],
                                    objects->objects[wallType::outlet][i + 1],
                                    0.05));
        }
    }
    if(!objects->objects[wallType::source].empty()){
        for(int i = 0; i < objects->objects[wallType::source].size(); i+=2){
            viz_objects.push_back(
                    rectFromPoints(objects->objects[wallType::source][i],
                                    objects->objects[wallType::source][i + 1],
                                    0.05));
        }
    }
}

void ObjectsViz::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for (const auto& r: viz_objects){
        target.draw(r, states);
    }
}

PrimitiveWrapper::PrimitiveWrapper(sf::Drawable *obj) {
    this->obj = obj;
}

void PrimitiveWrapper::draw(sf::RenderTarget &target,
                            sf::RenderStates states) const {
    target.draw(*obj, states);
}
