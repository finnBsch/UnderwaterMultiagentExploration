//
// Created by finn on 10/7/23.
//
#include <fstream>
#include "visual.h"
#include "configs.h"
#include "rrt_configs.h"
#include "sim_interface.h"
#include "logger.h"


ValueMap::ValueMap(
        const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
                *values_ptr, double sx, double sy, std::string colormap_name) {
    this->values_ptr = values_ptr;
    this->sx = sx;
    this->sy = sy;
    colormap = colormaps::getColorMap(colormap_name);
    pixel_values.resize(values_ptr->size() * subsampling_fac * subsampling_fac * 4);
}

void ValueMap::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    states.transform *= getTransform();
    sf::Texture texture;
    texture.create(values_ptr->rows()*subsampling_fac, values_ptr->cols()*subsampling_fac);
    texture.update(pixel_values.data());
    sf::Sprite sprite(texture);
    sprite.scale(sf::Vector2f(sx/(values_ptr->rows()*subsampling_fac - 1),
                              sy/(values_ptr->cols()*subsampling_fac - 1)));
    target.draw(sprite, states);
}

double ValueMap::interpolate_value(double x, double y) const {
    return interp(x, y);
//    int data_nx = values_ptr->rows();
//    int data_ny = values_ptr->cols();
//    int x0 = round(x / sx * (data_nx - 1));
//    int y0 = round(y / sy * (data_ny - 1));
//    return (*values_ptr)(x0, y0);
}


void ValueMap::updateData() {
    // Create flattened vector of x values corresponding to all entries for
    // value_ptr, then flatten it
    Eigen::Matrix<double, Eigen::Dynamic, 1> x_values = Eigen::VectorXd::LinSpaced(
            values_ptr->rows(), 0, sx);
    Eigen::Matrix<double, Eigen::Dynamic, 1> y_values = Eigen::VectorXd::LinSpaced(
            values_ptr->cols(), 0, sy);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x_ = x_values.replicate(
            1, y_values.size()).transpose();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y_ = y_values.transpose().replicate(
            x_values.size(), 1).transpose();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> z_values =
            (*values_ptr).transpose();

    interp.setData( x_.size(), x_.data() , x_.size(),  y_.data(), x_.size(),  z_values.data());
    int nx = values_ptr->rows()*subsampling_fac;
    int ny = values_ptr->cols()*subsampling_fac;
    double max_value = z_values.maxCoeff();
    double min_value = z_values.minCoeff();
    // TODO Implement other interpolations
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x = sx * i / (nx - 1);
            double y = sy * j / (ny - 1);
            double value = interpolate_value(x, y);
            value = (value - min_value)/(max_value - min_value);
            // Clip value between 0 and 1
            value = std::max(0.0, std::min(1.0, value));
            if (max_value == min_value){
                value = 0;
            }
            int id_floor = floor(value*255);
            int r = round(255.0*(double)((colormap[id_floor + 1][0] - colormap[id_floor][0])*(value*255 - id_floor) + colormap[id_floor][0]));
            int g = round(255.0*(double)((colormap[id_floor + 1][1] - colormap[id_floor][1])*(value*255 - id_floor) + colormap[id_floor][1]));
            int b = round(255.0*(double)((colormap[id_floor + 1][2] - colormap[id_floor][2])*(value*255 - id_floor) + colormap[id_floor][2]));
            // Clip r, g, and b
            r = std::max(0, std::min(255, r));
            g = std::max(0, std::min(255, g));
            b = std::max(0, std::min(255, b));
            pixel_values[4 * (i + j * values_ptr->rows()*subsampling_fac)] = r;
            pixel_values[4 * (i + j * values_ptr->rows()*subsampling_fac) + 1] = g;
            pixel_values[4 * (i + j * values_ptr->rows()*subsampling_fac) + 2] = b;
            pixel_values[4 * (i + j * values_ptr->rows()*subsampling_fac) + 3] = 255;
            // TODO Add colormap
        }
    }
}

Agent::Agent():
    line_strip(0.02)
    {
    agent_shape.setPointCount(3);
    agent_shape.setPoint(0, sf::Vector2f(0, 0.2));
    agent_shape.setPoint(1, sf::Vector2f(0.0, -0.2));
    agent_shape.setPoint(2, sf::Vector2f(0.45, 0.0));
    agent_shape.setFillColor(sf::Color(120, 120, 120));
    agent_shape.setOutlineColor(sf::Color::Black);
    agent_shape.setOutlineThickness(0.06);
}

void Agent::setState(double x, double y, double phi) {
    x_coords.push_back((float) x);
    y_coords.push_back((float) y);
    line_strip.updatePoints(x_coords, y_coords);
    this->x = x;
    this->y = y;
    this->phi = phi;
}

void Agent::updateData() {
    agent_shape.setPosition((float)x, (float)y);
    agent_shape.setRotation((float) (phi * 180 / M_PI));
}

void Agent::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    states.transform *= getTransform();
    target.draw(line_strip, states);
    target.draw(agent_shape, states);
}

double Agent::get_x() const {
    return x;
}

double Agent::get_y() const {
    return y;
}

/**
 * Export one frame to base path
 */
void Visual::exportFrame() {
    sf::Texture frame;
    frame.setSmooth(true);
    frame.create(conf.renderWidth, conf.renderHeight);
    frame.update(window);
    frame.setSmooth(true);
    auto content = frame.copyToImage();
    std::string filename = conf.base_path + "/" + std::to_string(frame_counter) + ".jpg";
    content.saveToFile(filename);
    frame_counter++;
}

Visual::Visual(VizConf conf, GMRF *gmrf,
               const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *gt_ptr,
               const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *est_ptr,
               const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *source_ptr,
               const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> *variance_ptr,
               GmrfParams gmrf_params) : view(sf::Vector2f((double)
               gmrf_params.size_x/2, (double)gmrf_params.size_y/2), sf::Vector2f((double)gmrf_params.size_x, -(double)gmrf_params.size_y)) {
    this->conf = conf;
    sx = (float)gmrf_params.size_x;
    sy = (float)gmrf_params.size_y;
    sf::ContextSettings settings;
    settings.antialiasingLevel = 16;
    window.create(sf::VideoMode(conf.renderWidth, conf.renderHeight), "x", sf::Style::None, settings);
    window.setView(view);
    window.setPosition(sf::Vector2i(0, 0));
    value_maps.push_back(new ValueMap(est_ptr, gmrf_params.size_x, gmrf_params.size_y, "plasma"));
    value_maps.push_back(new ValueMap(gt_ptr, gmrf_params.size_x, gmrf_params.size_y, "plasma"));
    value_maps.push_back(new ValueMap(source_ptr, gmrf_params.size_x, gmrf_params.size_y, "viridis"));
    value_maps.push_back(new ValueMap(variance_ptr, gmrf_params.size_x, gmrf_params.size_y, "hot"));
    if (conf.draw_gmrf_grid) {
        drawables.push_back(new GmrfGrid(&gmrf_params));
    }
    if (conf.draw_connectivity){
        Connectivity* con = new Connectivity(&agents);
        addDrawable(con);
    }
}

void Visual::addDrawable(MyDrawable *drawable) {
    drawables.push_back(drawable);
}

void Visual::draw() {
    window.clear(sf::Color::Black);
    sf::Event ev;
    while (window.pollEvent(ev)) {
        if(!conf.renderVideo){
            // Catch Keys, space means pause, 1 - 4 change color map
            if (ev.type == sf::Event::KeyPressed) {
                if (ev.key.code == sf::Keyboard::Space) {
                    paused = !paused;
                }
                if (ev.key.code == sf::Keyboard::Num1) {
                    conf.colorMapMode = ColorMapMode::estimates;
                }
                if (ev.key.code == sf::Keyboard::Num2) {
                    conf.colorMapMode = ColorMapMode::ground_truth;
                }
                if (ev.key.code == sf::Keyboard::Num3) {
                    conf.colorMapMode = ColorMapMode::source_predictions;
                }
                if (ev.key.code == sf::Keyboard::Num4) {
                    conf.colorMapMode = ColorMapMode::variance;
                }
                if (ev.key.code == sf::Keyboard::Enter) {
                    conf.close_up = !conf.close_up;
                }
            }
        }
    }
    if (conf.close_up) {
        is_closed_up = true;
        float min_x = std::numeric_limits<float>::max();
        float max_x = 0;
        float min_y = std::numeric_limits<float>::max();
        float max_y = 0;
        for (auto& agent : agents){
            min_x = std::min(min_x, (float)agent->get_x());
            max_x = std::max(max_x, (float)agent->get_x());
            min_y = std::min(min_y, (float)agent->get_y());
            max_y = std::max(max_y, (float)agent->get_y());
        }
        view.setSize(sx / 2, -sy / 2);
        float xc = (min_x + max_x)/2;
        float yc = (min_y + max_y)/2;
        if (xc < sx/4){
            xc = sx/4;
        }
        if (xc > sx - sx/4){
            xc = sx - sx/4;
        }
        if (yc < sy/4){
            yc = sy/4;
        }
        if (yc > sy - sy/4){
            yc = sy - sy/4;
        }
        view.setCenter(xc, yc);
        window.setView(view);
    }
    else if(is_closed_up){
        is_closed_up = false;
        view.setSize(sx, -sy);
        view.setCenter(sx/2, sy/2);
        window.setView(view);
    }
    if (conf.colorMapMode == ColorMapMode::estimates){
        value_maps[0]->updateData();
        window.draw(*value_maps[0]);
    }
    else if (conf.colorMapMode == ColorMapMode::ground_truth){
        value_maps[1]->updateData();
        window.draw(*value_maps[1]);
    }
    else if (conf.colorMapMode == ColorMapMode::source_predictions){
        value_maps[2]->updateData();
        window.draw(*value_maps[2]);
    }
    else if (conf.colorMapMode == ColorMapMode::variance){
        value_maps[3]->updateData();
        window.draw(*value_maps[3]);
    }
    for (auto &drawable : drawables) {
        drawable->updateData();
        window.draw(*drawable);
    }
    for (auto &agent : agents) {
        agent->updateData();
        window.draw(*agent);
    }
    window.display();
    if (conf.renderVideo) {
        exportFrame();
    }
}

void Visual::addAgent(Agent *agent) {
    agents.push_back(agent);
}

CrossShape::CrossShape(double x, double y, double size, double line_thickness) {
    this->x = x;
    this->y = y;
    this->size = size;
    this->line_thickness = line_thickness;
    this->setPosition(x, y);
}

void CrossShape::draw(sf::RenderTarget &target, sf::RenderStates states) const {

    states.transform *= getTransform();
    sf::RectangleShape rect(sf::Vector2f(line_thickness, size));
    rect.setFillColor(sf::Color(255, 255, 255, 100));
    rect.setPosition(sf::Vector2f(-line_thickness/2, -size/2));
    target.draw(rect, states);
    rect.setSize(sf::Vector2f(size, line_thickness));
    rect.setPosition(sf::Vector2f(-size/2, -line_thickness/2));
    target.draw(rect, states);
}

GmrfGrid::GmrfGrid(const GmrfParams *gmrf_params) {
    this->gmrf_params = gmrf_params;
    double dx = (double)gmrf_params->size_x / ((double)gmrf_params->N_X - 1);
    double dy = (double)gmrf_params->size_y / ((double)gmrf_params->N_Y - 1);
    for (int i = 0; i < gmrf_params->N_X; i++) {
        for (int j = 0; j < gmrf_params->N_Y; j++) {
            double x = i * dx;
            double y = j * dy;
            double size = dx/3;
            double line_thickness = size / 10;
            crosses.emplace_back(x, y, size, line_thickness);
        }
    }
}

void GmrfGrid::updateData() {

}

void GmrfGrid::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    states.transform *= getTransform();
    for (auto &cross : crosses) {
        target.draw(cross, states);
    }
}

sf::RectangleShape rectFromPoints2(sf::Vertex v1, sf::Vertex v2, double thickness){
    double angle = atan2(v2.position.y - v1.position.y, v2.position.x - v1.position.x);
    double length = sqrt(pow(v2.position.x - v1.position.x, 2) + pow(v2.position.y - v1.position.y, 2));
    sf::RectangleShape r(sf::Vector2f(length, thickness));
    double x0 = v1.position.x;
    double y0 = v1.position.y;
    r.rotate(angle * 180/M_PI);
    r.move(x0, y0);
    return r;

}

ScenarioViz::ScenarioViz(Objects *objects):objects(objects){
    if(!objects->objects[wallType::wall].empty()){
        for(int i = 0; i < objects->objects[wallType::wall].size(); i+=2){
            obstacles.push_back(
                    rectFromPoints2(objects->objects[wallType::wall][i],
                                    objects->objects[wallType::wall][i + 1],
                                    0.05));
        }
    }
    if(!objects->objects[wallType::inlet].empty()){
        for(int i = 0; i < objects->objects[wallType::inlet].size(); i+=2){
            obstacles.push_back(
                    rectFromPoints2(objects->objects[wallType::inlet][i],
                                    objects->objects[wallType::inlet][i + 1],
                                    0.05));
        }
    }
    if(!objects->objects[wallType::outlet].empty()){
        for(int i = 0; i < objects->objects[wallType::outlet].size(); i+=2){
            obstacles.push_back(
                    rectFromPoints2(objects->objects[wallType::outlet][i],
                                    objects->objects[wallType::outlet][i + 1],
                                    0.05));
        }
    }
    if(!objects->objects[wallType::source].empty()){
        for(int i = 0; i < objects->objects[wallType::source].size(); i+=2){
            obstacles.push_back(
                    rectFromPoints2(objects->objects[wallType::source][i],
                                    objects->objects[wallType::source][i + 1],
                                    0.05));
        }
    }
}

void
ScenarioViz::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for (const auto& r: obstacles){
        target.draw(r, states);
    }
}


void renderVideo(std::string path){
    // Load the scenario.
    std::ifstream config_file;
    config_file.open(path + "/config.txt");
    std::string scenario_name;
    config_file >> scenario_name;
    ConfigFileScenario scenario(scenario_name);

    // Load the module params
    SQPParams sqp_p(path + "/configs/sqp_conf.toml");
    GmrfParams gmrf_p(path + "/configs/gmrf_conf.toml");
    RRTParams rrt_p(path + "/configs/rrt_conf.toml");

    // Adjust params
    sqp_p.size_x = scenario.get_sx();
    sqp_p.size_y = scenario.get_sy();
    gmrf_p.size_x = scenario.get_sx();
    gmrf_p.size_y = scenario.get_sy();
    rrt_p.size_x = (float)scenario.get_sx();
    rrt_p.size_y = (float)scenario.get_sy();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> estimates;
    estimates.resize(gmrf_p.N_X, gmrf_p.N_Y);
    estimates.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> source_predictions;
    source_predictions.resize(gmrf_p.N_X, gmrf_p.N_Y);
    source_predictions.setZero();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> variances;
    variances.resize(gmrf_p.N_X, gmrf_p.N_Y);
    variances.setZero();

    // Load the modules
    GMRF gmrf(gmrf_p, &scenario);

    // Viz
    ScenarioViz scen_viz(scenario.getObjects());
    VizConf viz_conf;
    if (!viz_conf.renderVideo){
        viz_conf.renderHeight = sf::VideoMode::getDesktopMode().height;
        viz_conf.renderWidth = sf::VideoMode::getDesktopMode().width;
    }
    double aspect_ratio = (double)viz_conf.renderHeight / (double)viz_conf.renderWidth;
    double scenario_ar = scenario.get_sy() / scenario.get_sx();
    double ratio = scenario_ar / aspect_ratio;
    if(ratio > 1){
        viz_conf.renderWidth = viz_conf.renderWidth / ratio;
        std::cout << "Ratio: " << ratio << std::endl;
        std::cout << "Render width: " << viz_conf.renderWidth << std::endl;
    } else {
        viz_conf.renderHeight = viz_conf.renderHeight * ratio;
        std::cout << "Ratio: " << ratio << std::endl;
        std::cout << "Render height: " << viz_conf.renderHeight << std::endl;
    }
    Visual viz(viz_conf, &gmrf, scenario.getField(), &estimates, &source_predictions, &variances, gmrf_p);
    viz.addDrawable(&scen_viz);
    std::vector<Agent*> agents;
    for(int i = 0; i < sqp_p.num_agents; i++){
        agents.push_back(new Agent());
        viz.addAgent(agents[i]);
    }
    // Load the states
    std::ifstream state_file;
    std::ifstream gmrf_file;
    EigenVectorReader state_reader(path + "/state.csv");
    EigenVectorReader gmrf_reader(path + "/gmrf_estimates.csv");
    EigenMatrixReader field_cov_reader(path + "/field_covariances.csv");
    EigenMatrixReader flowx_cov_reader(path + "/flowx_covariances.csv");
    EigenMatrixReader flowy_cov_reader(path + "/flowy_covariances.csv");
    EigenMatrixReader source_pred_reader(path + "/source_estimates.csv");

    int num_nodes = gmrf_p.N_X * gmrf_p.N_Y;

    if(state_reader.getRows() != sqp_p.num_states*sqp_p.num_agents){
        std::cout << "WARNING: Number of states in "
                     "replay "
                     "file does "
                     "not "
                     "match number of states in config file" <<
                  std::endl;

    }
    int num_states = state_reader.getRows()/sqp_p.num_agents;
    double t = 0;
    double dt = sqp_p.dt / viz_conf.speed;
    if (viz_conf.renderVideo){
        dt = 0;
    }
    sf::Clock clock;
    clock.restart();
    // Start reading the actual data.
    while (true) {
        bool success = true;
        if(!viz.paused && clock.getElapsedTime().asSeconds() > dt) {
            clock.restart();
            auto state = state_reader.readVector(success);
            if (!success) {
                break;
            }
            t = state_reader.getTime();

            source_predictions = source_pred_reader.readMatrix(success);
            variances = field_cov_reader.readMatrix(success);


            // Slice num_nodes out of gmrf_row and reshape to N_X, N_Y
            auto gmrf_row = gmrf_reader.readVector(success);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> sliced = gmrf_row
                    .head(num_nodes);
            estimates = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(
                    sliced.data(), gmrf_p.N_X, gmrf_p.N_Y);
            // Clip estimates between 0 and 1
            estimates = estimates.cwiseMax(0);
            estimates = estimates.cwiseMin(1);

            for (int i = 0; i < sqp_p.num_agents; i++) {
                agents[i]->setState(state(i * num_states, 0),
                                    state(i * num_states +
                                          1, 0), state(i * num_states + 2, 0));
            }
            double xs;
            double ys;
            bool succ;
        }
        viz.draw();
    }
}

Connectivity::Connectivity(std::vector<Agent*>* agents) {
    this->agents = agents;
    // Create one line for each agent pair, this results
    for(int i = 0; i < agents->size(); i++){
        for(int j = i+1; j < agents->size(); j++){
            connectivity_lines.emplace_back(0.01, 0.2, 0.2);
            connected.push_back(false);
        }
    }
}

void Connectivity::updateData() {
    if(agents->size() * (agents->size() - 1)/2 != connectivity_lines.size()){
        connectivity_lines.clear();
        connected.clear();
        std::cout << "Oopsie" << std::endl;
        for (int i = 0; i < agents->size(); i++) {
            for (int j = i + 1; j < agents->size(); j++) {
                connectivity_lines.emplace_back(0.01, 0.2, 0.2);
                connected.push_back(false);
            }
        }
    }
    int k = 0;
    for(int i = 0; i < agents->size(); i++){
        for(int j = i+1; j < agents->size(); j++){
            double x0 = agents->at(i)->get_x();
            double y0 = agents->at(i)->get_y();
            double x1 = agents->at(j)->get_x();
            double y1 = agents->at(j)->get_y();
            if (sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2)) <= max_dist){
                connected[k] = true;
            } else {
                connected[k] = false;
            }
            std::vector<float> x_coords = {(float)x0, (float)x1};
            std::vector<float> y_coords = {(float)y0, (float)y1};

            connectivity_lines[k].setPoints(x_coords, y_coords);
            k++;
        }
    }
}

void
Connectivity::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    for(int i = 0; i < connectivity_lines.size(); i++){
        if(connected[i]){
            target.draw(connectivity_lines[i], states);
        }
    }
}
