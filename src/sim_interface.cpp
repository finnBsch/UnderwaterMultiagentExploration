//
// Created by finn on 7/2/23.
//
#include <sys/stat.h>

#include "sim_interface.h"

/**
 * Runs the simulation until stopping criteria is met.
 * @return bool indicating successful simulation.
 */
bool SimInterface::runFull() {
    while(!done()){
        bool successful_step = step();
        if(!successful_step){
            postProcess();
            return false;
        }
    }
    postProcess();
    return true;
}

FullInterface::FullInterface(std::string log_base_path, SQPParams *sqp_params, RRTParams *rrt_params,
                             GmrfParams *gmrf_params, ConfigFileScenario* scenario,
                             StoppingCritereon stop_critereon,
                             bool visualize, bool log, double x0, double y0):
                             sqp_params(sqp_params), rrt_params(rrt_params),
                             gmrf_params(gmrf_params) {
    if(log) {
        struct stat sb;
        if (stat(log_base_path.c_str(), &sb) == 0) {
            std::cout << "Log directory already exists. Exiting." << std::endl;
            exit(1);
        } else {
            mkdir(log_base_path.c_str(), 0777);
            std::cout << "Created Experiment Directory in " << log_base_path
                      << std::endl;
        }
        if (!stat((log_base_path + "/configs").c_str(), &sb) == 0) {
            mkdir((log_base_path + "/configs").c_str(), 0777);
        }
    }
    this->visualize = visualize;
    this->log = log;
    this->scenario = scenario;
    this->log_base_path = log_base_path;

    // Resize the runtime members
    current_state.resize(sqp_params->num_states * sqp_params->num_agents, 1);
    current_state.setZero();
    if(sqp_params->num_agents % 2 == 1){
        for(int i = 0; i < sqp_params->num_agents; i++){
            current_state(i * sqp_params->num_states + 0) = x0 + (i - (sqp_params->num_agents - 1)/2) * 0.7;
            current_state(i * sqp_params->num_states + 1) = y0;
//            std::cout << "Agent " << i << " at " << current_state(i * sqp_params->num_states + 0) << ", " << current_state(i * sqp_params->num_states + 1) << std::endl;
        }
    }
    else {
        for(int i = 0; i < sqp_params->num_agents; i++){
            current_state(i * sqp_params->num_states + 0) = x0 + (i - (sqp_params->num_agents)/2) * 0.5;
            current_state(i * sqp_params->num_states + 1) = y0;
        }
    }
    trajectory_state.resize(sqp_params->num_states * sqp_params->num_agents,
                            1000 + 1);
    trajectory_state.setZero();
    trajectory_state.col(0) = current_state;

    // Parse the correct map size
    this->gmrf_params->size_x = scenario->get_sx();
    this->gmrf_params->size_y = scenario->get_sy();

    this->sqp_params->size_x = scenario->get_sx();
    this->sqp_params->size_y = scenario->get_sy();

    this->rrt_params->size_x = (float)scenario->get_sx();
    this->rrt_params->size_y = (float)scenario->get_sy();

    // Initialize the obstacles
    auto objects = scenario->getObjects()->objects;
    auto walls = objects[wallType::wall];
    auto inlets = objects[wallType::inlet];
    auto outlets = objects[wallType::outlet];
    auto sources = objects[wallType::source];

    for(int i = 0; i < walls.size(); i+=2){
        rrt_obstacles.push_back(new RRTObstacle(walls[i].position.x, walls[i]
                .position
                .y, walls[i + 1].position.x, walls[i + 1].position.y));
    }
    for(int i = 0; i < inlets.size(); i+=2){
        rrt_obstacles.push_back(new RRTObstacle(inlets[i].position.x, inlets[i]
                .position
                .y, inlets[i + 1].position.x, inlets[i + 1].position.y));
    }
    for(int i = 0; i < outlets.size(); i+=2) {
        rrt_obstacles.push_back(
                new RRTObstacle(outlets[i].position.x, outlets[i]
                                        .position
                                        .y, outlets[i + 1].position.x,
                                outlets[i + 1].position.y));
    }
    for(int i = 0; i < sources.size(); i+=2) {
        rrt_obstacles.push_back(
                new RRTObstacle(sources[i].position.x, sources[i]
                                        .position
                                        .y, sources[i + 1].position.x,
                                sources[i + 1].position.y));
    }
    gmrf = new GMRF(*this->gmrf_params, this->scenario);
    obs_interface = new ObstacleInterface(&rrt_obstacles);
    field_interface = new FieldInterface(gmrf);
    rrt = new RRT(this->rrt_params, &rrt_obstacles, field_interface, x0, y0);
    rrt->buildTree();
    while(!rrt->rewireDone()){}
    goal_node = rrt->getBestPath();
    rrt->buildPathWithCorridor(goal_node);

    rrt_spline_interface = new RRTSplinePathInterface(rrt,
                                                      obs_interface->getObstacles());
    sqp = new SQPMultiAgent(*sqp_params, rrt_spline_interface, field_interface,
                            current_state, obs_interface->getObstacles());
    sqp->solve(current_state, true);
    if(visualize){
        ObjectsViz* obj_viz = new ObjectsViz(scenario->getObjects());
        Viz_Params viz_params;
        viz_params.offset_y = 100;
        viz = new RealtimeViz(gmrf, scenario->getField(), scenario->getFlowFieldX(), scenario->getFlowFieldY(), viz_params, *gmrf_params);
        viz->addVizObject(obj_viz,  map_types::gt_map);
        viz->addVizObject(obj_viz,  map_types::field_map);
        gmrf->addToViz(*viz, 0, 4, 0, 5);
        viz->addVizObject(rrt, map_types::std_map);
        viz->addVizObject(sqp, map_types::std_map);
        viz->addVizObject(rrt_spline_interface, map_types::std_map);
        for(int i = 0; i < sqp_params->num_agents; i++){
            agent_positions.push_back(new sf::ConvexShape(3));
            agent_positions[i]->setPoint(0, sf::Vector2f(0, 0.3));
            agent_positions[i]->setPoint(1, sf::Vector2f(0.0, -0.3));
            agent_positions[i]->setPoint(2, sf::Vector2f(0.6, 0.0));
            agent_positions[i]->setFillColor(*sqp_params->colors[i]);
            agent_positions[i]->setOutlineColor(sf::Color::Black);
            agent_positions[i]->setOutlineThickness(0.07);
            viz->addVizObject(agent_positions[i], map_types::std_map);
        }
    }
    // Enable automatic tree updates
    rrt->updateTree();
    rrt->offsetRoot(rrt_params->offset_replan);
    replan_angle = rrt_spline_interface->getAngle(rrt_params->offset_replan);
    timer.addModuleToLog("gmrf");
    timer.addModuleToLog("mpc");
    timer.addModuleToLog("source_localization");
    measurement_locations.resize(0, 2);
    // Prepare Logger
    if(this->log) {
        // open base config file to save scenario name
        std::ofstream config_file;
        config_file.open(log_base_path + "/config.txt");
        config_file << scenario->getName();
        config_file.close();
        sqp_params->exportConfig(log_base_path);
        gmrf_params->exportConfig(log_base_path);
        rrt_params->exportConfig(log_base_path);
        logger = new SimLogger(log_base_path);
        logger->addLoggable(new LoggableEigenVector("state", &current_state));
        logger->addLoggable(&timer);
        logger->addLoggable(new LoggableEigenVector("gmrf_estimates",
                                                    &gmrf->getAllEstimates()));
        logger->addLoggable(new LoggableEigen("source_estimates",
                                              &gmrf->getSourceLikelihoods()));
        logger->addLoggable(new LoggableEigen("field_covariances",
                            &gmrf->getCovariances()));

        logger->addLoggable(new LoggableDynamicEigen("measurement_locations",
                                              &measurement_locations));
        logger->addLoggable(new LoggableEigen("flowx_covariances",
                                              &gmrf->getCovariances
                                              (measurement_type::flow_x)));
        logger->addLoggable(new LoggableEigen("flowy_covariances",
                                              &gmrf->getCovariances
                                              (measurement_type::flow_x)));

        logger->addMetric(new RMSEMetric(gmrf, scenario, rrt_obstacles, 0.2));
        logger->addMetric(new FieldLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new FieldLogLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new SourceDistMetric(gmrf, scenario));
        logger->addLoggable(new LoggablePoint("Source", &gmrf->source_x, &gmrf->source_y));
    }
}

bool FullInterface::done() {
    return current_step >= 500;
}

bool FullInterface::step() {
    if(log) {
        logger->log(t);
    }
    // RRT Part
    double theta = sqp->getTheta();
    auto thetas = rrt_spline_interface->getThetas();
    auto root_id = rrt_spline_interface->getRootId();
    if((*thetas)[root_id + rrt_params->offset_replan] < theta) {
        int keep_id = rrt_params->offset_replan;
        while(!rrt->rewireDone()){}
            goal_node = rrt->getBestPath();
        rrt->buildPathWithCorridor(goal_node);
        rrt_spline_interface->update(keep_id, 0);
        replan_angle = rrt_spline_interface->getAngle(keep_id);
        rrt->offsetRoot(rrt_params->offset_replan);
    }
    else{
        if(rrt->rewireDone()){
            // TODO Wrong angle?
            rrt->updateFromRoot(replan_angle);
        }
    }
    // Take measurements
    measurement_locations.conservativeResize(measurement_locations.rows() + sqp_params->num_agents, 2);
    for(int i = 0; i < sqp_params->num_agents; i++) {
        double x = current_state(i * sqp_params->num_states);
        double y = current_state(i * sqp_params->num_states + 1);
        if(x < 0 ){
            x = 0.0001;
        }
        if(x > scenario->get_sx()){
            x = scenario->get_sx() - 0.0001;
        }
        if(y < 0 ){
            y = 0.0001;
        }
        if(y > scenario->get_sy()){
            y = scenario->get_sy() - 0.0001;
        }
        measurement_locations(measurement_locations.rows() - sqp_params->num_agents + i, 0) = x;
        measurement_locations(measurement_locations.rows() - sqp_params->num_agents + i, 1) = y;
        gmrf->addMeasurement(x, y, scenario->getFieldValue(x, y, 0.05), t,
                            field);
        double flow_x_meas;
        double flow_y_meas;
        scenario->getFlowValue(x, y, flow_x_meas, flow_y_meas, 0.02);
        gmrf->addMeasurement(x, y, flow_x_meas, t, flow_x);
        gmrf->addMeasurement(x, y, flow_y_meas, t, flow_y);
    }
    if(gmrf->near_source){
//        std::cout << "FOUND SOURCE!!!!" << std::endl;
//        return false;
    }
    // Update estimates
    if(current_step % 6 == 0) {
        timer.startTimer("gmrf");
        gmrf->updateEstimates(t);
        timer.stopTimer("gmrf");
        timer.startTimer("source_localization");
        gmrf->sourceHypothesisNeighbor(source_x, source_y,
                                       found_hypothesis);

        timer.stopTimer("source_localization");

        field_interface->update();
    }
    // Solve the optimization problem
    timer.startTimer("mpc");
    sqp->solve(current_state, false);
    auto solution = sqp->getSol();
    timer.stopTimer("mpc");
    if (current_step % 100 == 0) {
        std::cout << "Finished step " << current_step << std::endl;
    }
    // Update the new state
    current_state(Eigen::seq(0, sqp_params->num_states * sqp_params->num_agents - 1))
        = solution(Eigen::seq(0, sqp_params->num_states * sqp_params->num_agents - 1), 1);
    t += sqp_params->dt;
    current_step++;
    trajectory_state.col(current_step) = current_state;
    if(visualize){
        for(int i = 0; i < sqp_params->num_agents; i++){
            agent_positions[i]->setPosition(solution
            (i*sqp_params->num_states, 0), solution(i*sqp_params->num_states + 1, 0));
            agent_positions[i]->setRotation(solution(i * sqp_params->num_states
            + 2, 0)
            * 180 / M_PI);
        }
        draw();
    }
    return true;
}

void FullInterface::postProcess() {
//    std::cout << "[";
    for(int i = 0; i < trajectory_state.rows(); i++){
//        std::cout << "[";
//        std::cout << trajectory_state(i, 0);
        for(int j = 1; j < trajectory_state.cols(); j++){
//            std::cout << ", " << trajectory_state(i, j);
        }
        if(i < trajectory_state.rows() - 1) {
//            std::cout << "],\n";
        }

    }
//    std::cout << "]]\n";
    delete rrt;
//    delete gmrf;
//    delete field_interface;
//    delete rrt_spline_interface;
//    delete sqp;
//    delete obs_interface;
//    for(auto obs : rrt_obstacles){
//        delete obs;
//    }
//    if(visualize){
//        delete viz;
//    }

}

void FullInterface::draw() {
    viz->draw();
}

/**
 * @brief Construct a replay interface
 * @param replay_path
 */
FullInterface::FullInterface(std::string replay_path)
{
    visualize = true;

    // Load the scenario.
    std::ifstream config_file;
    config_file.open(replay_path + "/config.txt");
    std::string scenario_name;
    config_file >> scenario_name;
    scenario = new ConfigFileScenario(scenario_name);

    // Load the module params
    sqp_params = new SQPParams(replay_path + "/configs/sqp_conf.toml");
    gmrf_params = new GmrfParams(replay_path + "/configs/gmrf_conf.toml");
    rrt_params = new RRTParams(replay_path + "/configs/rrt_conf.toml");

    // Load the modules
    gmrf = new GMRF(*this->gmrf_params, this->scenario);

    // Viz
    ObjectsViz* obj_viz = new ObjectsViz(scenario->getObjects());
    Viz_Params viz_params;
    viz_params.offset_y = 100;
    viz = new RealtimeViz(gmrf, scenario->getField(), scenario->getFlowFieldX(), scenario->getFlowFieldY(), viz_params, *gmrf_params);
    viz->addVizObject(obj_viz,  map_types::gt_map);
    viz->addVizObject(obj_viz,  map_types::field_map);
    viz->addVizObject(obj_viz,  map_types::std_map);
    gmrf->addToViz(*viz, 0, 4, 0, 5);
//    viz->addVizObject(rrt, map_types::std_map);
//    viz->addVizObject(sqp, map_types::std_map);
//    viz->addVizObject(rrt_spline_interface, map_types::std_map);
    for(int i = 0; i < sqp_params->num_agents; i++){
        agent_positions.push_back(new sf::ConvexShape(3));
        agent_positions[i]->setPoint(0, sf::Vector2f(0, 0.3));
        agent_positions[i]->setPoint(1, sf::Vector2f(0.0, -0.3));
        agent_positions[i]->setPoint(2, sf::Vector2f(0.9, 0.0));
        agent_positions[i]->setFillColor(*sqp_params->colors[i]);
        agent_positions[i]->setOutlineColor(sf::Color::Black);
        agent_positions[i]->setOutlineThickness(0.07);
        viz->addVizObject(agent_positions[i], map_types::field_map);
    }
    // Load the states
    std::ifstream state_file;
    std::ifstream gmrf_file;
    EigenVectorReader state_reader(replay_path + "/state.csv");
    EigenVectorReader gmrf_reader(replay_path + "/gmrf_estimates.csv");
    EigenMatrixReader field_cov_reader(replay_path + "/field_covariances.csv");
    EigenMatrixReader flowx_cov_reader(replay_path + "/flowx_covariances.csv");
    EigenMatrixReader flowy_cov_reader(replay_path + "/flowy_covariances.csv");




    if(state_reader.getRows() != sqp_params->num_states*sqp_params->num_agents){
        std::cout << "WARNING: Number of states in "
                                                    "replay "
                                               "file does "
                                      "not "
                                "match number of states in config file" <<
                                std::endl;

    }
    int num_states = state_reader.getRows()/sqp_params->num_agents;
    // Start reading the actual data.
    while (true) {
        bool success = true;
        auto state = state_reader.readVector(success);
        if(!success){
            break;
        }
        t = state_reader.getTime();
        gmrf->setEstimates(gmrf_reader.readVector(success), t);
        gmrf->setCovariance(field_cov_reader.readMatrix(success), t);
        for(int i = 0; i < sqp_params->num_agents; i++){
            agent_positions[i]->setPosition(state
                                                    (i*num_states, 0), state
                                                    (i*num_states + 1, 0));
            agent_positions[i]->setRotation(state(i * num_states
                                                     + 2, 0)
                                            * 180 / M_PI);
        }
        double xs;
        double ys;
        bool succ;
//        gmrf->sourceHypothesisNeighbor( xs, ys, succ);
        draw();
        sf::sleep(sf::milliseconds(100));
    }

}

FullInterface::~FullInterface() {
    delete gmrf;
    for (auto & obs : rrt_obstacles) {
        delete obs;
    }
    delete field_interface;

    delete obs_interface;
    delete rrt_spline_interface;
    delete sqp;
    if(this->visualize){
        delete viz;
    }
    if(log) {
        delete logger;
        // Write exit status into file
        std::ofstream exit_file;
        exit_file.open(log_base_path + "/exit_status.txt");
        exit_file << "success!";
        exit_file.close();
    }
    std::cout << "Exiting" << std::endl;
}

EigenVectorReader::EigenVectorReader(std::string path) {
    file.open(path);
    std::string line;
    // Read header
    std::getline(file, line);
    if(line != "EigenVector"){
        throw std::runtime_error("File is not an EigenVector file");
    }
    // Read dimensions
    std::getline(file, line);
    std::stringstream ss(line);
    std::string substr;
    getline(ss, substr, ',' );
    int rows = std::stoi(substr);
    getline(ss, substr, ',' );
    int cols = std::stoi(substr);
    data.resize(rows);
    data.setZero();
}

Eigen::VectorXd& EigenVectorReader::readVector(bool & success) {
    std::string line;
    std::getline(file, line);
    if(line.empty()){
        success = false;
        return data;
    }
    std::stringstream ss(line);
    std::string substr;
    getline(ss, substr, ',' );
    t = std::stod(substr);
    for(int i = 0; i < data.rows(); i++){
        std::getline(ss, substr, ',');
        data(i) = std::stod(substr);
    }
    return data;
}

int EigenVectorReader::getRows() {
    return data.rows();
}

double EigenVectorReader::getTime() {
    return t;
}

bool EigenVectorReader::eof() {
    return file.eof() || file.peek() == '\n';
}

EigenMatrixReader::EigenMatrixReader(std::string path) {
    file.open(path);
    std::string line;
    // Read header
    std::getline(file, line);
    if(line != "Eigen"){
        throw std::runtime_error("File is not an EigenMatrix file");
    }
    // Read dimensions
    std::getline(file, line);
    std::stringstream ss(line);
    std::string substr;
    getline(ss, substr, ',' );
    int rows = std::stoi(substr);
    getline(ss, substr, ',' );
    int cols = std::stoi(substr);
    data.resize(rows, cols);
    data.setZero();

}

Eigen::MatrixXd &EigenMatrixReader::readMatrix(bool & success) {
    std::string line;
    std::getline(file, line);
    if(line.empty()){
        success = false;
        return data;
    }
    std::stringstream ss(line);
    std::string substr;
    getline(ss, substr, ',' );
    double t = std::stod(substr);
    // Column major, so iterate over columns first
    for(int j = 0; j < data.cols(); j++){
        for(int i = 0; i < data.rows(); i++){
            std::getline(ss, substr, ',');
            data(i, j) = std::stod(substr);
        }
    }
    return data;
}

BenchmarkInterface::BenchmarkInterface(BenchmarkPlanner planner_type, std::string log_base_path,
                                       GmrfParams *gmrf_params,
                                       SQPParams *sqp_params,
                                       ConfigFileScenario *scenario,
                                       bool visualize, bool log, double x0,
                                       double y0) {
    this->planner_type = planner_type;
    this->log_base_path = log_base_path;
    if(log) {
        struct stat sb;
        if (stat(log_base_path.c_str(), &sb) == 0) {
            std::cout << "Log directory already exists. Exiting." << std::endl;
            exit(1);
        } else {
            mkdir(log_base_path.c_str(), 0777);
            std::cout << "Created Experiment Directory in " << log_base_path
                      << std::endl;
        }
        if (!stat((log_base_path + "/configs").c_str(), &sb) == 0) {
            mkdir((log_base_path + "/configs").c_str(), 0777);
        }
    }
    this->visualize = visualize;
    this->log = log;
    this->scenario = scenario;

    this->gmrf_params = gmrf_params;
    this->sqp_params = sqp_params;
    num_agents = sqp_params->num_agents;

    current_state.resize(3 * num_agents, 1);
    current_state.setZero();
    if(sqp_params->num_agents % 2 == 1){
        for(int i = 0; i < sqp_params->num_agents; i++){
            current_state(i * 3 + 0) = x0 + (i - (sqp_params->num_agents - 1)/2) * 0.7;
            current_state(i * 3 + 1) = y0;
//            std::cout << "Agent " << i << " at " << current_state(i * 3 + 0) << ", " << current_state(i * 3 + 1) << std::endl;
        }
    }
    else {
        for(int i = 0; i < sqp_params->num_agents; i++){
            current_state(i * 3 + 0) = x0 + (i - (sqp_params->num_agents)/2) * 0.5;
            current_state(i * 3 + 1) = y0;
        }
    }
    // Parse the correct map size
    this->gmrf_params->size_x = scenario->get_sx();
    this->gmrf_params->size_y = scenario->get_sy();

    this->sqp_params->size_x = scenario->get_sx();
    this->sqp_params->size_y = scenario->get_sy();
    // Initialize the obstacles
    auto objects = scenario->getObjects()->objects;
    auto walls = objects[wallType::wall];
    auto inlets = objects[wallType::inlet];
    auto outlets = objects[wallType::outlet];
    auto sources = objects[wallType::source];

    for(int i = 0; i < walls.size(); i+=2){
        obstacles.push_back(new RandomWalkObstacle({walls[i].position.x, walls[i]
                .position
                .y, walls[i + 1].position.x, walls[i + 1].position.y}));
    }
    for(int i = 0; i < inlets.size(); i+=2){
        obstacles.push_back(new RandomWalkObstacle({inlets[i].position.x, inlets[i]
                .position
                .y, inlets[i + 1].position.x, inlets[i + 1].position.y}));
    }
    for(int i = 0; i < outlets.size(); i+=2) {
        obstacles.push_back(
                new RandomWalkObstacle({outlets[i].position.x, outlets[i]
                                        .position
                                        .y, outlets[i + 1].position.x,
                                outlets[i + 1].position.y}));
    }
    for(int i = 0; i < sources.size(); i+=2) {
        obstacles.push_back(
                new RandomWalkObstacle({sources[i].position.x, sources[i]
                                        .position
                                        .y, sources[i + 1].position.x,
                                sources[i + 1].position.y}));
    }
    for(int i = 0; i < walls.size(); i+=2){
        rrt_obstacles.push_back(new RRTObstacle(walls[i].position.x, walls[i]
                .position
                .y, walls[i + 1].position.x, walls[i + 1].position.y));
    }
    for(int i = 0; i < inlets.size(); i+=2){
        rrt_obstacles.push_back(new RRTObstacle(inlets[i].position.x, inlets[i]
                .position
                .y, inlets[i + 1].position.x, inlets[i + 1].position.y));
    }
    for(int i = 0; i < outlets.size(); i+=2) {
        rrt_obstacles.push_back(
                new RRTObstacle(outlets[i].position.x, outlets[i]
                                        .position
                                        .y, outlets[i + 1].position.x,
                                outlets[i + 1].position.y));
    }
    for(int i = 0; i < sources.size(); i+=2) {
        rrt_obstacles.push_back(
                new RRTObstacle(sources[i].position.x, sources[i]
                                        .position
                                        .y, sources[i + 1].position.x,
                                sources[i + 1].position.y));
    }
    gmrf = new GMRF(*this->gmrf_params, this->scenario);
    if(planner_type == BenchmarkPlanner::Ballistic){
        random_walk = new BallisticMotion(current_state, sqp_params->dt, obstacles,
                                             sqp_params->num_agents);
    }
    else if(planner_type == BenchmarkPlanner::Brownian){
        random_walk = new BrownianRandomWalk(current_state, sqp_params->dt, obstacles,
                                             sqp_params->num_agents);
    }

    if(visualize){
        ObjectsViz* obj_viz = new ObjectsViz(scenario->getObjects());
        Viz_Params viz_params;
        viz_params.offset_y = 100;
        viz = new RealtimeViz(gmrf, scenario->getField(), scenario->getFlowFieldX(), scenario->getFlowFieldY(), viz_params, *gmrf_params);
        viz->addVizObject(obj_viz,  map_types::gt_map);
        viz->addVizObject(obj_viz,  map_types::field_map);
        gmrf->addToViz(*viz, 0, 4, 0, 5);
        for(int i = 0; i < sqp_params->num_agents; i++){
            agent_positions.push_back(new sf::ConvexShape(3));
            agent_positions[i]->setPoint(0, sf::Vector2f(0, 0.3));
            agent_positions[i]->setPoint(1, sf::Vector2f(0.0, -0.3));
            agent_positions[i]->setPoint(2, sf::Vector2f(0.6, 0.0));
            agent_positions[i]->setFillColor(*sqp_params->colors[i]);
            agent_positions[i]->setOutlineColor(sf::Color::Black);
            agent_positions[i]->setOutlineThickness(0.07);
            viz->addVizObject(agent_positions[i], map_types::std_map);
        }
    }
    // Prepare Logger
    if(this->log) {
        // open base config file to save scenario name
        std::ofstream config_file;
        config_file.open(log_base_path + "/config.txt");
        config_file << scenario->getName();
        config_file.close();
        sqp_params->exportConfig(log_base_path);
        gmrf_params->exportConfig(log_base_path);
        logger = new SimLogger(log_base_path);
        logger->addLoggable(new LoggableEigenVector("state", &current_state));
        logger->addLoggable(new LoggableEigenVector("gmrf_estimates",
                                                    &gmrf->getAllEstimates()));
        logger->addLoggable(new LoggableEigen("source_estimates",
                                              &gmrf->getSourceLikelihoods()));
        logger->addLoggable(new LoggableEigen("field_covariances",
                                              &gmrf->getCovariances()));
        logger->addLoggable(new LoggableEigen("flowx_covariances",
                                              &gmrf->getCovariances
                                                      (measurement_type::flow_x)));
        logger->addLoggable(new LoggableEigen("flowy_covariances",
                                              &gmrf->getCovariances
                                                      (measurement_type::flow_x)));
        logger->addMetric(new RMSEMetric(gmrf, scenario, rrt_obstacles, 0.2));
        logger->addMetric(new FieldLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new FieldLogLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new SourceDistMetric(gmrf, scenario));
        logger->addLoggable(new LoggablePoint("Source", &gmrf->source_x, &gmrf->source_y));



    }
}

bool BenchmarkInterface::done() {
    return current_step >= 1000;
}

bool BenchmarkInterface::step() {
    if(log) {
        logger->log(t);
    }

    if(gmrf->near_source){
//        std::cout << "FOUND SOURCE!!!!" << std::endl;
//        return false;
    }
    t += sqp_params->dt;
    current_step++;
    for(int i = 0; i < sqp_params->num_agents; i++) {
        double x = current_state(i * 3);
        double y = current_state(i * 3 + 1);
        if(x < 0 ){
            x = 0.0001;
        }
        if(x > scenario->get_sx()){
            x = scenario->get_sx() - 0.0001;
        }
        if(y < 0 ){
            y = 0.0001;
        }
        if(y > scenario->get_sy()){
            y = scenario->get_sy() - 0.0001;
        }
        gmrf->addMeasurement(x, y, scenario->getFieldValue(x, y), t,
                             field);
        double flow_x_meas;
        double flow_y_meas;
        scenario->getFlowValue(x, y, flow_x_meas, flow_y_meas);
        gmrf->addMeasurement(x, y, flow_x_meas, t, flow_x);
        gmrf->addMeasurement(x, y, flow_y_meas, t, flow_y);
    }
    if(current_step % 6 == 0) {
        gmrf->updateEstimates(t);
        gmrf->sourceHypothesisNeighbor(source_x, source_y,
                                       found_hypothesis);
    }
    random_walk->step();
    current_state = random_walk->getStates();
    if(visualize){
        for(int i = 0; i < sqp_params->num_agents; i++){
            agent_positions[i]->setPosition(current_state
                                                    (i*3, 0), current_state(i*3 + 1, 0));
            agent_positions[i]->setRotation(current_state(i * 3
                                                     + 2, 0)
                                            * 180 / M_PI);
        }
        draw();
    }
    return true;

}

void BenchmarkInterface::postProcess() {

}

void BenchmarkInterface::draw() {
    viz->draw();
}

BenchmarkInterface::~BenchmarkInterface() {
    delete gmrf;
    for (auto& obs : obstacles) {
        delete obs;
    }
    if(this->visualize){
        delete viz;
    }
    delete random_walk;
    if(log) {
        delete logger;
        // Write exit status into file
        std::ofstream exit_file;
        exit_file.open(log_base_path + "/exit_status.txt");
        exit_file << "success!";
        exit_file.close();
    }
    std::cout << "Exiting" << std::endl;
}

GMRFInterface::GMRFInterface(std::string log_base_path, GmrfParams* gmrf_params,
                             ConfigFileScenario* scenario,
                             MeasurementStrategy meas_strat,
                             bool visualize, bool log,
                             int max_meas, int batch_size):
                gmrf_params(gmrf_params){
    this->batch_size = batch_size;
    this->max_meas = max_meas;
    this->meas_strat = meas_strat;
    if(log) {
        struct stat sb;
        if (stat(log_base_path.c_str(), &sb) == 0) {
            std::cout << "Log directory already exists. Exiting." << std::endl;
            exit(1);
        } else {
            mkdir(log_base_path.c_str(), 0777);
            std::cout << "Created Experiment Directory in " << log_base_path
                      << std::endl;
        }
        if (!stat((log_base_path + "/configs").c_str(), &sb) == 0) {
            mkdir((log_base_path + "/configs").c_str(), 0777);
        }
    }
    this->visualize = visualize;
    this->log = log;
    this->scenario = scenario;

    // Parse the correct map size
    this->gmrf_params->size_x = scenario->get_sx();
    this->gmrf_params->size_y = scenario->get_sy();


    // Initialize the obstacles
    auto objects = scenario->getObjects()->objects;
    auto walls = objects[wallType::wall];
    auto inlets = objects[wallType::inlet];
    auto outlets = objects[wallType::outlet];
    auto sources = objects[wallType::source];

    gmrf = new GMRF(*this->gmrf_params, this->scenario);
    if(visualize) {
        ObjectsViz *obj_viz = new ObjectsViz(scenario->getObjects());
        Viz_Params viz_params;
        viz_params.offset_y = 100;
        viz = new RealtimeViz(gmrf, scenario->getField(),
                              scenario->getFlowFieldX(),
                              scenario->getFlowFieldY(), viz_params,
                              *gmrf_params);
        viz->addVizObject(obj_viz, map_types::gt_map);
        viz->addVizObject(obj_viz, map_types::field_map);
        gmrf->addToViz(*viz, 0, 4, 0, 5);
    }
    for(int i = 0; i < walls.size(); i+=2){
        rrt_obstacles.push_back(new RRTObstacle(walls[i].position.x, walls[i]
                .position
                .y, walls[i + 1].position.x, walls[i + 1].position.y));
    }
    for(int i = 0; i < inlets.size(); i+=2){
        rrt_obstacles.push_back(new RRTObstacle(inlets[i].position.x, inlets[i]
                .position
                .y, inlets[i + 1].position.x, inlets[i + 1].position.y));
    }
    for(int i = 0; i < outlets.size(); i+=2) {
        rrt_obstacles.push_back(
                new RRTObstacle(outlets[i].position.x, outlets[i]
                                        .position
                                        .y, outlets[i + 1].position.x,
                                outlets[i + 1].position.y));
    }
    for(int i = 0; i < sources.size(); i+=2) {
        rrt_obstacles.push_back(
                new RRTObstacle(sources[i].position.x, sources[i]
                                        .position
                                        .y, sources[i + 1].position.x,
                                sources[i + 1].position.y));
    }
    timer.addModuleToLog("gmrf");
    timer.addModuleToLog("source_localization");
    // Prepare Logger
    measurement_locations.resize(0, 2);
    if(this->log) {
        // open base config file to save scenario name
        std::ofstream config_file;
        config_file.open(log_base_path + "/config.txt");
        config_file << scenario->getName();
        config_file.close();
        gmrf_params->exportConfig(log_base_path);
        logger = new SimLogger(log_base_path);
        logger->addLoggable(&timer);
        logger->addLoggable(new LoggableEigenVector("gmrf_estimates",
                                                    &gmrf->getAllEstimates()));
        logger->addLoggable(new LoggableEigen("source_estimates",
                                              &gmrf->getSourceLikelihoods()));
        logger->addLoggable(new LoggableEigen("field_covariances",
                                              &gmrf->getCovariances()));
        logger->addLoggable(new LoggableEigen("flowx_covariances",
                                              &gmrf->getCovariances
                                                      (measurement_type::flow_x)));
        logger->addLoggable(new LoggableEigen("flowy_covariances",
                                              &gmrf->getCovariances
                                                      (measurement_type::flow_x)));
        logger->addLoggable(new LoggableDynamicEigen("measurement_locations",
                                              &measurement_locations));
        logger->addLoggable(new LoggablePoint("Source", &gmrf->source_x, &gmrf->source_y));
        logger->addMetric(new RMSEMetric(gmrf, scenario, rrt_obstacles,0.2 ));
        logger->addMetric(new FieldLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new FieldLogLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new SourceDistMetric(gmrf, scenario));
    }
}

bool GMRFInterface::step() {
    static thread_local std::mt19937 generator(std::random_device{}());
    if(log) {
        logger->log(t);
    }

    // Sample random x and y uniform
    double x = 0;
    double y = 0;
    for(int k = 0; k < batch_size; k++) {
        measurement_locations.conservativeResize(
                measurement_locations.rows() + 1, 2);
        if (meas_strat == MeasurementStrategy::RandomUniform) {
            // distribution
            std::uniform_real_distribution<double> distribution_x(0,
                                                                  gmrf_params->size_x);
            std::uniform_real_distribution<double> distribution_y(0,
                                                                  gmrf_params->size_y);
            x = distribution_x(generator);
            y = distribution_y(generator);
        } else if (meas_strat == MeasurementStrategy::RandomizedSwipe) {
            double mean =
                    (gmrf_params->size_x - 1.0) * (max_meas - t) / max_meas +
                    1.0;
            std::uniform_real_distribution<double> distribution_x(mean - 1.0,
                                                                  mean + 1.0);
            std::uniform_real_distribution<double> distribution_y(0,
                                                                  gmrf_params->size_y);
            x = distribution_x(generator);
            y = distribution_y(generator);
        } else if (meas_strat == MeasurementStrategy::SourceLocalization){
            // Sample next measurement location from the discrete source likelihoods
            Eigen::MatrixXd source_likelihoods = gmrf->getSourceLikelihoods();
//            source_likelihoods =source_likelihoods.array()*0.5 +  0.5*source_likelihoods.maxCoeff();
            source_likelihoods = source_likelihoods.array() / source_likelihoods.sum();
            std::discrete_distribution<int> distribution(source_likelihoods.data(),
                                                         source_likelihoods.data() +
                                                         source_likelihoods.size());
            int idx = distribution(generator);

            x = (double)(idx % gmrf_params->N_X) * gmrf->node_dist[0];
            y = (int)(idx / gmrf_params->N_X) * gmrf->node_dist[1];
        } else if (meas_strat == MeasurementStrategy::UncertaintyReduction){
            Eigen::MatrixXd uncertainties = gmrf->getCovariances();
            uncertainties = uncertainties.array() / uncertainties.sum();
            std::discrete_distribution<int> distribution(uncertainties.data(),
                                                         uncertainties.data() +
                                                         uncertainties.size());
            int idx = distribution(generator);
            x = (double)(idx % gmrf_params->N_X) * gmrf->node_dist[0];
            y = (int)(idx / gmrf_params->N_X) * gmrf->node_dist[1];
        }
        measurement_locations(measurement_locations.rows() - 1, 0) = x;
        measurement_locations(measurement_locations.rows() - 1, 1) = y;
        gmrf->addMeasurement(x, y, std::min(1.0, std::max(0.0, scenario->getFieldValue(x, y, 0.05))), t,
                             field);
        double flow_x_meas;
        double flow_y_meas;
        scenario->getFlowValue(x, y, flow_x_meas, flow_y_meas, 0.0);
        gmrf->addMeasurement(x, y, flow_x_meas, t, flow_x);
        gmrf->addMeasurement(x, y, flow_y_meas, t, flow_y);
    }
    double x_;
    double y_;
    bool fnd;
    timer.startTimer("gmrf");
    gmrf->updateEstimates(t);
    timer.stopTimer("gmrf");
    timer.startTimer("source_localization");
    gmrf->sourceHypothesisNeighbor(x_, y_,
                                   fnd);

    timer.stopTimer("source_localization");
    if(visualize){
        draw();
    }
    t += batch_size;
    return true;
}

bool GMRFInterface::done() {
    return t > max_meas;
}

void GMRFInterface::postProcess() {
    if (!log) {
        // calculate final metrics
        std::cout << "RMSE: " << RMSEMetric(gmrf, scenario, rrt_obstacles, 1.0).updateMetric()
                  << std::endl;
        std::cout << "Field Likelihood: "
                    << FieldLikelihoodMetric(gmrf, scenario, 0.05).updateMetric()
                    << std::endl;
        std::cout << "Field Log Likelihood: "
                    << FieldLogLikelihoodMetric(gmrf, scenario, 0.05).updateMetric()
                    << std::endl;
        std::cout << "Source Distance: "
                    << SourceDistMetric(gmrf, scenario).updateMetric()
                    << std::endl;
        return;
    }
}

void GMRFInterface::draw() {
    viz->draw();
}

GMRFInterface::GMRFInterface(std::string replay_path) {
    visualize = true;

    // Load the scenario.
    std::ifstream config_file;
    config_file.open(replay_path + "/config.txt");
    std::string scenario_name;
    config_file >> scenario_name;
    scenario = new ConfigFileScenario(scenario_name);

    // Load the module params
    gmrf_params = new GmrfParams(replay_path + "/configs/gmrf_conf.toml");

    // Load the modules
    gmrf = new GMRF(*this->gmrf_params, this->scenario);

    // Viz
    ObjectsViz* obj_viz = new ObjectsViz(scenario->getObjects());
    Viz_Params viz_params;
    viz_params.offset_y = 100;
    viz = new RealtimeViz(gmrf, scenario->getField(), scenario->getFlowFieldX(), scenario->getFlowFieldY(), viz_params, *gmrf_params);
    viz->addVizObject(obj_viz,  map_types::gt_map);
    viz->addVizObject(obj_viz,  map_types::field_map);
    viz->addVizObject(obj_viz,  map_types::std_map);
    gmrf->addToViz(*viz, 0, 4, 0, 5);
    // Load the states
    std::ifstream state_file;
    std::ifstream gmrf_file;
    EigenVectorReader gmrf_reader(replay_path + "/gmrf_estimates.csv");
    EigenMatrixReader field_cov_reader(replay_path + "/field_covariances.csv");
    EigenMatrixReader flowx_cov_reader(replay_path + "/flowx_covariances.csv");
    EigenMatrixReader flowy_cov_reader(replay_path + "/flowy_covariances.csv");


    // Start reading the actual data.
    while (true) {
        bool success = true;
        if(!success){
        }else {
            t+=1;
            gmrf->setEstimates(gmrf_reader.readVector(success), t);
            gmrf->setCovariance(field_cov_reader.readMatrix(success), t);
            gmrf->setCovariance(flowx_cov_reader.readMatrix(success), t, measurement_type::flow_x);
            gmrf->setCovariance(flowy_cov_reader.readMatrix(success), t,  measurement_type::flow_y);
        }
        draw();
    }
}

GMRFInterface::~GMRFInterface() {
    std::cout << "Deleting GMRFInterface" << std::endl;
    delete gmrf;

    if (visualize) {
        delete viz;
    }
    if (log) {
        delete logger;
    }

}

BFSInterface::BFSInterface(std::string log_base_path, GmrfParams *gmrf_params,
                           ConfigFileScenario *scenario, bool visualize,
                           bool log, double x0, double y0) {
    if(log) {
        struct stat sb;
        if (stat(log_base_path.c_str(), &sb) == 0) {
            std::cout << "Log directory already exists. Exiting." << std::endl;
            exit(1);
        } else {
            mkdir(log_base_path.c_str(), 0777);
            std::cout << "Created Experiment Directory in " << log_base_path
                      << std::endl;
        }
        if (!stat((log_base_path + "/configs").c_str(), &sb) == 0) {
            mkdir((log_base_path + "/configs").c_str(), 0777);
        }
    }
    this->visualize = visualize;
    this->log = log;
    this->scenario = scenario;

    this->gmrf_params = gmrf_params;

    current_state.resize(2, 1);
    current_state.setZero();
    // Parse the correct map size
    this->gmrf_params->size_x = scenario->get_sx();
    this->gmrf_params->size_y = scenario->get_sy();

    // Initialize the obstacles
    auto objects = scenario->getObjects()->objects;
    auto walls = objects[wallType::wall];
    auto inlets = objects[wallType::inlet];
    auto outlets = objects[wallType::outlet];
    auto sources = objects[wallType::source];

    for(int i = 0; i < walls.size(); i+=2){
        obstacles.push_back(new RRTObstacle({walls[i].position.x, walls[i]
                .position
                .y, walls[i + 1].position.x, walls[i + 1].position.y}));
    }
    for(int i = 0; i < inlets.size(); i+=2){
        obstacles.push_back(new RRTObstacle({inlets[i].position.x, inlets[i]
                .position
                .y, inlets[i + 1].position.x, inlets[i + 1].position.y}));
    }
    for(int i = 0; i < outlets.size(); i+=2) {
        obstacles.push_back(
                new RRTObstacle({outlets[i].position.x, outlets[i]
                        .position
                        .y, outlets[i + 1].position.x,
                                        outlets[i + 1].position.y}));
    }
    for(int i = 0; i < sources.size(); i+=2) {
        obstacles.push_back(
                new RRTObstacle({sources[i].position.x, sources[i]
                        .position
                        .y, sources[i + 1].position.x,
                                        sources[i + 1].position.y}));
    }
    gmrf = new GMRF(*this->gmrf_params, this->scenario);
    bfs = new BFSPlanner(&obstacles, gmrf_params->N_X, gmrf_params->N_Y, x0, y0, gmrf_params->size_x, gmrf_params->size_y);

    if(visualize){
        ObjectsViz* obj_viz = new ObjectsViz(scenario->getObjects());
        Viz_Params viz_params;
        viz_params.offset_y = 100;
        viz = new RealtimeViz(gmrf, scenario->getField(), scenario->getFlowFieldX(), scenario->getFlowFieldY(), viz_params, *gmrf_params);
        viz->addVizObject(obj_viz,  map_types::gt_map);
        viz->addVizObject(obj_viz,  map_types::field_map);
        gmrf->addToViz(*viz, 0, 4, 0, 5);
        agent_position = new sf::ConvexShape(3);
        agent_position->setPoint(0, sf::Vector2f(0, 0.3));
        agent_position->setPoint(1, sf::Vector2f(0.0, -0.3));
        agent_position->setPoint(2, sf::Vector2f(0.6, 0.0));
        agent_position->setOutlineColor(sf::Color::Black);
        agent_position->setOutlineThickness(0.07);
        viz->addVizObject(agent_position, map_types::std_map);
    }
    // Prepare Logger
    if(this->log) {
        // open base config file to save scenario name
        std::ofstream config_file;
        config_file.open(log_base_path + "/config.txt");
        config_file << scenario->getName();
        config_file.close();
        gmrf_params->exportConfig(log_base_path);
        logger = new SimLogger(log_base_path);
        logger->addLoggable(new LoggableEigenVector("state", &current_state));
        logger->addLoggable(new LoggableEigenVector("gmrf_estimates",
                                                    &gmrf->getAllEstimates()));
        logger->addLoggable(new LoggableEigen("source_estimates",
                                              &gmrf->getSourceLikelihoods()));
        logger->addLoggable(new LoggableEigen("field_covariances",
                                              &gmrf->getCovariances()));
        logger->addLoggable(new LoggableEigen("flowx_covariances",
                                              &gmrf->getCovariances
                                                      (measurement_type::flow_x)));
        logger->addLoggable(new LoggableEigen("flowy_covariances",
                                              &gmrf->getCovariances
                                                      (measurement_type::flow_x)));
        logger->addLoggable(new LoggablePoint("Source", &gmrf->source_x, &gmrf->source_y));
        logger->addMetric(new RMSEMetric(gmrf, scenario, obstacles, 0.2));
        logger->addMetric(new FieldLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new FieldLogLikelihoodMetric(gmrf, scenario, 0.2));
        logger->addMetric(new SourceDistMetric(gmrf, scenario));

    }
}

bool BFSInterface::done() {
    return done_;
}

bool BFSInterface::step() {
    if(log) {
        logger->log(t);
    }

    t += 1;
    current_step++;
    // Sample random variation of x and y within 0.5m using mt19937
    static std::random_device rd;
    static std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-0.5, 0.5);
    double x = dist(mt) + bfs->getX();
    double y = dist(mt) + bfs->getY();
    x = std::max(x, 0.0);
    x = std::min(x, gmrf_params->size_x);
    y = std::max(y, 0.0);
    y = std::min(y, gmrf_params->size_y);
    double x_c = (x + last_x)/2;
    double y_c = (y + last_y)/2;
    double v = scenario->getFieldValue(x, y, 0.05);
    double v_c = scenario->getFieldValue(x_c, y_c, 0.05);
    gmrf->addMeasurement(x, y, v, t,
                         field);
    gmrf->addMeasurement(x_c, y_c, v_c, t,
                            field);
    double flow_x_meas;
    double flow_y_meas;
    double flow_x_meas_c;
    double flow_y_meas_c;

    scenario->getFlowValue(x, y, flow_x_meas, flow_y_meas, 0.05);
    scenario->getFlowValue(x_c, y_c, flow_x_meas_c, flow_y_meas_c, 0.05);
    gmrf->addMeasurement(x, y, flow_x_meas, t, flow_x);
    gmrf->addMeasurement(x, y, flow_y_meas, t, flow_y);
    gmrf->addMeasurement(x_c, y_c, flow_x_meas_c, t, flow_x);
    gmrf->addMeasurement(x_c, y_c, flow_y_meas_c, t, flow_y);
    double x_ = gmrf->source_x;
    double y_ = gmrf->source_y;
    bool success = false;
    if((int)round(t) % 3 ==0){
        gmrf->updateEstimates(t);
        gmrf->sourceHypothesisNeighbor(x_, y_, success);
    }
    double dx = gmrf_params->size_x/(gmrf_params->N_X - 1.0);
    double dy = gmrf_params->size_y/(gmrf_params->N_Y - 1.0);
    bfs->step(x_/dx, y_/dy);
    current_state(0, 0) = bfs->getX();
    current_state(1, 0) = bfs->getY();
    if(current_state(0,0) == last_x && current_state(1,0) == last_y){
        num_repeats++;
    }
    else {
        num_repeats = 0;
    }
    if (num_repeats >= 5){
        done_ = true;
    }
    last_x = current_state(0,0);
    last_y = current_state(1,0);
    if(visualize){
        agent_position->setPosition(current_state
                                                (0, 0), current_state(1 , 0));
        draw();
    }
    return true;
}

void BFSInterface::postProcess() {

}

void BFSInterface::draw() {
    viz->draw();
}

BFSInterface::~BFSInterface() {
    std::cout << "Deleting BFSInterface" << std::endl;
    delete gmrf;
    delete bfs;
    for (auto &obj : obstacles) {
        delete obj;
    }

    if (visualize) {
        delete viz;
        delete agent_position;
    }
    if (log) {
        delete logger;
    }

}
