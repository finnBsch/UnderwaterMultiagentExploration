//
// Created by finn on 7/2/23.
//

#ifndef MASTERMAIN_SIM_INTERFACE_H
#define MASTERMAIN_SIM_INTERFACE_H

#include "rrt.h"
#include "gmrf/gmrf.h"
#include "sqp_multiagent.h"
#include "scenario.h"
#include "field_interface.h"
#include "obstacle_interface.h"
#include "rrt_interface.h"
#include "logger.h"
#include "metrics.h"
#include "random_walk.h"
#include "brownian_random_walk.h"
#include "ballistic_motion.h"
#include "bfs_planner.h"

/**
 * Abstract sim interface class. Implements some basic functionality that all
 * sims need.
 */
class SimInterface{
protected:
    SimLogger* logger;
    bool visualize = false;  // Whether to visualize the simulation
    bool log = false;  // Whether to run in real time or step in fixed
                             // time steps as provided by the sqp params.
public:
    virtual bool step() = 0;
    virtual bool done() = 0;
    virtual void postProcess() = 0;
    virtual void draw() = 0;
    bool runFull();
};

enum class StoppingCritereon{
    N_steps = 0,
    Time = 1,
    SourceFound = 2,
};
enum class BenchmarkPlanner{
    Brownian = 0,
    Ballistic = 1,
};
enum class MeasurementStrategy{
    RandomUniform = 0,
    RandomizedSwipe = 1,
    SourceLocalization= 2,
    UncertaintyReduction = 3,
};

class BenchmarkInterface : public SimInterface {
private:
    std::string log_base_path;
    bool found_hypothesis = false;
    double t = 0;
    std::vector<RandomWalkObstacle*> obstacles;
    std::vector<RRTObstacle*> rrt_obstacles;
    GmrfParams *gmrf_params;
    SQPParams *sqp_params;
    ConfigFileScenario *scenario;
    ObstacleInterface *obs_interface;
    BenchmarkPlanner planner_type;

    GMRF *gmrf;
    RandomWalk* random_walk;
    // Viz
    std::vector<sf::ConvexShape *> agent_positions;
    RealtimeViz *viz;

    Eigen::Matrix<double, Eigen::Dynamic, 1> current_state;
    int current_step = 0;
    double source_x = 0;
    double source_y = 0;
    int num_agents;
public:
    BenchmarkInterface(BenchmarkPlanner planner_type, std::string log_base_path, GmrfParams *gmrf_params,  SQPParams* sqp_params,
                       ConfigFileScenario *scenario, bool visualize, bool log, double x0, double y0);
    bool done() override;
    bool step() override;

    void postProcess() override;
    void draw() override;
    ~BenchmarkInterface();
};

class BFSInterface : public SimInterface {
private:

    double last_x = 0;
    double last_y = 0;
    int num_repeats = 0;
    double t = 0;
    std::vector<RRTObstacle*> obstacles;
    GmrfParams *gmrf_params;
    ConfigFileScenario *scenario;
    ObstacleInterface *obs_interface;
    bool done_ = false;
    GMRF *gmrf;
    BFSPlanner* bfs;
    // Viz
    sf::ConvexShape * agent_position;
    RealtimeViz *viz;

    Eigen::Matrix<double, Eigen::Dynamic, 1> current_state;
    int current_step = 0;

    int num_agents;
public:
    BFSInterface(std::string log_base_path, GmrfParams *gmrf_params,
                 ConfigFileScenario *scenario, bool visualize, bool log, double x0, double y0);
    bool done() override;
    bool step() override;

    void postProcess() override;
    void draw() override;

    ~BFSInterface();
};


/**
 * All modules combined in closed-loop.
 */
class FullInterface : public SimInterface{
private:
    std::string log_base_path;
    bool found_hypothesis = false;
    SQPParams* sqp_params;
    RRTParams* rrt_params;
    GmrfParams* gmrf_params;

    // Interfaces
    ConfigFileScenario* scenario;
    ObstacleInterface* obs_interface;
    FieldInterface* field_interface;
    RRTSplinePathInterface* rrt_spline_interface;
    std::vector<RRTObstacle*> rrt_obstacles;

    // Modules
    GMRF* gmrf;
    RRT* rrt;
    SQPMultiAgent* sqp;

    // Extra Viz
    std::vector<sf::ConvexShape*> agent_positions;
    RealtimeViz* viz;

    // Runtime members
    Eigen::Matrix<double, Eigen::Dynamic, 1> current_state;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> trajectory_state;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> trajectory_input;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> measurement_locations;
    int current_step = 0;
    Node* goal_node = nullptr;
    double source_x = 0;
    double source_y = 0;
    double t = 0;
    double replan_angle = 0;
    LoggableTimer timer;
public:
    FullInterface(std::string log_base_path, SQPParams* sqp_params, RRTParams*
    rrt_params,
                  GmrfParams*
    gmrf_params, ConfigFileScenario* scenario, StoppingCritereon stop_critereon,
    bool visualize, bool log, double x0, double y0);
    FullInterface(std::string replay_path);
    bool done();
    bool step();


    void postProcess() override;
    void draw() override;
    ~FullInterface();
};

/**
 * Just the GMRF
 */
class GMRFInterface : public SimInterface{
private:
    int batch_size = 1;
    ConfigFileScenario* scenario;
    GmrfParams* gmrf_params;
    MeasurementStrategy meas_strat;
    GMRF* gmrf;
    LoggableTimer timer;
    int max_meas;
    int t = 0;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> measurement_locations;
    std::vector<RRTObstacle*> rrt_obstacles;
public:
    bool step() override;

    bool done() override;

    void postProcess() override;

    void draw() override;

    ~GMRFInterface();

private:
    // Viz
    RealtimeViz* viz;
public:
    GMRFInterface(std::string log_base_path, GmrfParams* gmrf_params, ConfigFileScenario* scenario, MeasurementStrategy meas_strat, bool visualize, bool log, int max_meas, int batch_size);
    GMRFInterface(std::string replay_path);
};

/**
 * Just the planning (RRT + SQP)
 */
class PlanningInterface : public SimInterface {
private:
    SQPParams* sqp_params;
    RRTParams* rrt_params;

    RRT* rrt;
    SQPMultiAgent* sqp;
public:
    PlanningInterface(SQPParams* sqp_params, RRTParams* rrt_params, bool visualize);
};


/**
 * Just the RRT
 */
class RRTInterface : public SimInterface {
private:
    RRTParams* rrt_params;

    RRT* rrt;
public:
    RRTInterface(RRTParams* rrt_params, bool visualize);

};

class SQPInterface : public SimInterface {
private:
    SQPParams* sqp_params;

    SQPMultiAgent* sqp;
public:
    SQPInterface(SQPParams* sqp_params, bool visualize);
};


class EigenVectorReader {
private:
    double t;
    std::ifstream file;
    Eigen::VectorXd data;
public:
    EigenVectorReader(std::string path);
    Eigen::VectorXd& readVector(bool & success);
    int getRows();
    double getTime();
    bool eof();
};

class EigenMatrixReader {
private:
    std::ifstream file;
    Eigen::MatrixXd data;
public:
    EigenMatrixReader(std::string path);
    Eigen::MatrixXd& readMatrix(bool & success);
};
#endif //MASTERMAIN_SIM_INTERFACE_H
