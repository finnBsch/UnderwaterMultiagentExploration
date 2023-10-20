//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
// Implementation of https://link.springer.com/content/pdf/10.1007/s10514-015-9437-0.pdf?pdf=button
// Infos on covariance calculation using upper triagonal cholesky: https://mathoverflow.net/questions/163414/numerical-trace-of-inverse-matrix-from-cholesky
// Reference implementation https://github.com/MAPIRlab/gmrf_gas_mapping

#ifndef GMRF_GMRF_H
#define GMRF_GMRF_H

#include <iostream>
#include <stdexcept>
#include <unordered_map>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseCore>

#include <eigen3/Eigen/SparseCholesky>

#include <SFML/Graphics.hpp>

#include "gmrf_viz.h"
#include "settings.h"
#include "gmrf_configs.h"

#include <set>
#include <random>

class GMRF;
typedef Eigen::Map<const Eigen::Matrix<double, Dynamic, Dynamic>> VecToMat;
typedef Eigen::Map<Eigen::Matrix<double, Dynamic, 1>> EstToOneMap;

enum measurement_type{
    field = 0,
    flow_x = 1,
    flow_y = 2
};

struct Edge{
    int N0;
    int N1;
    bool horiz;
    int x0;
    int y0;
    int x1;
    int y1;
    bool is_cut = false;
    double min_cut_pos = std::numeric_limits <double>::infinity();
    double min_cut_max_depth = -std::numeric_limits <double>::infinity();
    double min_cut_min_depth = std::numeric_limits <double>::infinity();
    double max_cut_pos = -std::numeric_limits <double>::infinity();
    double max_cut_max_depth = -std::numeric_limits <double>::infinity();
    double max_cut_min_depth = std::numeric_limits <double>::infinity();
    friend std::ostream& operator <<(std::ostream& os, const Edge& edge);

};
class GMRFNode{
    friend class Edge;
    friend class GMRF;
private:

public:
    int node_id;
    GMRFNode();
    int num_meas = 0;
    std::array<double, 2> dists;
    std::vector<Edge*> connected_edges;
    int x;
    int y;
    bool x_wall = false;
    bool y_wall = false;
    bool source = false;
    bool reachable = false;

    void addMeasurement(double& x, double &y, double &v, double &s, double & t);
};

struct Obstacle {
    double x0;
    double x1;
    double y0;
    double y1;
};

class EnergyFunc{
protected:
    GMRF* gmrf;
public:
    virtual ~EnergyFunc() = 0;
    int num_energies;
    std::array<int, 2> idx;
    virtual void initializeJ() = 0;
    virtual void initializeInfomat() = 0;
    virtual void updateJ() = 0;
    virtual void updateF() = 0;
};

class FieldBoundary_EnergyFunc : public EnergyFunc {
private:
public:
    explicit FieldBoundary_EnergyFunc(GMRF* gmrf, int idx0);
    void initializeJ() override;
    void initializeInfomat() override;
    void updateJ() override;
    void updateF() override;
};


class FieldN_EnergyFunc : public EnergyFunc {
private:
public:
    explicit FieldN_EnergyFunc(GMRF* gmrf, int idx0);
    void initializeJ() override;
    void initializeInfomat() override;
    void updateJ() override;
    void updateF() override;
};

class FlowN_EnergyFunc : public EnergyFunc {
private:
public:
    explicit FlowN_EnergyFunc(GMRF* gmrf, int idx0);
    void initializeJ() override;
    void initializeInfomat() override;
    void updateJ() override;
    void updateF() override;
};

class FlowC_EnergyFunc: public EnergyFunc{
private:
public:
    explicit FlowC_EnergyFunc(GMRF* gmrf, int idx0);
    void initializeJ() override;
    void initializeInfomat() override;
    void updateJ() override;
    void updateF() override;
};

class FlowField_EnergyFunc: public EnergyFunc{
private:
public:
    explicit FlowField_EnergyFunc(GMRF* gmrf, int idx0);
    void initializeJ() override;
    void initializeInfomat() override;
    void updateJ() override;
    void updateF() override;
};

class ObstacleFlow_EnergyFunc: public EnergyFunc{
private:
public:
    explicit ObstacleFlow_EnergyFunc(GMRF* gmrf, int idx0);
    void initializeJ() override;
    void initializeInfomat() override;
    void updateJ() override;
    void updateF() override;
};


class GMRF: public DrawableObject {
friend class EnergyFunc;
friend class FieldN_EnergyFunc;
friend class FieldBoundary_EnergyFunc;
friend class FlowN_EnergyFunc;
friend class FlowC_EnergyFunc;
friend class FlowField_EnergyFunc;
friend class ObstacleFlow_EnergyFunc;
private:

private:
    // Temporary
    double max_est = 1;
    double min_est = 0;
    double max_cov = 1;
    double min_cov = 0;

    // Algorithm Attributes
    GmrfParams params;
    int num_meas = 0;
    int num_nodes;
    int num_edges;
    int num_energies;  // Number of energies without measurements
    std::vector<EnergyFunc*> energy_funs;



    std::unordered_map<measurement_type, Eigen::Matrix<double, Dynamic, Dynamic>> covariances;
    std::vector<measurement_type> meas_types;
    Eigen::Matrix<double, Dynamic, 1> measurement_age;
    Eigen::Matrix<double, Dynamic, 1> f;
    Eigen::Matrix<double, Dynamic, 3> y_meas;
    Eigen::Matrix<double, Dynamic, 1> all_estimates;
    Eigen::Matrix<double, Dynamic, Dynamic> source_likelihoods;
    std::unordered_map<measurement_type, EstToOneMap> estimates;
    std::vector<Edge> edges;
    Eigen::MatrixXd Sigma;
    std::vector<GMRFNode> nodes;
    double sigma_p = 0.0;
    double ff_sigma_p = 0.0;

    SparseMatrix<double> information_matrix;
    SimplicialLLT<SparseMatrix<double>> solver;
    std::vector<int> parameterized_edges;
    SparseMatrix<double> J;

    // Viz Attributes
    bool viz = false;
    sf::VertexArray lines;
    RealtimeViz* viz_obj;
    sf::CircleShape source;


    // Runtime Methods
//    void updateParameters(double &t);
    void updateInfomat(double t);
    void update_J();
    void updateCorrelationInfomat(int L_id, double prob);
    void updateMeasurementInfomat(double t);
    void updateF();
    void updateViz(double t);
    void determineReachableNodes();
    void checkChildren(GMRFNode* node, std::set<GMRFNode*>& checked_nodes);


    // Initialization methods
    void initializeInformationMatrix();
    void initializeJ();
    void initializeNodes();
    void initializeEdges();
    void initializePrior();

    // Helper Methods
    std::array<double, 6> getShapeFacs(double& x, double& y, bool normalized
    = false) const;

    inline int fieldEstimateID(int x, int y) const{
        return xyToN(x, y);
    }
    inline int flowEstimateID(int x, int y) const{
        return num_nodes + xyToN(x, y);
    }
    int nodepairToL(int x0, int x1, int y0, int y1) const;

public:
    GMRF(GmrfParams params, Scenario* scenario);

    // Runtime Methods
    double interpolateValue(double x, double y, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &field) const;
    std::array<double, 2> uncertaintyGradient(double x, double y, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& field) const;
    void updateEstimates(double t);
    void addMeasurement(double x, double y, double v, double t, measurement_type meas_t);
    double getEstimate(double x, double y) const;
    double getEstimate(double x, double y, measurement_type m_type) const;
    double getLikelihood(double x, double y, double v) const;
    double getLogLikelihood(double x, double y, double v) const;
    double getCovariance(double x, double y) const;
    double getCovariance(double x, double y, measurement_type m_type) const;

    void sourceHypothesisNeighbor(double& x_s, double& y_s, bool& success);

    // Obstacle
    void addLineObstacle(double x0, double y0, double x1, double y1);

    // Viz Methods
    VecToMat getEstimates() const;
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&
            getSourceLikelihoods() const;
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&
    getCovariances() const;
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&
    getCovariances(measurement_type meas_type) const;
    Eigen::Matrix<double, Dynamic, Dynamic> estimates_map;
    Eigen::Matrix<double, Dynamic, Dynamic> stds_map;
    Eigen::Matrix<double, Dynamic, Dynamic> f_stds_map;
    void addToViz(RealtimeViz& viz_obj, double min_est, double max_est, double min_std, double max_std);
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    void prepareViz() override;
    void setEstimates(const Eigen::Matrix<double, Eigen::Dynamic, 1>&
            estimates, double t);
    void setCovariance(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&
            covariances, double t, measurement_type meas_type = measurement_type::field);

    // Temporary
    std::vector<sf::CircleShape*> particle_viz;
    std::vector<sf::CircleShape*> source_viz;

    // Getter Methods
    const EstToOneMap* getFlowXMap() const;
    const EstToOneMap* getFlowYMap() const;

    // Log and init
    const Eigen::Matrix<double, Dynamic, 1>& getAllEstimates() const;

    bool near_source = false;
    bool isReachable(double x, double y) const;
    // Source
    double source_x = 0;
    double source_y = 0;
    bool source_found = false;

    // Destructor
    ~GMRF();
    inline int xyToN(int x, int y) const{
        return y * params.N_X + x;
    }
    std::array<double, 2> node_dist;
};

#endif //GMRF_GMRF_H
