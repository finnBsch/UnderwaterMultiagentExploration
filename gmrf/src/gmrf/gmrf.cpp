//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//

#include "gmrf/gmrf.h"
#include <Eigen/Eigenvalues>
#include "msg_assert.h"
#include <chrono>

GMRF::GMRF(GmrfParams params, Scenario* scenario): params(params), source(0.3){
    source.setOrigin(0.3, 0.3);
    sigma_p = params.sigma_p_per_m * pow(params.size_x / (params.N_X - 1), 2);
    ff_sigma_p = params.ff_sigma_p_per_m * pow(params.size_x / (params.N_X - 1), 2);
    // Problem Dimension.
    num_nodes = params.N_X * params.N_Y;
    num_edges = (params.N_X - 1) * params.N_Y + (params.N_Y - 1) * params.N_X;
    Sigma = MatrixXd(num_nodes * 3,num_nodes * 3);
    num_energies = 0;
    if(params.field_neighbors){
        energy_funs.push_back(new FieldN_EnergyFunc(this, num_energies));
        num_energies += energy_funs[energy_funs.size() - 1]->num_energies;
    }
    if(params.flow_neighbors){
        energy_funs.push_back(new FlowN_EnergyFunc(this, num_energies));
        num_energies += energy_funs[energy_funs.size() - 1]->num_energies;
    }
    if(params.flow_conservation){
        energy_funs.push_back(new FlowC_EnergyFunc(this, num_energies));
        num_energies += energy_funs[energy_funs.size() - 1]->num_energies;
    }
    if(params.flow_field){
        energy_funs.push_back(new FlowField_EnergyFunc(this, num_energies));
        num_energies += energy_funs[energy_funs.size() - 1]->num_energies;
    }
    if(params.obstacle_flow){
        energy_funs.push_back(new ObstacleFlow_EnergyFunc(this, num_energies));
        num_energies += energy_funs[energy_funs.size() - 1]->num_energies;
    }
//    energy_funs.push_back(new FieldBoundary_EnergyFunc(this, num_energies));
//    num_energies += energy_funs[energy_funs.size() - 1]->num_energies;
//    num_energies = num_edges // Neighboring concentrations
//            + num_edges * 2  // Neighboring flows
//            + num_nodes  // Conservation of Flow Mass
//            + num_nodes  // Influence of SQPObstacle on Flow Vector TODO Can't this be implicitly in the conservation?!
//            + num_nodes;  // Flow-Quantity relationship

//    std::cout << "GMRF Problem Dimensions: \n\tNumber of Nodes: " << num_nodes << "\n\tNumber of Edges: " << num_edges <<
//              "\n\tNumber of Energies (not including measurements): " << num_energies << "\n";
    for(int i = 0; i < num_nodes; i++){
        nodes.emplace_back();
    }
    for(int i = 0; i < num_edges; i++){
        edges.emplace_back();
    }
    node_dist[0] = params.size_x / (double)(params.N_X - 1);
    node_dist[1] = params.size_y / (double)(params.N_Y - 1);
    all_estimates.resize(num_nodes * 3, 1);

    estimates.emplace(field, EstToOneMap(all_estimates.data(), num_nodes));

    estimates.emplace(flow_x, EstToOneMap(all_estimates.data() + num_nodes, num_nodes));

    estimates.emplace(flow_y, EstToOneMap(all_estimates.data() + 2 * num_nodes, num_nodes));

    covariances.emplace(field, Matrix<double, Dynamic, Dynamic>());
    covariances[field].resize(params.N_X, params.N_Y);
    covariances[field].setZero();

    covariances.emplace(flow_x, Matrix<double, Dynamic, Dynamic>());
    covariances[flow_x].resize(params.N_X, params.N_Y);
    covariances[field].setZero();

    covariances.emplace(flow_y, Matrix<double, Dynamic, Dynamic>());
    covariances[flow_y].resize(params.N_X, params.N_Y);
    covariances[field].setZero();

    source_likelihoods.resize(params.N_X, params.N_Y);
    source_likelihoods.setZero();

    estimates_map.resize(MAP_X, MAP_Y);
    stds_map.resize(MAP_X, MAP_Y);
    f_stds_map.resize(params.N_X, params.N_Y);

    // Resize information matrix to num_energies
    information_matrix.resize(num_energies, num_energies);
    f.resize(num_energies, f.cols());
    y_meas.resize(num_energies, y_meas.cols());

    // Default to zero
    y_meas.setZero();
    f.setZero();
    all_estimates.setZero();

    // Fill initial information matrix
    initializeNodes();
    initializeEdges();
    initializeInformationMatrix();
    auto walls = scenario->getObjects()->objects[wallType::wall];
    for(int i = 0; i < walls.size(); i+=2){
        addLineObstacle(walls[i].position.x, walls[i]
                .position
                .y, walls[i + 1].position.x, walls[i + 1].position.y);
    }
    determineReachableNodes();
    auto sources = scenario->getObjects()->objects[wallType::source];
    for(int i = 0; i < sources.size(); i+=2){
        double x0 = sources[i].position.x;
        double y0 = sources[i].position.y;
        double x1 = sources[i + 1].position.x;
        double y1 = sources[i + 1].position.y;
//        std::cout << "Adding source at (" << x0 << ", " << y0 << ") to (" << x1 << ", " << y1 << ")\n";
        double length_squared = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
        if(length_squared > 0) {
            for (auto& node: nodes) {
                double x = node.x * node_dist[0];
                double y = node.y * node_dist[1];
                double t = std::max(0.0, std::min(1.0, ((x - x0) * (x1 - x0) +
                                                        (y - y0) * (y1 - y0)) /
                                                       length_squared));
                double proj_x = x0 + t * (x1 - x0);
                double proj_y = y0 + t * (y1 - y0);
                double dist = std::sqrt((x - proj_x) * (x - proj_x) + (y - proj_y) * (y - proj_y));
                if(dist < 0.05){
                    node.source = true;
//                    std::cout << "node at (" << x << ", " << y << ") is a source\n";
                }

            }
        }
    }
    initializeJ();
    initializePrior();

}

/**
 * Add a measurement v at a given position x, y at time t.
 * @param x measurement position x
 * @param y measurement position y
 * @param v measurement value
 * @param t time measurement was taken
 */
void GMRF::addMeasurement(double x, double y, double v, double t, measurement_type meas_t) {
    M_Assert(x <= params.size_x && y <= params.size_y, "Measurement location out of bounds!");
    num_meas += 1;
    meas_types.push_back(meas_t);
    auto shape_factors = getShapeFacs(x, y, true);

    f.conservativeResize(num_energies + num_meas, f.cols());
    y_meas.conservativeResize(num_energies + num_meas, y_meas.cols());
    y_meas(num_meas - 1 + num_energies, 0) = x;
    y_meas(num_meas - 1 + num_energies, 1) = y;
    y_meas(num_meas - 1 + num_energies, 2) = v;
    measurement_age.conservativeResize(num_meas, measurement_age.cols());
    measurement_age(num_meas -1, 0) = t; // TODO Parse correct time here in the future.
    information_matrix.conservativeResize(num_energies + num_meas, num_energies + num_meas);
    // TODO Make this for the right measurement type!
    information_matrix.insert(num_energies + num_meas - 1, num_energies + num_meas - 1) =
            1/( pow(params.sigma_s, 2) + pow(params.sigma_zeta, 2) * (t - measurement_age(num_meas - 1, 0)));

    // Update Jacobian for the 4 closest nodes
    int offset;
    switch(meas_t){
        case field:
            offset = 0;
            break;
        case flow_x:
            offset = num_nodes;
            break;
        case flow_y:
            offset = num_nodes*2;
            break;
    }
    J.conservativeResize(num_energies + num_meas, J.cols());
    int x_id = shape_factors[0];
    int y_id = shape_factors[1];
    int node_id = fieldEstimateID(x_id, y_id);
    nodes[node_id].num_meas += 1;
    if(nodes[node_id].source){
        near_source = true;
    }
    J.insert(num_energies + num_meas - 1, params.N_X * y_id + x_id + offset) = shape_factors[2 + 2];

    x_id = shape_factors[0] + 1;
    node_id = fieldEstimateID(x_id, y_id);
    nodes[node_id].num_meas += 1;
    if(nodes[node_id].source){
        near_source = true;
    }
    J.insert(num_energies + num_meas - 1, params.N_X * y_id + x_id + offset) = shape_factors[2 + 3];

    y_id = shape_factors[1] + 1;
    node_id = fieldEstimateID(x_id, y_id);
    nodes[node_id].num_meas += 1;
    if(nodes[node_id].source){
        near_source = true;
    }
    J.insert(num_energies + num_meas - 1, params.N_X * y_id + x_id + offset) = shape_factors[2 + 1];
    x_id = shape_factors[0];
    node_id = fieldEstimateID(x_id, y_id);
    nodes[node_id].num_meas += 1;
    if(nodes[node_id].source){
        near_source = true;
    }
    J.insert(num_energies + num_meas - 1, params.N_X * y_id + x_id + offset) = shape_factors[2 + 0];

//    updateParameters(t);
    //updateEstimates(t);
}

// Update Procedure


/**
 * Updates the actual estimates + uncertainties
 * @param t current time
 */
void GMRF::updateEstimates(double t) {
    auto time = std::chrono::high_resolution_clock::now();
    updateInfomat(t);
    updateF();
    bool converged = false;
    Matrix<double, Dynamic, 1> x1;
    double step_size = 1.0;
    double g0 = (f - y_meas.col(2)).transpose() * information_matrix*(f - y_meas.col(2));
    double g = g0;
    bool first_run = true;
    Matrix<double, Dynamic, 1> best_sol;
    while(!converged){
        if(g < g0 || first_run) {
            step_size = 1.0;
            for(auto &e: energy_funs){
                e->updateJ();
            }
            best_sol = all_estimates;
//            std::cout << "Successful step, Cost Fun: " << g << ", step size: " << step_size << ", x1 norm: " << x1.norm() << ", relative change: " << abs((g - g0)/g0) << std::endl;
//            MatrixXf temp;
//            temp.resize(num_nodes*3, num_nodes*3);
//            temp.setZero();
//            temp.setIdentity();
//            MatrixXd tempo = -J.transpose() * information_matrix * (f -
//                    y_meas.col(2));
//            std::cout << tempo.maxCoeff() << std::endl;
//            std::cout << tempo.minCoeff() << std::endl;
            x1 = solver.compute(J.transpose() * information_matrix * J ).solve(
                    -J.transpose() * information_matrix * (f - y_meas.col(2)));
            if((x1.norm() < 1e-5|| (abs((g - g0)/g0) < 1e-5) && !first_run)){
//                std::cout << "Converged..\n";
                converged = true;
                break;
            }
            first_run = false;
            g0 = g;
            if (solver.info() != Eigen::Success) {
                MatrixXd temp = J.transpose() * information_matrix * J;
//                std::cout << temp.eigenvalues() << std::endl;
//                std::cout << solver.info() << std::endl;
//                std::cout << temp << std::endl;
//                std::cout << temp.inverse() << std::endl;
//                std::cout << (temp * 25).inverse() * 25 << std::endl;
                throw std::runtime_error("Eigen Solver failed!");
            }
            all_estimates += step_size * x1;
        }
        else{
            step_size /= 2.0;
//            std::cout << "Unsuccessful step, Cost Fun: " << g << ", new step size: " << step_size << ", x1 norm: " << x1.norm() << std::endl;
            all_estimates -= step_size * x1;
            if(step_size < 0.0001){// || (g - g0)/g0 < 1e-5){
                converged = true;
            }
        }
        updateF();
        g = (f - y_meas.col(2)).transpose() * information_matrix*(f - y_meas.col(2));
    }
//    std::cout << "Converged" << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
//    std::cout << "Time to update estimates: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - time).count() << std::endl;
    all_estimates = best_sol;
    time = std::chrono::high_resolution_clock::now();
    MatrixXd upperTriang = solver.matrixU();
    MatrixXi indices = solver.permutationP().indices();
    Sigma.setZero();
    // Get the covariances
    for (int l=num_nodes*3-1; l>=0; l--)
    {
        //Computes variances in the inferior submatrix of "l"
        double subSigmas = 0.0;
        for(int j=l+1; j<num_nodes*3; j++)
        {
            if (upperTriang(l,j) != 0)
            {
                //Compute off-diagonal variances Sigma(j,l) = Sigma(l,j);

                //SUM 1
                double sum = 0.0;
                for(int i=l+1; i<=j; i++)
                {
                    if( upperTriang(l,i) !=0 )
                    {
                        sum += upperTriang(l,i) * Sigma(i,j);
                    }
                }
                //SUM 2
                for(int i=j+1; i<num_nodes*3; ++i)
                {
                    if( upperTriang(l,i) !=0 )
                    {
                        sum += upperTriang(l,i) * Sigma(j,i);
                    }
                }
                //Save off-diagonal variance (only Upper triangular)
                Sigma(l,j) = ( -sum/upperTriang(l,l) );
                subSigmas += upperTriang(l,j) * Sigma(l,j);
            }
        }

        Sigma(l,l) = (1/upperTriang(l,l)) * ( 1/upperTriang(l,l) - subSigmas );
    }
    stop = std::chrono::high_resolution_clock::now();
    int counter = 0;
    for (int i = 0; i < params.N_Y; i++){
        for(int j = 0; j< params.N_X; j++) {
            // Recover the diagonal covariance values, undoing the permutation:
            int idx = indices.coeff(counter);
            covariances[field](j, i) = Sigma.coeff(idx, idx);
            counter += 1;
        }
    }
    for (int i = 0; i < params.N_Y; i++){
        for(int j = 0; j< params.N_X; j++) {
            // Recover the diagonal covariance values, undoing the permutation:
            int idx = indices.coeff(counter);
            covariances[flow_x](j, i) = Sigma.coeff(idx, idx);
            counter += 1;
        }
    }
    for (int i = 0; i < params.N_Y; i++){
        for(int j = 0; j< params.N_X; j++) {
            // Recover the diagonal covariance values, undoing the permutation:
            int idx = indices.coeff(counter);
            covariances[flow_y](j, i) = Sigma.coeff(idx, idx);
            counter += 1;
        }
    }
//    std::cout << "Time to update covariances: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - time).count() << std::endl;
    updateViz(t);
//    std::cout << "Max est: " << estimates.at(field).maxCoeff() << std::endl;
}

void GMRF::updateInfomat(double t) {
    updateMeasurementInfomat(t);
}

/**
 * Update infomatrix based on edge id.
 * @param L_id edge id
 * @param prob prob
 */
void GMRF::updateCorrelationInfomat(int L_id, double prob) {
//    std::cout << "Updated edge at (" << edges[L_id].x0 << ", " << edges[L_id].y0 << "), (" << edges[L_id].x1 << ", " << edges[L_id].y1 << ") to " << prob << std::endl;
    information_matrix.coeffRef(L_id, L_id) = pow((1 - prob)/sigma_p, 2);
    information_matrix.coeffRef(L_id + num_edges, L_id + num_edges) = pow((1
                                                                           - prob)/ff_sigma_p, 2);
    information_matrix.coeffRef(L_id + num_edges*2, L_id + num_edges*2) = pow((1
                                                                               - prob)/ff_sigma_p, 2);

    if(prob == 1.0){
        edges[L_id].is_cut = true;
    }
    if(viz) {
        std::cout << "VIZZZZ" << std::endl;
        lines[L_id * 2].color = sf::Color(255, 255, 255, 255 * (1 - prob));
        lines[L_id * 2 + 1].color = sf::Color(255, 255, 255, 255 * (1 - prob));
    }
}

void GMRF::updateMeasurementInfomat(double t) {
    // starting after the prior which consists of num_nodes*3 measurements
    for(int i = num_nodes*3; i < num_meas; i++){
        if(meas_types[i] == field) {
            information_matrix.coeffRef(num_energies + i, num_energies + i) =
                    1 / (pow(params.sigma_s, 2) + pow(params.sigma_zeta, 2) * (t - measurement_age(i, 0)));
        }
        else{
            information_matrix.coeffRef(num_energies + i, num_energies + i) =
                    1 / (pow(params.ff_sigma_s, 2) + pow(params.ff_sigma_zeta, 2) * (t - measurement_age(i, 0)));
        }
    }
}


void GMRF::update_J() {

}


void GMRF::updateF() {
    // Field neighborhood
    f.setZero();
    for(auto& e : energy_funs){
        e->updateF();
    }
    // From here on observations
    // Prior
    int counter = 0;
    for(int j = 0; j < params.N_Y; j++){
        for(int i = 0; i < params.N_X; i++){
            f(counter + num_energies, 0) = estimates.at(field)(xyToN(i, j));
            counter += 1;
            f(counter + num_energies, 0) = estimates.at(flow_x)(xyToN(i, j));
            counter += 1;
            f(counter + num_energies, 0) = estimates.at(flow_y)(xyToN(i, j));
            counter += 1;
        }
    }

    // True Measurements
    for(int i = counter; i < num_meas; i++){
        f(i + num_energies, 0) = 0;
        auto shape_factors = getShapeFacs(y_meas(i + num_energies, 0), y_meas(i + num_energies, 1), true);
        int x_id = shape_factors[0];
        int y_id = shape_factors[1];
        f(i + num_energies, 0) += shape_factors[2 + 2] * estimates.at(meas_types[i])(xyToN(x_id, y_id));
        x_id = shape_factors[0] + 1;
        f(i + num_energies, 0) += shape_factors[2 + 3] * estimates.at(meas_types[i])(xyToN(x_id, y_id));
        y_id = shape_factors[1] + 1;
        f(i + num_energies, 0) += shape_factors[2 + 1] * estimates.at(meas_types[i])(xyToN(x_id, y_id));
        x_id = shape_factors[0];
        f(i + num_energies, 0) += shape_factors[2 + 0] * estimates.at(meas_types[i])(xyToN(x_id, y_id));
    }
}

// Initialization Procedure

/**
 * Add prior as observations fully defined by prior mean and prior variance.
 * Each node gets one observation.
 */
void GMRF::initializePrior(){
    num_meas += 3*params.N_X*params.N_Y;  // we will have N_X*N_Y observations for
    // the estimate + 2 * N_X*N_Y observations for the flowfield

    // Resize the corresponding matrices and vectors
    f.conservativeResize(num_energies + num_meas, f.cols());
    y_meas.conservativeResize(num_energies + num_meas, y_meas.cols());
    measurement_age.conservativeResize(num_meas, measurement_age.cols());
    information_matrix.conservativeResize(num_energies + num_meas, num_energies + num_meas);
    J.conservativeResize(num_energies + num_meas, J.cols());
    double t = 0.0;
    // Fill the matrices + vector
    int counter = 0;
    for(int j = 0; j < params.N_Y; j++){
        for(int i = 0; i < params.N_X; i++) {
            meas_types.emplace_back(field);
            double x = (double) i * node_dist[0];
            double y = (double) j * node_dist[1];
            y_meas(counter + num_energies, 0) = x;
            y_meas(counter + num_energies, 1) = y;
            y_meas(counter + num_energies,
                   2) = params.prior_mean;  // Prior mean
            measurement_age(counter, 0) = 0;
            information_matrix.insert(num_energies + counter,
                                      num_energies + counter) =
                    1 / (pow(params.prior_sigma, 2));  // Prior variance
            J.insert(num_energies + counter, fieldEstimateID(i,
                                                             j)) = 1;  // No shape functions here.
            double s = 1.0;
            nodes[fieldEstimateID(i, j)].num_meas += 1;
            counter += 1;

            meas_types.emplace_back(flow_x);
            y_meas(counter + num_energies, 0) = x;
            y_meas(counter + num_energies, 1) = y;
            y_meas(counter + num_energies,
                   2) = params.ff_prior_mean_x;  // Prior mean
            measurement_age(counter, 0) = 0;
            information_matrix.insert(num_energies + counter,
                                      num_energies + counter) =
                    1 / (pow(params.ff_prior_sigma, 2));  // Prior variance
            J.insert(num_energies + counter,
                     flowEstimateID(i, j)) = 1;  // No shape functions here.
            counter += 1;

            meas_types.emplace_back(flow_y);
            y_meas(counter + num_energies, 0) = x;
            y_meas(counter + num_energies, 1) = y;
            y_meas(counter + num_energies,
                   2) = params.ff_prior_mean_y;  // Prior mean
            measurement_age(counter, 0) = 0;
            information_matrix.insert(num_energies + counter,
                                      num_energies + counter) =
                    1 / (pow(params.ff_prior_sigma, 2));  // Prior variance
            J.insert(num_energies + counter, flowEstimateID(i, j) +
                                             num_nodes) = 1;  // No shape functions here.
            counter += 1;
        }
    }
    updateEstimates(0.0);  // Generate initial estimates.
}

/**
 * Initialize information matrix content for prior
 */
void GMRF::initializeInformationMatrix() {
    for(auto &e: energy_funs){
        e->initializeInfomat();
    }
}

void GMRF::initializeNodes() {
    int counter = 0;

    for(int j = 0; j < params.N_Y; j++){
        for(int i = 0; i < params.N_X; i++){
            nodes[counter].dists = node_dist;
            nodes[counter].x = i;
            nodes[counter].y = j;
            nodes[counter].node_id = fieldEstimateID(i, j);
            counter += 1;
        }
    }
}

void GMRF::initializeEdges() {
    int counter = 0;
    for (int i = 0; i < params.N_X - 1; i++) {
        for(int j = 0; j < params.N_Y; j++) {
            edges[counter].horiz = true;
            edges[counter].N0 = fieldEstimateID(i, j);
            edges[counter].x0 = i;
            edges[counter].y0 = j;
            edges[counter].N1 = fieldEstimateID(i + 1, j);
            edges[counter].x1 = i + 1;
            edges[counter].y1 = j;
            counter += 1;
        }
    }
    for (int i = 0; i < params.N_X; i++) {
        for(int j = 0; j < params.N_Y - 1; j++) {
            edges[counter].horiz = false;
            edges[counter].N0 = fieldEstimateID(i, j);
            edges[counter].x0 = i;
            edges[counter].y0 = j;
            edges[counter].N1 = fieldEstimateID(i, j + 1);
            edges[counter].x1 = i;
            edges[counter].y1 = j + 1;
            counter += 1;
        }
    }
    for(auto& edge: edges){
        nodes[xyToN(edge.x0, edge.y0)].connected_edges.push_back(&edge);
        nodes[xyToN(edge.x1, edge.y1)].connected_edges.push_back(&edge);
    }

}
void GMRF:: initializeJ() {
    J.resize(num_energies, num_nodes*3);
    for(auto &e: energy_funs) {
        e->initializeJ();
    }
}

void GMRF::prepareViz() {
    viz = true;
    lines.resize(num_edges * 2);
    lines.setPrimitiveType(sf::Lines);
    int counter = 0;
    for(int i = 0; i < 100 * num_nodes; i++){
        particle_viz.push_back(new sf::CircleShape(0.1));
        particle_viz[i]->setFillColor(sf::Color::White);
        particle_viz[i]->setOrigin(0.1, 0.1);
    }
    for(int i = 0; i < num_nodes; i++){
        source_viz.push_back(new sf::CircleShape(0.1));
        source_viz[i]->setOrigin(0.1, 0.1);
        source_viz[i]->setFillColor(sf::Color::Blue);
        source_viz[i]->setPosition(nodes[i].x * node_dist[0], nodes[i].y * node_dist[1]);
    }
    for(int i = 0; i < num_edges; i++){
        lines[counter].position = sf::Vector2f(edges[i].x0 * node_dist[0], edges[i].y0 * node_dist[1]);
        lines[counter+1].position = sf::Vector2f(edges[i].x1 * node_dist[0], edges[i].y1 * node_dist[1]);
        lines[counter].color=sf::Color(255, 255, 255, 255 * (1 - edges[i]
                .is_cut));
        lines[counter+1].color=sf::Color(255, 255, 255, 255 * (1 - edges[i]
                .is_cut));
        counter += 2;
    }
}

// Helper Functions

/**
 * Get the 4 shape function weights and x, y id of point 2.
 * ordering:
 * 0 --- 1
 * 2 --- 3
 * @param dx
 * @param dy
 * @return
 */
std::array<double, 6> GMRF::getShapeFacs(double& x, double& y, bool normalized)
const{
    M_Assert(x <= params.size_x && y <= params.size_y, "Requested value out of bounds for the field!");
    // x ID of closest node (bottom left)
    double id_x_bottom_left = std::min(floor(x / node_dist[0]), (double)params
            .N_X -
                                                                2);
    // y ID of closest
    // node (bottom left)
    double id_y_bottom_left = std::min(floor(y / node_dist[1]), (double)params
            .N_Y -2);
    // Cell center point coordinate.
    double dx = x - (id_x_bottom_left * node_dist[0] + node_dist[0] / 2);
    double dy = y - (id_y_bottom_left * node_dist[1] + node_dist[1] / 2);
    bool bottom_left_used = nodes[xyToN(id_x_bottom_left, id_y_bottom_left)]
            .reachable;
    bool bottom_right_used = nodes[xyToN(id_x_bottom_left + 1, id_y_bottom_left)]
            .reachable;
    bool top_left_used = nodes[xyToN(id_x_bottom_left, id_y_bottom_left + 1)]
            .reachable;
    bool top_right_used = nodes[xyToN(id_x_bottom_left + 1, id_y_bottom_left + 1)]
            .reachable;
    auto top = edges[nodepairToL(id_x_bottom_left, id_x_bottom_left + 1,
                                 id_y_bottom_left + 1, id_y_bottom_left + 1)];
    auto bottom = edges[nodepairToL(id_x_bottom_left, id_x_bottom_left + 1,
                                    id_y_bottom_left, id_y_bottom_left)];
    auto left = edges[nodepairToL(id_x_bottom_left, id_x_bottom_left,
                                  id_y_bottom_left, id_y_bottom_left + 1)];
    auto right = edges[nodepairToL(id_x_bottom_left + 1, id_x_bottom_left + 1,
                                   id_y_bottom_left, id_y_bottom_left + 1)];
    if(left.is_cut){
        if(y < left.max_cut_pos && x < left.max_cut_max_depth) {
            top_left_used = false;
            top_right_used = false;
        }
        if(y > left.min_cut_pos && x < left.min_cut_max_depth){
            bottom_left_used = false;
            bottom_right_used = false;
        }
    }
    if(right.is_cut){
        if(y < right.max_cut_pos && x > right.max_cut_min_depth){
            top_right_used = false;
            top_left_used = false;
        }
        if( y > right.min_cut_pos && x > right.min_cut_min_depth){
            bottom_right_used = false;
            bottom_left_used = false;
        }
    }
    if(top.is_cut){
        if(x < top.max_cut_pos && y > top.max_cut_min_depth) {
            top_right_used = false;
            bottom_right_used = false;
        }
        if (x > top.min_cut_pos && y > top.min_cut_min_depth){
            top_left_used = false;
            bottom_left_used = false;
        }
    }
    if(bottom.is_cut){
        if(x < bottom.max_cut_pos && y < bottom.max_cut_max_depth) {
            bottom_right_used = false;
            top_right_used = false;
        }
        if(x > bottom.min_cut_pos && y < bottom.min_cut_max_depth){
            bottom_left_used = false;
            top_left_used = false;
        }
    }
    double val_top_left = -1 / (node_dist[0] * node_dist[1]) * (dx - node_dist[0] / 2) *
                          (dy + node_dist[1] / 2) * top_left_used;
    double val_top_right = 1 / (node_dist[0] * node_dist[1]) * (dx + node_dist[0] / 2) *
                           (dy + node_dist[1] / 2) * top_right_used;
    double val_bottom_left = 1 / (node_dist[0] * node_dist[1]) * (dx - node_dist[0] / 2) *
                             (dy - node_dist[1] / 2) * bottom_left_used;
    double val_bottom_right = -1 / (node_dist[0] * node_dist[1]) * (dx + node_dist[0] / 2) *
                              (dy - node_dist[1] / 2) * bottom_right_used;
    double normalization_factor = 1;
    if(normalized or true){
        normalization_factor = val_top_left + val_top_right + val_bottom_left + val_bottom_right;
        if(normalization_factor <= 0){
            normalization_factor = 1;
        }
    }
    return std::array<double, 6>({
                                         id_x_bottom_left,
                                         id_y_bottom_left,
                                         val_top_left/normalization_factor,
                                         val_top_right/normalization_factor,
                                         val_bottom_left/normalization_factor,
                                         val_bottom_right/normalization_factor
                                 });
}


std::array<double, 2> GMRF::uncertaintyGradient(double x, double y, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& field) const {
    M_Assert(x <= params.size_x && y <= params.size_y && x >= 0 && y >= 0,
             "Requested value out of bounds for the field!");
    // TODO Handle out of bounds requests
    if(x < 0){
        return std::array<double, 2>({1, 0});
    }
    if(x > params.size_x){
        return std::array<double, 2>({-1, 0});
    }
    if(y < 0){
        return std::array<double, 2>({0, 1});
    }
    if(y > params.size_y){
        return std::array<double, 2>({0, -1});
    }

    double id_x_bottom_left = floor(x / node_dist[0]); // x ID of closest node (bottom left)
    double id_y_bottom_left = floor(y / node_dist[1]); // y ID of closest node (bottom left)
    double m0 = field.coeff(id_x_bottom_left, id_y_bottom_left + 1);
    double m1 = field.coeff(id_x_bottom_left + 1, id_y_bottom_left + 1);
    double m2 = field.coeff(id_x_bottom_left, id_y_bottom_left);
    double m3 = field.coeff(id_x_bottom_left + 1, id_y_bottom_left);
    double lx = node_dist[0];
    double ly = node_dist[1];
    double x0 = id_x_bottom_left * lx;
    double y0 = id_y_bottom_left * ly;
    double gradient_x = 0.5*(-m0*(ly*y0 - y) + m1*(ly*y0 - y) + m2*(ly*y0 +
                                                                    ly - y) - m3*(ly*y0 + ly - y))/(m0*(ly*y0 - y)*(lx*x0 + lx - x) - m1*(lx*x0 - x)*(ly*y0 - y) - m2*(lx*x0 + lx - x)*(ly*y0 + ly - y) + m3*(lx*x0 - x)*(ly*y0 + ly - y));
    double gradient_y = 0.5*(-m0*(lx*x0 + lx - x) + m1*(lx*x0 - x) + m2*
                                                                     (lx*x0 + lx - x) - m3*(lx*x0 - x))/(m0*(ly*y0 - y)*(lx*x0 + lx - x) - m1*(lx*x0 - x)*(ly*y0 - y) - m2*(lx*x0 + lx - x)*(ly*y0 + ly - y) + m3*(lx*x0 - x)*(ly*y0 + ly - y));
    return std::array<double, 2>({gradient_x, gradient_y});
}

/**
 * returns the edge id in L of a pair of neighboring nodes
 * @param x0 x position node 0
 * @param x1 x position node 1
 * @param y0 y position node 0
 * @param y1 y position node 1
 * @return
 */
int GMRF::nodepairToL(int x0, int x1, int y0, int y1) const{
    if(x0 == x1){
        if(y1 == y0){
            throw std::runtime_error("Supplied the same node twice, no edge can be determined!");
        }
        // vertical edge, comes after all horizontal edge.
        if(y1 > y0){
            // then, node 0 is the "root" node of the edge
            return (params.N_X-1)*params.N_Y + (params.N_Y - 1) * x0 + y0;
        }
        else{
            // then, node 1 is the "root" node of the edge
            return (params.N_X-1)*params.N_Y + (params.N_Y - 1) * x1 + y1;
        }
    }
    else if(y0 == y1){
        if(x1 == x0){
            throw std::runtime_error("Supplied the same node twice, no edge can be determined!");
        }
        // horizontal edge
        if(x1 > x0){
            return params.N_Y * x0 + y0;
            // then, node 0 is the "root" node of the edge
        }
        else{
            return params.N_Y * x1 + y1;
            // then, node 1 is the "root" node of the edge
        }
    }
    else{
        throw std::runtime_error("The two nodes provided are not neighboring nodes, got node 0 located at (" +
                                 std::to_string(x0) + ", " + std::to_string(y0) + ") and node 1 located at (" +
                                 std::to_string(x1) + ", " + std::to_string(y1) + ")!");
    }
}

double GMRF::getEstimate(double x, double y) const {
    return std::min(std::max(getEstimate(x, y, field), 0.0), 1.0);
}

double GMRF::getEstimate(double x, double y, measurement_type m_type) const {
    VecToMat estimates_mapped(estimates.at(m_type).data(), params.N_X, params
            .N_Y);
    double ret = 0;
    auto shape_factors = getShapeFacs(x, y, true);
    int x_id = shape_factors[0];
    int y_id = shape_factors[1];
    ret += shape_factors[2 + 2] * estimates_mapped(x_id, y_id);
    x_id = shape_factors[0] + 1;
    ret += shape_factors[2 + 3] * estimates_mapped(x_id, y_id);
    y_id = shape_factors[1] + 1;
    ret += shape_factors[2 + 1] * estimates_mapped(x_id, y_id);
    x_id = shape_factors[0];
    ret += shape_factors[2 + 0] * estimates_mapped(x_id, y_id);
    return ret;
}


double GMRF::getCovariance(double x, double y) const{
    return getCovariance(x, y, field);
}


double GMRF::getCovariance(double x, double y, measurement_type m_type) const {
    double ret = 0;
    auto shape_factors = getShapeFacs(x, y, true);
    int x_id = shape_factors[0];
    int y_id = shape_factors[1];
    ret += shape_factors[2 + 2] * covariances.at(m_type)(x_id, y_id);
    x_id = shape_factors[0] + 1;
    ret += shape_factors[2 + 3] * covariances.at(m_type)(x_id, y_id);
    y_id = shape_factors[1] + 1;
    ret += shape_factors[2 + 1] * covariances.at(m_type)(x_id, y_id);
    x_id = shape_factors[0];
    ret += shape_factors[2 + 0] * covariances.at(m_type)(x_id, y_id);
    return ret;
}


void GMRF::addToViz(RealtimeViz& viz_obj, double min_est, double max_est, double min_std, double max_std) {
    this->viz_obj = &viz_obj;
    viz = true;
    this->min_est = min_est;
    this->max_est = max_est;
    this->min_cov = min_std;
    this->max_cov = max_std;
    viz_obj.addVizObject(this, map_types::std_map);
}


void GMRF::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    if(!viz){
        throw std::runtime_error("GMRF not configured for visualization. This shouldn't ever happen!");
    }
    target.draw(lines, states);
    for(auto c: source_viz){
        target.draw(*c, states);
    }
    target.draw(source, states);

    for(auto part:particle_viz){
        target.draw(*part, states);
    }
}


const EstToOneMap * GMRF::getFlowXMap() const {
    return &estimates.at(flow_x);
}

const EstToOneMap * GMRF::getFlowYMap() const {
    return &estimates.at(flow_y);
}

double GMRF::interpolateValue(double x, double y, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &field) const {
    auto shape_factors = getShapeFacs(x, y, true);
    double ret = 0;
    int x_id = shape_factors[0];
    int y_id = shape_factors[1];
    ret += shape_factors[2 + 2] * field(x_id, y_id);
    x_id = shape_factors[0] + 1;
    ret += shape_factors[2 + 3] * field(x_id, y_id);
    y_id = shape_factors[1] + 1;
    ret += shape_factors[2 + 1] * field(x_id, y_id);
    x_id = shape_factors[0];
    ret += shape_factors[2 + 0] * field(x_id, y_id);
    return ret;
}

VecToMat GMRF::getEstimates() const {
    return {estimates.at(field).data(), params.N_X, params.N_Y};
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&
GMRF::getCovariances() const {
    return covariances.at(field);
}

void GMRF::sourceHypothesisNeighbor(double &x_s, double &y_s, bool &success) {
    int num_parts = num_nodes; // one particle per node
    VecToMat estimates_mapped(estimates.at(field).data(), params.N_X, params.N_Y);
    VecToMat vel_x_mapped(estimates.at(flow_x).data(), params.N_X, params.N_Y);
    VecToMat vel_y_mapped(estimates.at(flow_y).data(), params.N_X, params.N_Y);

    double sum_ests = ((estimates_mapped.array() > 0)
                               .cast<double>()*estimates_mapped.array()).matrix().sum();
    if(sum_ests == 0) {
//        std::cout << "Need more measurements" << "\n";
        success = false;
        return;
    }

    // initialize the weights for each particle
    Eigen::Matrix<double, Eigen::Dynamic, 1> weights;
    Eigen::Matrix<double, Eigen::Dynamic, 1> weights_old;
    weights.resize(num_parts, 1);
    for(int i = 0; i < num_nodes; i++){
        int x = i % params.N_X;
        int y = i / params.N_X;

        weights(i, 0) = 1 - 0.5 * (1 + erf((0.0 - estimates_mapped(x, y)) /
                                           sqrt(2* covariances.at(field)(x, y))));
        if(weights(i, 0) > 1){
            throw "xd";
        }
    }
    Eigen::Matrix<double, Eigen::Dynamic, 1> weights_0 = weights;
    weights.setZero();
    weights_old = weights;
    // discrete simulation for each particle
    static thread_local std::mt19937 generator(std::random_device{}());
    std::normal_distribution<double> distribution(0.0,1.0);
    // Get transition probabilities for each node.
    Eigen::Matrix<double, Eigen::Dynamic, 8> transition_probs; // 0: right,
    // 1: up, 2: left, 3: down
    transition_probs.resize(num_nodes, 8);
    transition_probs.setZero();
    int num_parts_node = 200;
    for(int i = 0; i < num_nodes; i++){
        int x = i % params.N_X;
        int y = i / params.N_X;
        for(int p = 0; p < num_parts_node; p++){
            double std_x = distribution(generator)*0.1;
            double std_y = distribution(generator)*0.1;
            double x_ = vel_x_mapped(x, y) + sqrt(std::max(covariances.at(flow_x)
                                                                   .operator()
                                                                           (x,
                                                                            y), 0.0)) *
                                             std_x;
            double y_ = vel_y_mapped(x, y) + sqrt(std::max(covariances.at(flow_y)
                                                                   .operator()
                                                                           (x,
                                                                            y), 0.0)) *
                                             std_y;
            double angle = atan2(y_, x_);
            int dx = 0;
            int dy = 0;
            // Determine which of the 8 neighbors is the most likely
            // TODO Adjust angle to cell size
            // TODO Maybe consider 8 neighbors instead of the 4
            if(params.all_neighbors) {
                int neighbor_id = round((angle + M_PI)/ (M_PI / 4));
                if(neighbor_id == 8){
                    neighbor_id = 0;
                }
                transition_probs(i, neighbor_id) += 1.0;
            }
            else {
                if (angle > M_PI / 4 && angle < M_PI - M_PI / 4) {
                    dy = 1;
                } else if (angle < -M_PI / 4 && angle > -M_PI + M_PI / 4) {
                    dy = -1;
                } else if (angle < M_PI / 4 && angle > -M_PI / 4) {
                    dx = 1;
                } else if (angle > 3 * M_PI / 4 ||
                           angle < -M_PI / 2 - M_PI / 4) {
                    dx = -1;
                }
                if (dx != 0) {
                    if (x + dx >= 0 && x + dx < params.N_X) {
                        if (!edges[nodepairToL(x, x + dx, y, y)].is_cut) {
                            if (dx < 0) {
                                transition_probs(i, 2) += 1.0;
                            } else {
                                transition_probs(i, 0) += 1.0;
                            }
                        }
                    }
                }
                if (dy != 0) {
                    if (y + dy >= 0 && y + dy < params.N_Y) {
                        if (!edges[nodepairToL(x, x, y, y + dy)].is_cut) {
                            if (dy < 0) {
                                transition_probs(i, 3) += 1.0;
                            } else {
                                transition_probs(i, 1) += 1.0;
                            }

                        }
                    }
                }
            }
        }
        if(params.all_neighbors) {
            if(!(x - 1 >= 0 && !edges[nodepairToL(x, x - 1, y, y)].is_cut)){
                transition_probs(i, 0) = 0;
            }
            if(!(y - 1 >= 0 && x - 1 >= 0 && !edges[nodepairToL(x, x - 1, y, y)].is_cut && !edges[nodepairToL(x, x, y, y - 1)].is_cut)){
                transition_probs(i, 1) = 0;
            }
            if(!(y - 1 >= 0 && !edges[nodepairToL(x, x, y, y - 1)].is_cut)){
                transition_probs(i, 2) = 0;
            }
            if(!(y - 1 >= 0 && x + 1 < params.N_X && !edges[nodepairToL(x, x + 1, y, y)].is_cut && !edges[nodepairToL(x, x, y, y - 1)].is_cut)){
                transition_probs(i, 3) = 0;
            }
            if(!(x + 1 < params.N_X && !edges[nodepairToL(x, x + 1, y, y)].is_cut)){
                transition_probs(i, 4) = 0;
            }
            if(!(y + 1 < params.N_Y && x + 1 < params.N_X && !edges[nodepairToL(x, x + 1, y, y)].is_cut && !edges[nodepairToL(x, x, y, y + 1)].is_cut)){
                transition_probs(i, 5) = 0;
            }
            if(!(y + 1 < params.N_Y  && !edges[nodepairToL(x, x, y, y + 1)].is_cut)){
                transition_probs(i, 6) = 0;
            }
            if(!(y + 1 < params.N_Y && x - 1 >= 0 && !edges[nodepairToL(x, x - 1, y, y)].is_cut && !edges[nodepairToL(x, x, y, y + 1)].is_cut)){
                transition_probs(i, 7) = 0;
            }
        }
        transition_probs.row(i) /= num_parts_node;
        // Get the CDF for > 0.0
        double max_est_ = std::min(estimates_mapped.maxCoeff(), 1.0);
        double cdf = 1 - 0.5 * (1 + erf((0.0  - std::max(std::min(estimates_mapped(x, y), 1.0), 0.0)/max_est_) /
                                        sqrt(covariances.at(field)(x, y))));
        transition_probs.row(i) *= cdf;
    }

    while(true){
        for (int i = 0; i < num_nodes; i++) {
            double weight_val = 0;
            int x = i % params.N_X;
            int y = i / params.N_X;
            double reward_val = 0;
            if(params.all_neighbors){
                if(transition_probs(i, 0) != 0){
                    weight_val += transition_probs(i, 0) * (1 +
                                                            weights_old[xyToN(x - 1, y)]);
                }
                if(transition_probs(i, 1) != 0){
                    weight_val += transition_probs(i, 1) * (1 +
                                                            weights_old[xyToN(x - 1, y - 1)]);
                }
                if(transition_probs(i, 2) != 0){
                    weight_val += transition_probs(i, 2) * (1 +
                                                            weights_old[xyToN(x, y - 1)]);
                }
                if(transition_probs(i, 3) != 0){
                    weight_val += transition_probs(i, 3) * (1 +
                                                            weights_old[xyToN(x + 1, y - 1)]);
                }
                if(transition_probs(i, 4) != 0){
                    weight_val += transition_probs(i, 4) * (1 +
                                                            weights_old[xyToN(x + 1, y)]);
                }
                if(transition_probs(i, 5) != 0){
                    weight_val += transition_probs(i, 5) * (1 +
                                                            weights_old[xyToN(x + 1, y + 1)]);
                }
                if(transition_probs(i, 6) != 0){
                    weight_val += transition_probs(i, 6) * (1 +
                                                            weights_old[xyToN(x, y + 1)]);
                }
                if(transition_probs(i, 7) != 0){
                    weight_val += transition_probs(i, 7) * (1 +
                                                            weights_old[xyToN(x - 1, y + 1)]);
                }
            }
            else {
                if (transition_probs(i, 0) != 0) {
                    weight_val += transition_probs(i, 0)
                                  * (1 + weights_old(xyToN(x + 1, y), 0));
                    reward_val += transition_probs(i, 0)
                                  * 1;
                }
                if (transition_probs(i, 1) != 0) {
                    weight_val += transition_probs(i, 1)
                                  * (1 + weights_old(xyToN(x, y + 1), 0));
                    reward_val += 1;
                }
                if (transition_probs(i, 2) != 0) {
                    weight_val += transition_probs(i, 2)
                                  * (1 + weights_old(xyToN(x - 1, y), 0));
                    reward_val += 1;
                }
                if (transition_probs(i, 3) != 0) {
                    weight_val += transition_probs(i, 3)
                                  * (1 + weights_old(xyToN(x, y - 1), 0));
                    reward_val += 1;
                }
            }
            weights(i, 0) = std::min(weight_val, 100.0);
        }

        double delta_change = (weights - weights_old).norm()/weights.norm();
        if(delta_change < 0.001){
            break;
        }
        weights_old = weights;
    }
    int id;
    for(int i = 0; i < num_nodes; i++){
        if(weights(i, 0) < 0){
            weights(i, 0) = 0;
        }
    }
    double max_ = weights.maxCoeff(&id);
    weights /=  max_;
    source.setPosition(id % params.N_X * node_dist[0], id / params.N_X *
                                                       node_dist[1]);
    for(int i = 0; i < num_nodes; i++){
        int x = i % params.N_X;
        int y = i / params.N_X;
        source_likelihoods(x, y) = weights(i, 0);
    }
    if(viz) {
        for (int i = 0; i < num_nodes; i++) {
            source_viz[i]->setScale(weights(i, 0) * 6, weights(i, 0) * 6);
        }
    }
    x_s = id % params.N_X * node_dist[0];
    y_s = id / params.N_X * node_dist[1];
    source_x = x_s;
    source_y = y_s;
    source_found = true;
    success = true;
}


void GMRF::addLineObstacle(double x0, double y0, double x1, double y1) {
    // First the horizontal cuts.
    if(y1 < y0){
        std::swap(y0, y1);
        std::swap(x0, x1);
    }
    double ymin = y0;
    double ymax = y1;
    if(y0 != y1) {
        int y_id_min = (int)std::ceil((ymin - 0.01)/node_dist[1]);
        int y_id_max = (int)std::floor((ymax + 0.01)/node_dist[1]);
        for (int i = y_id_min; i <= y_id_max; i++) {
            double x;
            if (x0 == x1) {
                x = x0;
            } else {
                x = x0 + (x1 - x0) / (y1 - y0) * (i * node_dist[1] - y0);
            }
            if (x0 == x1 && (x0 == 0 || x0 == params.size_x)){
                int node_id;
                if(x0 == 0) {
                    node_id = xyToN(0, i);
                }
                else {
                    node_id = xyToN(params.N_X - 1, i);
                }
                nodes[node_id].x_wall = true;
            }
            else{
                int x_id = (int) std::floor(x / node_dist[0]);
                if (x_id >= 0 && x_id + 1 < params.N_X) {
                    int id = nodepairToL(x_id, x_id + 1, i, i);
                    edges[id].min_cut_pos = std::min(edges[id].min_cut_pos, x);
                    edges[id].min_cut_max_depth = std::max(edges[id]
                                                                   .min_cut_max_depth, y1);
                    edges[id].min_cut_min_depth = std::min(edges[id]
                                                                   .min_cut_min_depth, y0);
                    edges[id].max_cut_pos = std::max(edges[id].max_cut_pos, x);
                    edges[id].max_cut_max_depth = std::max(edges[id].max_cut_max_depth, y1);
                    edges[id].max_cut_min_depth = std::min(edges[id]
                                                                   .max_cut_min_depth, y0);
                    updateCorrelationInfomat(id, 1.0);
                } else {
                    int node_id = xyToN(x_id, i);
                    nodes[node_id].x_wall = true;
                }
            }
        }
    }
    // Then the vertical cuts.
    if(x1 < x0){
        std::swap(y0, y1);
        std::swap(x0, x1);
    }
    double xmin = x0;
    double xmax = x1;
    if(x0 != x1) {
        int x_id_min = (int)std::ceil((xmin - 0.01)/node_dist[0]);
        int x_id_max = (int)std::floor((xmax + 0.01)/node_dist[0]);
        for (int i = x_id_min; i <= x_id_max; i++) {
            double y;
            if (y0 == y1) {
                y = y0;
            } else {
                y = y0 + (y1 - y0) / (x1 - x0) * (i * node_dist[0] - x0);
            }
            if (y0 == y1 && (y0 == 0 || y0 == params.size_y)){
                int node_id;
                if(y0 == 0) {
                    node_id = xyToN(i, 0);
                }
                else {
                    node_id = xyToN(i, params.N_Y - 1);
                }
                nodes[node_id].y_wall = true;
            }
            else {
                int y_id = (int) std::floor(y / node_dist[1]);
                if (y_id >= 0 && y_id + 1 < params.N_Y) {
                    int id = nodepairToL(i, i, y_id, y_id + 1);
                    edges[id].min_cut_pos = std::min(edges[id].min_cut_pos, y);
                    edges[id].min_cut_max_depth = std::max(edges[id]
                                                                   .min_cut_max_depth, x1);
                    edges[id].min_cut_min_depth = std::min(edges[id]
                                                                   .min_cut_min_depth, x0);

                    edges[id].max_cut_pos = std::max(edges[id].max_cut_pos, y);
                    edges[id].max_cut_max_depth = std::max(edges[id]
                                                                   .max_cut_max_depth, x1);
                    edges[id].max_cut_min_depth = std::min(edges[id]
                                                                   .max_cut_min_depth, x0);
                    updateCorrelationInfomat(id, 1.0);
                } else {
                    int node_id = xyToN(i, y_id);
                    nodes[node_id].y_wall = true;
                }
            }
        }
    }
}

const Eigen::Matrix<double, Dynamic, 1> &GMRF::getAllEstimates() const {
    return all_estimates;
}

double GMRF::getLikelihood(double x, double y, double v) const {
    double estimate = std::min(std::max(getEstimate(x, y), 0.0), 1.0);
    double cov = getCovariance(x, y);
    return 1 / sqrt(2 * M_PI * cov) * exp(-0.5 * (pow(v - estimate,2)/cov));
}

double GMRF::getLogLikelihood(double x, double y, double v) const {
    double estimate = std::min(std::max(getEstimate(x, y), 0.0), 1.0);
    double cov = getCovariance(x, y);
    return - log(sqrt(2 * M_PI * cov)) - 1 / (2 * cov) * pow(estimate - v, 2);
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &
GMRF::getSourceLikelihoods() const {
    return source_likelihoods;
}

void GMRF::updateViz(double t) {
    if (viz) {
        for (int i = 0; i < params.N_X; i++) {
            for (int j = 0; j < params.N_Y; j++) {
                f_stds_map(i, j) = sqrt(covariances[flow_x](i, j) +
                                        covariances[flow_y](i, j));
            }
        }
        for (int i = 0; i < MAP_X; i++) {
            for (int j = 0; j < MAP_Y; j++) {
                stds_map(i, j) = sqrt(
                        getCovariance(params.size_x / MAP_X * (double) i,
                                      params.size_y / MAP_Y * (double) j));
                estimates_map(i, j) = getEstimate(params.size_x / MAP_X * (double) i,
                                                  params.size_y / MAP_Y * (double) j);
                if (i == 0 && j == 0) {
                    min_est = estimates_map(i, j);
                    max_est = estimates_map(i, j);
                    min_cov = stds_map(i, j);
                    max_cov = stds_map(i, j);
                } else {
                    max_est = std::max(estimates_map(i, j), max_est);
                    min_est = std::min(estimates_map(i, j), min_est);
                    max_cov = std::max(stds_map(i, j), max_cov);
                    min_cov = std::min(stds_map(i, j), min_cov);
                }
            }
        }
        viz_obj->update();
        viz_obj->setTime(t);
    }
}

void GMRF::setEstimates(const Matrix<double, Eigen::Dynamic, 1> &estimates,
                        double t) {
    if(estimates.rows() == all_estimates.rows()) {
        all_estimates = estimates;
    }
    else {
        throw std::runtime_error("provided estimates do not fit GMRF config.");
    }
    updateViz(t);
}

void GMRF::determineReachableNodes() {
    std::set<GMRFNode*> checked_nodes;
    GMRFNode* start_node = &nodes[0];
    start_node->reachable = true;
    checked_nodes.insert(start_node);
    checkChildren(start_node, checked_nodes);
    for(int i = 0; i < edges.size(); i++){
        auto& edge = edges[i];
        if(!nodes[xyToN(edge.x0, edge.y0)].reachable && !nodes[xyToN(edge.x1,
                                                                     edge.y1)].reachable){
            edge.is_cut = true;
            updateCorrelationInfomat(i, 1.0);
        }
        if(edge.is_cut){
            updateCorrelationInfomat(i, 1.0);
        }
    }

}

void GMRF::checkChildren(GMRFNode *node, std::set<GMRFNode *> &checked_nodes) {
    for(Edge* edge : node->connected_edges){
        if(!edge->is_cut){
            GMRFNode* child;
            if(edge->x0 == node->x && edge->y0 == node->y){
                child = &nodes[xyToN(edge->x1, edge->y1)];
            }
            else{
                child = &nodes[xyToN(edge->x0, edge->y0)];
            }
            if(!checked_nodes.contains(child)) {
                child->reachable = true;
                checked_nodes.insert(child);
                checkChildren(child, checked_nodes);
            }
        }
    }
}

const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &
GMRF::getCovariances(measurement_type meas_type) const {
    return covariances.at(meas_type);
}

bool GMRF::isReachable(double x, double y) const {
    auto shape_facs = getShapeFacs(x, y);
    return shape_facs[2] != 0 || shape_facs[3] != 0 || shape_facs[4] != 0 ||
           shape_facs[5] != 0;
}

void GMRF::setCovariance(
        const Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &covariances,
        double t,  measurement_type meas_type) {
    this->covariances.at(meas_type) = covariances;

}

GMRF::~GMRF() {
    // energy functions
    for (auto& energy : energy_funs){
        delete energy;
    }
}
int xyToN2(int x, int y, int nx, int ny){
    return x * ny + y;
}

EnergyFunc::~EnergyFunc() {}

std::ostream &operator<<(std::ostream &os, const Edge & edge) {
    os << "Edge from (" << edge.x0 << "," << edge.y0 << ") to (" << edge.x1 << "," << edge.y1 << ") ";
    return os;
}

FieldN_EnergyFunc::FieldN_EnergyFunc(GMRF *gmrf, int idx0) {
    num_energies = gmrf->num_edges;
    this->gmrf = gmrf;
    idx[0] = idx0;
    idx[1] = idx0 + num_energies - 1;
}

void FieldN_EnergyFunc::initializeJ() {
    int counter = idx[0];
    for(int i = 0; i < gmrf->params.N_X - 1; i++){
        for(int j = 0; j < gmrf->params.N_Y; j++) {
            gmrf->J.insert(counter,  gmrf->fieldEstimateID(i, j)) = 1;
            gmrf->J.insert(counter, gmrf->fieldEstimateID(i+1, j)) = -1;
            counter += 1;
        }
    }
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for(int j = 0; j < gmrf->params.N_Y - 1; j++) {
            gmrf->J.insert(counter, gmrf->fieldEstimateID(i, j)) = 1;
            gmrf->J.insert(counter, gmrf->fieldEstimateID(i, j + 1)) = -1;
            counter += 1;
        }
    }
}

void FieldN_EnergyFunc::initializeInfomat() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->information_matrix.insert(i, i) = pow(1.0/gmrf->sigma_p, 2);
    }
}

void FieldN_EnergyFunc::updateJ() {

}

void FieldN_EnergyFunc::updateF() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->f(i, 0) = gmrf->estimates.at(field)(gmrf->edges[i - idx[0]].N0, 0) - gmrf->estimates.at(field)(gmrf->edges[i - idx[0]].N1, 0);
    }
}


FlowN_EnergyFunc::FlowN_EnergyFunc(GMRF *gmrf, int idx0) {
    num_energies = gmrf->num_edges * 2;
    this->gmrf = gmrf;
    idx[0] = idx0;
    idx[1] = idx0 + num_energies - 1;
}

void FlowN_EnergyFunc::initializeJ() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X - 1; i++) {
        for(int j = 0; j < gmrf->params.N_Y; j++) {
            gmrf->J.insert(counter,  gmrf->flowEstimateID(i,j)) = 1;
            gmrf->J.insert(counter, gmrf->flowEstimateID(i + 1,j)) = -1;
            counter += 1;
            gmrf->J.insert(counter,  gmrf->flowEstimateID(i, j) + gmrf->num_nodes)  = 1;
            gmrf->J.insert(counter, gmrf->flowEstimateID(i + 1,j) + gmrf->num_nodes) = -1;
            counter += 1;
        }
    }
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for(int j = 0; j < gmrf->params.N_Y - 1; j++) {
            gmrf->J.insert(counter,  gmrf->flowEstimateID(i,j)) = 1;
            gmrf->J.insert(counter, gmrf->flowEstimateID(i,j + 1)) = -1;
            counter += 1;
            gmrf->J.insert(counter,  gmrf->flowEstimateID(i,j) + gmrf->num_nodes)  = 1;
            gmrf->J.insert(counter, gmrf->flowEstimateID(i,j + 1) + gmrf->num_nodes) = -1;
            counter += 1;
        }
    }
}

void FlowN_EnergyFunc::initializeInfomat() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->information_matrix.insert(i, i) = pow(1./gmrf->ff_sigma_p, 2);
    }
}

void FlowN_EnergyFunc::updateJ() {

}

void FlowN_EnergyFunc::updateF() {
    int counter = idx[0];
    for(int i = idx[0]; i <= idx[1] - gmrf->num_edges; i++){
        gmrf->f(counter, 0) = gmrf->estimates.at(flow_x)(gmrf->edges[i - idx[0]].N0, 0) - gmrf->estimates.at(flow_x)(gmrf->edges[i - idx[0]].N1, 0);
        counter += 1;
        gmrf->f(counter, 0) = gmrf->estimates.at(flow_y)(gmrf->edges[i - idx[0]].N0, 0) - gmrf->estimates.at(flow_y)(gmrf->edges[i - idx[0]].N1, 0);
        counter += 1;
    }
}

FlowC_EnergyFunc::FlowC_EnergyFunc(GMRF *gmrf, int idx0) {
    num_energies = (gmrf->params.N_X - 1)* (gmrf->params.N_Y-1);
    this->gmrf = gmrf;
    idx[0] = idx0;
    idx[1] = idx0 + num_energies - 1;
}

void FlowC_EnergyFunc::initializeJ() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X- 1; i++) {
        for (int j = 0; j < gmrf->params.N_Y - 1; j++) {
            // For each node, we need the neighboring nodes in positive direction (x,y),(x+1,y),(x+1, y+1),(x, y+1)
            // and their values for the Jacobian. The corresponding entries of the Jacobian remain constant.
            // (x,y)
            gmrf->J.insert(counter, gmrf->flowEstimateID(i, j)) = 1; // w_x(x, y)
            gmrf->J.insert(counter, gmrf->flowEstimateID(i, j) + gmrf->num_nodes) = 1; // w_y(x, y)
            // (x+1, y)
            if(i + 1 < gmrf->params.N_X) {
                gmrf->J.insert(counter, gmrf->flowEstimateID(i + 1, j)) = -1;
                gmrf->J.insert(counter, gmrf->flowEstimateID(i + 1, j) + gmrf->num_nodes) = 1;
            }
            if(j + 1 < gmrf->params.N_Y){
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j + 1) ) = 1;
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j + 1) + gmrf->num_nodes) = -1;
            }
            if(i + 1 < gmrf->params.N_X && j + 1 < gmrf->params.N_Y) {
                gmrf->J.insert(counter, gmrf->flowEstimateID(i + 1, j + 1)) = -1;
                gmrf->J.insert(counter, gmrf->flowEstimateID(i + 1, j + 1) + gmrf->num_nodes) = -1;
            }
            counter += 1;
        }
    }
}

void FlowC_EnergyFunc::initializeInfomat() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->information_matrix.insert(i, i) = pow(1.0/gmrf->params.ff_sigma_c, 2);
    }
}

void FlowC_EnergyFunc::updateJ() {

}

void FlowC_EnergyFunc::updateF() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X - 1; i++) {
        for (int j = 0; j < gmrf->params.N_Y - 1; j++) {
            // For each node, we need the neighboring nodes in positive direction (x,y),(x+1,y),(x+1, y+1),(x, y+1)
            // and their values for the Jacobian. The corresponding entries of the Jacobian remain constant.
            // (x,y)
            gmrf->f(counter, 0) = gmrf->estimates.at(flow_x)(gmrf->xyToN(i, j)); // w_x(x, y)
            gmrf->f(counter, 0) += gmrf->estimates.at(flow_y)(gmrf->xyToN(i, j)); // w_y(x, y)
            // (x+1, y)
            if(i + 1 < gmrf->params.N_X) {
                gmrf->f(counter, 0) -= gmrf->estimates.at(flow_x)(gmrf->xyToN(i + 1, j)); // w_x(x + 1, y)
                gmrf->f(counter, 0) += gmrf->estimates.at(flow_y)(gmrf->xyToN(i + 1, j)); // w_y(x + 1, y)
            }
            if(j + 1 < gmrf->params.N_Y){
                gmrf->f(counter, 0) += gmrf->estimates.at(flow_x)(gmrf->xyToN(i, j + 1)); // w_x(x + 1, y)
                gmrf->f(counter, 0) -= gmrf->estimates.at(flow_y)(gmrf->xyToN(i, j + 1)); // w_y(x + 1, y)
            }
            if(i + 1 < gmrf->params.N_X && j + 1 < gmrf->params.N_Y) {
                gmrf->f(counter, 0) -= gmrf->estimates.at(flow_x)(gmrf->xyToN(i + 1, j + 1)); // w_x(x + 1, y + 1)
                gmrf->f(counter, 0) -= gmrf->estimates.at(flow_y)(gmrf->xyToN(i + 1, j + 1)); // w_y(x + 1, y + 1)
            }
            counter += 1;
        }
    }
}

FlowField_EnergyFunc::FlowField_EnergyFunc(GMRF *gmrf, int idx0) {
    num_energies = gmrf->num_nodes;
    this->gmrf = gmrf;
    idx[0] = idx0;
    idx[1] = idx0 + num_energies - 1;
}

void FlowField_EnergyFunc::initializeJ() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for (int j = 0; j < gmrf->params.N_Y; j++) {
            if(i + 1 < gmrf->params.N_X && !gmrf->edges[
                    gmrf->nodepairToL(i, i+1, j, j)].is_cut){
                // Then x gradient exists.
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j)) = gmrf->estimates.at(field)(gmrf->xyToN(i + 1, j))
                                                                      - gmrf->estimates.at(field)(gmrf->fieldEstimateID(i, j));
                gmrf->J.insert(counter, gmrf->fieldEstimateID(i, j)) = -gmrf->estimates.at(flow_x)(gmrf->xyToN(i, j));
                gmrf->J.insert(counter, gmrf->fieldEstimateID(i + 1, j)) = gmrf->estimates.at(flow_x)(gmrf->xyToN(i, j));
            }
            if(j + 1 < gmrf->params.N_Y && !gmrf->edges[
                    gmrf->nodepairToL(i, i, j, j+1)].is_cut){
                // Then y gradient exists.
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j) + gmrf->num_nodes) = gmrf->estimates.at(field)(
                        gmrf->xyToN(i, j + 1))
                                                                                        - gmrf->estimates.at(field)(gmrf->fieldEstimateID(i, j));
                gmrf->J.insert(counter, gmrf->fieldEstimateID(i, j + 1)) = gmrf->estimates.at(flow_y)(gmrf->xyToN(i, j));
                // Some entries in J were already created in the case that i+1 < N_X.
                if(i + 1 < gmrf->params.N_X && !gmrf->edges[
                        gmrf->nodepairToL(i, i+1, j, j)].is_cut){
                    gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i, j)) += -gmrf->estimates.at(flow_y)(
                            gmrf->xyToN(i, j));
                }
                    // otherwise just create them.
                else{
                    gmrf->J.insert(counter, gmrf->fieldEstimateID(i, j)) = -gmrf->estimates.at(flow_y)(gmrf->xyToN(i, j));
                }
            }
            counter += 1;
        }
    }
}

void FlowField_EnergyFunc::initializeInfomat() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->information_matrix.insert(i, i) = pow(1.0/gmrf->params.sigma_flow, 2);
    }
}

void FlowField_EnergyFunc::updateJ() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for (int j = 0; j < gmrf->params.N_Y; j++) {
            if(i + 1 < gmrf->params.N_X && !gmrf->edges[
                    gmrf->nodepairToL(i, i+1, j, j)].is_cut){
                // Then x gradient exists.
                gmrf->J.coeffRef(counter, gmrf->flowEstimateID(i, j)) = gmrf->estimates.at(field)(gmrf->xyToN(i + 1, j))
                                                                        - gmrf->estimates.at(field)(gmrf->fieldEstimateID(i, j));
                gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i, j)) = -gmrf->estimates.at(flow_x)(gmrf->xyToN(i, j));
                gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i + 1, j)) = gmrf->estimates.at(flow_x)(
                        gmrf->xyToN(i, j));
            }
            if(j + 1 < gmrf->params.N_Y && !gmrf->edges[
                    gmrf->nodepairToL(i, i, j, j+1)].is_cut){
                // Then y gradient exists.
                gmrf->J.coeffRef(counter, gmrf->flowEstimateID(i, j) + gmrf->num_nodes) = gmrf->estimates.at(field)(
                        gmrf->xyToN(i, j + 1))
                                                                                          - gmrf->estimates.at(field)(gmrf->fieldEstimateID(i, j));
                gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i, j + 1)) = gmrf->estimates.at(flow_y)(
                        gmrf->xyToN(i, j));
                if(i + 1 < gmrf->params.N_X && !gmrf->edges[
                        gmrf->nodepairToL(i, i+1, j, j)].is_cut){
                    gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i, j)) += -gmrf->estimates.at(flow_y)(
                            gmrf->xyToN(i, j));
                }
                else{
                    gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i, j)) = -gmrf->estimates.at(flow_y)(
                            gmrf->xyToN(i, j));
                }
            }
            counter += 1;
        }
    }
}

void FlowField_EnergyFunc::updateF() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for (int j = 0; j < gmrf->params.N_Y; j++) {
            gmrf->f(counter, 0) = 0.0;
            if(!gmrf -> edges[gmrf->nodepairToL(i, i + 1, j, j)].is_cut && !gmrf -> edges[gmrf->nodepairToL(i, i, j, j +1)].is_cut) {
                if (i + 1 < gmrf->params.N_X && !gmrf->edges[
                        gmrf->nodepairToL(i, i + 1, j, j)].is_cut) {
                    gmrf->f(counter, 0) +=
                            (gmrf->estimates.at(field)(gmrf->xyToN(i + 1, j))
                             - gmrf->estimates.at(field)(gmrf->xyToN(i, j)))
                            * gmrf->estimates.at(flow_x)(gmrf->xyToN(i, j));
                }
                if (j + 1 < gmrf->params.N_Y && !gmrf->edges[
                        gmrf->nodepairToL(i, i, j, j + 1)].is_cut) {
                    gmrf->f(counter, 0) +=
                            (gmrf->estimates.at(field)(gmrf->xyToN(i, j + 1))
                             - gmrf->estimates.at(field)(gmrf->xyToN(i, j)))
                            * gmrf->estimates.at(flow_y)(gmrf->xyToN(i, j));
                }
            }
            counter += 1;
        }
    }
}

GMRFNode::GMRFNode() {
    num_meas = 0;
}

ObstacleFlow_EnergyFunc::ObstacleFlow_EnergyFunc(GMRF *gmrf, int idx0) {
    num_energies = gmrf->num_nodes * 4;  // 4 energies per node, left right
    // up down
    this->gmrf = gmrf;
    idx[0] = idx0;
    idx[1] = idx0 + num_energies - 1;
}

void ObstacleFlow_EnergyFunc::initializeJ() {
    int counter = idx[0];
    double val_no_obs = pow(1/gmrf->sigma_p, 2);
    for (int i = 0; i< gmrf->params.N_X; i++){
        for (int j = 0; j < gmrf->params.N_Y; j++){
            if (i + 1 < gmrf->params.N_X){
                int L_id = gmrf->nodepairToL(i, i + 1, j, j);
                // Not outter wall in x ir
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j)) =
                        (gmrf->edges[L_id].is_cut);
            }
            else {
                // Outter wall in x dir
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].x_wall) {
                    gmrf->J.insert(counter, gmrf->flowEstimateID(i, j)) = 1.0;
                }
            }
            counter += 1;
            if(i > 0){
                int L_id = gmrf->nodepairToL(i, i - 1, j, j);
                // Not outter wall in x ir
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j)) =
                        (gmrf->edges[L_id].is_cut);
            }
            else{
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].x_wall) {
                    gmrf->J.insert(counter, gmrf->flowEstimateID(i, j)) = 1.0;
                }
            }
            counter += 1;
            if (j + 1 < gmrf->params.N_Y){
                // Not outter wall in y dir
                int L_id = gmrf->nodepairToL(i, i, j, j + 1);
                // Not outter wall in x ir
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j) +
                                        gmrf->num_nodes) =
                        (gmrf->edges[L_id].is_cut);
            }
            else {
                // Outter wall in y dir
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].y_wall) {
                    gmrf->J.insert(counter, gmrf->flowEstimateID(i, j) +
                                            gmrf->num_nodes) = 1.0;
                }
            }
            counter += 1;
            if (j > 0){
                int L_id = gmrf->nodepairToL(i, i, j, j - 1);
                // Not outter wall in x ir
                gmrf->J.insert(counter, gmrf->flowEstimateID(i, j) +
                                        gmrf->num_nodes) =
                        (gmrf->edges[L_id].is_cut);
            }
            else{
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].y_wall) {
                    gmrf->J.insert(counter, gmrf->flowEstimateID(i, j) +
                                            gmrf->num_nodes) = 1.0;
                }
            }
            counter += 1;
        }
    }
}

void ObstacleFlow_EnergyFunc::initializeInfomat() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->information_matrix.insert(i, i) = pow(1.0/gmrf->params
                .ff_sigma_o, 2);
    }
}

void ObstacleFlow_EnergyFunc::updateJ() {

}

void ObstacleFlow_EnergyFunc::updateF() {
    int counter = idx[0];
    double val_no_obs = pow(1/gmrf->sigma_p, 2);
    for (int i = 0; i< gmrf->params.N_X; i++){
        for (int j = 0; j < gmrf->params.N_Y; j++){
            if (i + 1 < gmrf->params.N_X){
                int L_id = gmrf->nodepairToL(i, i + 1, j, j);
                // Not outer wall in x dir
                gmrf->f(counter, 0) = gmrf->estimates.at(flow_x)(gmrf->xyToN
                        (i, j))
                                      * (gmrf->edges[L_id].is_cut);
            }
            else {
                // Outter wall in x dir
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].x_wall) {
                    gmrf->f(counter, 0) = gmrf->estimates.at(flow_x)
                            (gmrf->xyToN(i, j));
                }
            }
            counter += 1;
            if(i > 0){
                int L_id = gmrf->nodepairToL(i, i - 1, j, j);
                // Not outer wall in x ir
                gmrf->f(counter, 0) = gmrf->estimates.at(flow_x)(gmrf->xyToN
                        (i, j))
                                      * (gmrf->edges[L_id].is_cut);
            }
            else{
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].x_wall) {
                    gmrf->f(counter, 0) = gmrf->estimates.at(flow_x)
                            (gmrf->xyToN(i, j));
                }
            }
            counter += 1;
            if (j + 1 < gmrf->params.N_Y){
                // Not outer wall in y dir
                int L_id = gmrf->nodepairToL(i, i, j, j + 1);
                // Not outer wall in x ir
                gmrf->f(counter, 0) = gmrf->estimates.at(flow_y)(gmrf->xyToN
                        (i, j))
                                      * (gmrf->edges[L_id].is_cut);
            }
            else {
                // Outter wall in y dir
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].y_wall) {
                    gmrf->f(counter, 0) = gmrf->estimates.at(flow_y)(gmrf->xyToN
                            (i, j));
                }
            }
            counter += 1;
            if (j > 0){
                int L_id = gmrf->nodepairToL(i, i, j, j - 1);
                // Not outer wall in x ir
                gmrf->f(counter, 0) = gmrf->estimates.at(flow_y)(gmrf->xyToN
                        (i, j))
                                      * (gmrf->edges[L_id].is_cut);
            }
            else{
                int node_id = gmrf->xyToN(i, j);
                if(gmrf->nodes[node_id].y_wall) {
                    gmrf->f(counter, 0) = gmrf->estimates.at(flow_y)
                            (gmrf->xyToN(i, j));
                }
            }
            counter += 1;
        }
    }
}

FieldBoundary_EnergyFunc::FieldBoundary_EnergyFunc(GMRF *gmrf, int idx0) {
    num_energies = gmrf->num_nodes;
    this->gmrf = gmrf;
    idx[0] = idx0;
    idx[1] = idx0 + num_energies - 1;
}

void FieldBoundary_EnergyFunc::initializeJ() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for (int j = 0; j < gmrf->params.N_Y; j++) {
            gmrf->J.insert(counter, gmrf->fieldEstimateID(i, j)) = 1.0;
            counter += 1;
        }
    }
}

void FieldBoundary_EnergyFunc::initializeInfomat() {
    for(int i = idx[0]; i <= idx[1]; i++){
        gmrf->information_matrix.insert(i, i) = pow(1./0.1, 2);
    }
}

void FieldBoundary_EnergyFunc::updateJ() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for (int j = 0; j < gmrf->params.N_Y; j++) {
            gmrf->J.coeffRef(counter, gmrf->fieldEstimateID(i, j)) = 1.39766584353102e-13*exp(30* gmrf->estimates.at(field)(gmrf->xyToN(i, j))) - 1.49361205103592*exp(-30* gmrf->estimates.at(field)(gmrf->xyToN(i, j)));
            counter += 1;
        }
    }
}

void FieldBoundary_EnergyFunc::updateF() {
    int counter = idx[0];
    for (int i = 0; i < gmrf->params.N_X; i++) {
        for (int j = 0; j < gmrf->params.N_Y; j++) {
            gmrf->f(counter, 0) =  exp(-30 * (gmrf->estimates.at(field)(gmrf->xyToN(i, j)) + 0.1)) + exp(30 * (gmrf->estimates.at(field)(gmrf->xyToN(i, j)) - 1.1));
            counter += 1;
        }
    }
}
