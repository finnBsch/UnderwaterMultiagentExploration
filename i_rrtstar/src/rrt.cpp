//
// Created by finn on 6/6/23.
//
#include "rrt.h"
#include <iostream>
#include <random>
#include <fstream>

RRTObstacle::RRTObstacle(float x0, float y0, float x1, float y1): x0(x0), y0
        (y0),
                                                                  x1(x1), y1(y1), line(sf::Lines, 2) {
    line[0].position = sf::Vector2f(x0, y0);
    line[0].color = sf::Color(120, 120, 120);
    line[1].position = sf::Vector2f(x1, y1);
    line[1].color = sf::Color(120, 120, 120);
    if(this->x0 > this->x1){
        float temp = this->x0;
        this->x0 = this->x1;
        this->x1 = temp;
        temp = this->y0;
        this->y0 = this->y1;
        this->y1 = temp;
    }
    length = sqrtf(pow(x1 - x0, 2) + pow(y1 - y0, 2));
    angle = atan2f(y1 - y0, x1 - x0);
}

void RRTObstacle::draw(sf::RenderTarget &target, sf::RenderStates states)
const {
    target.draw(line, states);
}

double RRTObstacle::getDistance(float x, float y) const {

    // get obstacle dir
    double dx_o = x1 - x0;
    double dy_o = y1 - y0;

    double dx = x0 - x;
    double dy = y0 - y;
    double d_sq = dx_o * dx_o + dy_o * dy_o; // Lenght squared of the obstacle

    double t = std::max(0.0, std::min(1.0, (dx_o * dx + dy_o * dy)/d_sq));
    double px = x0 + t * dx_o;
    double py = y0 + t * dy_o;
    return sqrt((px - x) * (px - x) + (py - y) * (py - y));
}

float RRT::localCostNodePair(Node *node_parent, Node *node_child) const {
    return localCostNodePair(node_parent, node_child,
                             atan2_approximation1(node_child->getY() -
                                                  node_parent->getY
                                                          (),
                                                  node_child->getX() -
                                                  node_parent->getX()));
}



float RRT::localCostNodePair(Node *node_parent, Node *node_child, float angle)
const {
    float d_node_pair = distanceNodePair(node_parent, node_child);
    return d_node_pair*params->cost_p.weight_distance + params->cost_p
                                                                .weight_steering *
                                                        powf(getAbsoluteDiff2Angles(node_parent->getAngle(), angle), 2)
                                                        /d_node_pair;
}

float RRT::distanceNodePair(Node *node_0, Node *node_1) const {
    return euclideanDistance(node_0->getX(), node_0->getY(),
                             node_1->getX(), node_1->getY())/params->expand_radius;
}

const std::vector<Node *> &Node::getChildren() const {
    return children;
}


Node::Node(float x, float y, RRT *rrt, int node_id) : x(x), y(y), line(sf::Lines, 2), rrt
        (rrt), node_id(node_id) {
    line[0].position = sf::Vector2f(x, y);
}

Node::Node(float x, float y, const sf::Color &color, RRT *rrt, int node_id)
        : x(x), y(y), line(sf::Lines, 2), rrt
        (rrt), node_id(node_id) {
    this->color = color;
    line[0].color = color;
    line[1].color = color;
    line[0].position = sf::Vector2f(x, y);
}

void Node::setParent(Node *parent_node, bool propagate) {
    if(this->parent_node) {
        this->parent_node->removeChild(this);
    }

    this->parent_node = parent_node;
    updateLocalCost();
    line[1].position = sf::Vector2f(parent_node->getX(), parent_node->getY());
    for(auto &child : children) {
        child->updateLocalCost();
    }
    if(propagate) {
        propagateCosts();
    }
}

void Node::setParent(Node *parent_node) {
    setParent(parent_node, true);
}

void Node::addChild(Node *child_node) {
    children.push_back(child_node);
}

float Node::getX() const {
    return x;
}

float Node::getY() const {
    return y;
}

void Node::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    if(parent_node) {
        target.draw(line, states);
    }
}

void Node::setLocalCost(float cost) {
    node_cost = cost;
}

float Node::getCostToNode() const {
    return cost_to_node;
}

float Node::getLocalCost() const {
    return node_cost;
}

float Node::getAngle() const {
    return angle_to_parent;
}

void Node::propagateCosts() {
    if(parent_node) {
        node_info = rrt->localInfoFunction(parent_node, this);
        path_length = parent_node->getPathLength() + 1;
        cost_to_node = parent_node->cost_to_node + node_cost;
        info_to_node = parent_node->info_to_node +
                       pow(rrt->params->info_p.gamma,path_length) * node_info;
        util_to_node = rrt->cumulativeFunction(parent_node, this);

    }
    for(auto & node : children){
        node->propagateCosts();
    }
}

void Node::removeChild(Node *node) {
    auto node_pos = std::find(children.begin(), children.end(), node);
    if(node_pos != children.end()) {
        children.erase(node_pos);
    }
}

Node *Node::getParent() const {
    return parent_node;
}

void Node::move(float x, float y) {
    move(x, y, true);
}

void Node::move(float x, float y, bool propagate) {
    this->x = x;
    this->y = y;
    for(auto &child : children) {
        child->setParent(this, false);
    }
    if(propagate) {
        propagateCosts();
    }
}

void Node::setAngle(float angle) {
    setAngle(angle, true);
}


void Node::setAngle(float angle, bool propagate) {
    angle_to_parent = angle;
    for(auto &child : children) {
        child->updateLocalCost();
    }
    if(propagate) {
        propagateCosts();
    }
}


float Node::getInfoToNode() const {
    return info_to_node;
}

float Node::getUtilToNode() const {
    return util_to_node;
}

float Node::getLocalInfo() const {
    return node_info;
}

int Node::getPathLength() const {
    return path_length;
}

bool Node::isOnPath(Node * node) const {
    auto q_node = getParent();
    while(q_node){
        if(q_node == node){
            return true;
        }
        q_node = q_node->getParent();
    }
    return false;
}

void Node::clearChilden() {
    children.clear();
}

void Node::updateLocalCost() {
    angle_to_parent = atan2_approximation1( getY() - parent_node->getY(),  getX() -
                                                                           parent_node->getX
                                                                                   ());
    auto cost = rrt->localCostNodePair(parent_node, this, angle_to_parent);
    auto info = rrt->localInfoFunction(parent_node, this);
    node_cost = cost;
    node_info = info;
    path_length = parent_node->getPathLength() + 1;
    cost_to_node = parent_node->cost_to_node + node_cost;
    info_to_node = parent_node->info_to_node +
                   pow(rrt->params->info_p.gamma,path_length) * node_info;
    util_to_node = rrt->cumulativeFunction(parent_node, this);
}


void Node::makeRoot() {
    if(parent_node) {
        parent_node->removeChild(this);
    }
    parent_node = nullptr;
    path_length = 0;
    node_cost = 0.0f;
    cost_to_node = 0.0f;
    info_to_node = 0.0f;
    util_to_node = 0.0f;
    node_info = 0.0f;
    for(auto &child : children) {
        child->updateLocalCost();
    }
//    move(x, y, false);
}

int Node::getNodeId() const {
    return node_id;
}

std::string Node::exportNodeAsString() {
    // node_id, x, y, parent_id, cost_to_node, info_to_node
    std::string export_string = "";
    export_string += std::to_string(node_id) + ",";
    export_string += std::to_string(x) + ",";
    export_string += std::to_string(y) + ",";
    if(!parent_node){
        export_string += "-1";
    }
    else{
        export_string += std::to_string(parent_node->getNodeId());
    }
    export_string += ",";
    export_string += std::to_string(cost_to_node) + ",";
    export_string += std::to_string(info_to_node) + ",";
    return export_string;
}

void RRT::killThread() {
    kill_thread=true;
    thread.join();
}

RRT::RRT(RRTParams *params, std::vector<RRTObstacle *> *obstacles_ptr, RRTField*
field, float x0,
         float y0): rrt_nodes(params), field(field), root_circle(0.2) {
    if(params->maximize_util){
        cumulativeFunction = std::bind(&RRT::cumulativeUtil, this,
                                       std::placeholders::_1,
                                       std::placeholders::_2);
        localInfoFunction = std::bind(&RRT::localInfoNodePairUD, this,
                                      std::placeholders::_1,
                                      std::placeholders::_2);
        bestPathFunction = std::bind(&RRT::getBestPathUD, this);
    }
    else{
        cumulativeFunction = std::bind(&RRT::cumulativeCost, this,
                                       std::placeholders::_1,
                                       std::placeholders::_2);
        localInfoFunction = std::bind(&RRT::localInfoNodePair, this,
                                      std::placeholders::_1,
                                      std::placeholders::_2);
        bestPathFunction = std::bind(&RRT::getBestPathN, this);
    }
    int r = rand() % 255;
    int g = rand() % 255;
    int b = rand() % 255;
    root_circle.setOrigin(0.2, 0.2);
    kdtree = new my_kd_tree_t(2, rrt_nodes);
    color = sf::Color( r, g, b );
    this->obstacles_ptr = obstacles_ptr;
    this->params = params;
    root_node = new Node(x0, y0, this, 0);
//    root_node->setAngle(M_PI);
    rrt_nodes.push_back(root_node);
    kdtree->addPoints(rrt_nodes.size() - 1, rrt_nodes.size() - 1);
}

void RRT::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    if(params->draw_all || !goal_node) {
        for (const auto &node: rrt_nodes) {
            node->draw(target, states);
        }
    }
    else{
        const Node* node = goal_node;
        while(node){
            node->draw(target, states);
            node = node->getParent();
        }
    }
    for(auto & seg : segments){
        target.draw(seg, states);
    }
    if(circ){
        target.draw(*circ, states);
    }
    for(const auto& obs:*obstacles_ptr){
        obs->draw(target, states);
    }
    target.draw(root_circle, states);
}

std::vector<Node *> RRT::nearNodes(float &x, float &y) const {
    std::vector<Node * > near_nodes;


    float query_pt[2] = {x, y};
    float radius_sq = powf(params->rewire_radius, 2);
    std::vector<std::pair<size_t, float>> indices_dists;
    nanoflann::RadiusResultSet<float, size_t> result_set(radius_sq,
                                                         indices_dists);
    kdtree->findNeighbors(result_set, query_pt, {10});
    near_nodes.reserve(indices_dists.size());
    for(auto ress : indices_dists){
        near_nodes.push_back(rrt_nodes[ress.first]);
    }
    return near_nodes;
}

bool RRT::isPointInBounds(float x, float y) const {
    return x >= 0 && x <= params->size_x && y >= 0 && y <= params->size_y;
}
bool RRT::obstacleFree(float x0, float y0, float x1, float y1) const {
    if(!isPointInBounds(x0, y0) || !isPointInBounds(x1, y1)){
        return false;
    }
    for(const auto& obs: *obstacles_ptr){
        if(doIntersect(x0, y0, x1, y1, obs->x0, obs->y0, obs->x1, obs->y1)){
            return false;
        }
    }

    return true;
}

Node *RRT::nearestNode(float &x, float &y) const {
    nanoflann::KNNResultSet<float> res_set(1);
    size_t ret_id;
    float d;
    float query_pt[2] = {x, y};
    res_set.init(&ret_id, &d);
    kdtree->findNeighbors(res_set, query_pt, {10});

    return rrt_nodes[ret_id];
}

bool RRT::isGoalNode(Node *node) const {
    if(isGoalPoint(node->getX(), node->getY())){
        return true;
    }
    return false;
}

bool RRT::isGoalPoint(float x, float y) const {
    if(euclideanDistance(x, y, goal_x, goal_y) < 2.0){
        return true;
    }
    return false;
}

void RRT::expand_() {
    float x = sampleScalar(params->size_x);
    float y = sampleScalar(params->size_y);
    // Find the nearest node for steering first
    auto nearest_node = nearestNode(x, y);
    linearInterp(nearest_node->getX(), nearest_node->getY(), x, y,
                 params->expand_radius, x, y);
    if (obstacleFree(x, y, nearest_node->getX(), nearest_node->getY())) {
        auto near_nodes = nearNodes(x, y);
        Node *added_node = new Node(x, y, color, this, rrt_nodes.size());
        rrt_nodes.push_back(added_node);
        kdtree->addPoints(rrt_nodes.size() - 1, rrt_nodes.size() - 1);
        float min_cost = cumulativeFunction(nearest_node, added_node);
        // Find best parent
        for (auto &node: near_nodes) {
            if (obstacleFree(x, y, node->getX(), node->getY())) {
                float util = cumulativeFunction(node, added_node);
                if (util < min_cost) {
                    nearest_node = node;
                    min_cost = util;
                }
            }
        }
        random_queue.push_front(added_node);
        added_node->setParent(nearest_node);
        nearest_node->addChild(added_node);
    }
}


float RRT::getBestCost() const {
    return goal_node->getCostToNode();
}



RRT::~RRT() {
    kill_thread=true;
    thread.join();
    for(auto& node: rrt_nodes){
        delete node;
    }
    delete kdtree;
}

const Node *RRT::getGoalNode() const {
    return goal_node;
}

const Node *RRT::getRootNode() const {
    return root_node;
}

void RRT::buildTree() {
    while(rrt_nodes.size() < params->num_nodes) {
        expand_();
    }
    if(params->static_tree){
        for(auto& node: rrt_nodes){
            auto x = node->getX();
            auto y = node->getY();
            node->neighboring_nodes = nearNodes(x, y);
        }
    }
    thread =std::thread([&]{
        while(!kill_thread){
            std::lock_guard<std::mutex> lock(mutex);
            if(!random_queue.empty() || !root_queue.empty()){
                done = false;
                for(int i = 0; i < params->batch_size; i++) {
                    rewireFromRoot();
                    rewireRandom();
                    if (random_queue.empty() && root_queue.empty() &&!update_flag){
                        done = true;
                        break;
                    }
                }
            }
            else{
                if(update_flag) {
                    if(!root_offset_queue.empty()){
                        auto node = root_offset_queue.top();
                        root_offset_queue.pop();
                        makeRootNode(node);
                    }
                    else {
                        if (update_once) {
                            update_once = false;
                            rewired_set.clear();
//                            root_node->updateLocalCost();
                            root_node->propagateCosts();
                            root_queue.push_back(root_node);
                        }
                        else{
                            done = true;
                        }
                    }
                }
            }

        }
    });
}

void RRT::rewireRandom() {

    if(random_queue.empty()){
        return;
    }
    auto node = random_queue.pop_front();
    float x = node->getX();
    float y = node->getY();
    // Rewiring
    std::vector<Node *> near_nodes;
    if(params->static_tree){
        near_nodes = node->neighboring_nodes;
    } else {
        near_nodes = nearNodes(x, y);
    }
    for (auto &near_node: near_nodes) {
        if(!node->isOnPath(near_node)) {
            if (near_node != node->getParent() and near_node != root_node) {
                // Order the if statements to minimize the number of calls to
                // obstacleFree and localCostNodePair
                if (obstacleFree(x, y, near_node->getX(), near_node->getY())) {
                    if (cumulativeFunction(node, near_node) <
                        near_node->getUtilToNode()) {
                        near_node->setParent(node);
                        node->addChild(near_node);
                        random_queue.push_back(near_node);
                    }
                }
            }
        }
    }
}

void RRT::rewireFromRoot() {
    if(root_queue.empty()) {
        return;
    }
//    std::cout << "root: " << root_queue.size() << std::endl;
    auto node = root_queue.pop_front();
    rewired_set.insert(node);
    float x = node->getX();
    float y = node->getY();

    // Rewiring
    std::vector<Node*> near_nodes;
    std::vector<Node*> near_nodes2;
    if(params->static_tree){
        near_nodes = node->neighboring_nodes;
    } else {
        near_nodes = nearNodes(x, y);
    }
    for (auto &near_node: near_nodes) {
        if (near_node != node->getParent() && near_node != root_node &&
            near_node != node) {
            if (!node->isOnPath(near_node)) {
                if (obstacleFree(x, y, near_node->getX(),
                                 near_node->getY())) {
                    if (cumulativeFunction(node, near_node) <
                        near_node->getUtilToNode()) {
                        near_node->setParent(node);
                        node->addChild(near_node);
                    }
                }
            }

            if (!rewired_set.contains(near_node)) {
                root_queue.push_back(near_node);
            }
        }
    }
}

Node* RRT::getPathTo(float x, float y) {
    std::lock_guard<std::mutex> lock(mutex);
    Node* nearest_node = nearestNode(x, y);
    goal_node = nearest_node;
    return goal_node;
}

Node *RRT::getBestPathUD() {
    std::lock_guard<std::mutex> lock(mutex);
    float min_cost = std::numeric_limits<float>::infinity();
    for(auto & node: rrt_nodes){
        if(node->getPathLength() > params->min_tail_length + params->offset_replan){
            float util = node->getUtilToNode();
            if (util < min_cost) {
                min_cost = util;
                goal_node = node;
                root_move_counter = 0;
            }
        }
    }
    return goal_node;
}

Node* RRT::getBestPathN() {
    std::lock_guard<std::mutex> lock(mutex);
    float min_cost = std::numeric_limits<float>::infinity();
    for(auto & node: rrt_nodes){
        if(!node->getChildren().size()) {
            if(node->getPathLength() > params->min_tail_length + params->offset_replan){
                float util = node->getCostToNode()/node->getInfoToNode();
                if (util < min_cost) {
                    min_cost = util;
                    goal_node = node;
                    root_move_counter = 0;
                }
            }
        }
    }
    return goal_node;
}

void RRT::makeRootNode(Node *node) {
//    auto x = node->getX();
//    auto y = node->getY();
//    root_circle.setPosition(x, y);
//    std::cout << "moving root node to " << x << ", " << y << std::endl;
//    float angle = node->getAngle();
//    float root_x = root_node->getX();
//    float root_y = root_node->getY();
//    rewired_set.clear();
//    std::vector<Node*> old_neighbors = root_node->neighboring_nodes;
//    node->neighboring_nodes = old_neighbors;
////        new_node->neighboring_nodes.push_back(root_node);
//    for (auto &old_neighbor: old_neighbors) {
//        // TODO this might add them twice
//        old_neighbor->neighboring_nodes.push_back(node);
//    }
//    root_node->neighboring_nodes = nearNodes(x, y);  // Update the
//    // neighboring nodes
//    std::vector<Node*> children = root_node->getChildren();
//    std::vector<Node*> children_node = node->getChildren();
//    node->clearChilden();
//    root_node->clearChilden();
//    node->move(root_x, root_y, false);  // Move the node
//    root_node->move(x, y, false);  // Move the node
//    root_node->setAngle(angle, false);
//    node->setParent(root_node, false);
//    for (auto &child: children_node) {
//        child->setParent(root_node, false);
//        root_node->addChild(child);
//    }
//    for (auto &child: children) {
//        child->setParent(node, false);
//        node->addChild(child);
//    }
//    root_node->addChild(node);
//    root_node->propagateCosts();
//    root_queue.clear();
//    root_queue.push_back(root_node);

    root_move_counter += 1;
    root_circle.setPosition(node->getX(), node->getY()); // Move the root viz

    rewired_set.clear();
    root_node->removeChild(node);
    node->makeRoot(); // future root node now has root attributes!
    root_node->setParent(node, false);  // This makes the old root node a
    // child of the new root node, calculates the local costs of the old
    // root's children
    node->addChild(root_node);  // make the old root node a child of the new
    std::swap(root_node, node);  // Swap the root node and the new root node
//    std::cout << "COSTS ROOT " << root_node->getCostToNode() << std::endl;
//    std::cout << "COSTS OLD ROOT " << node->getCostToNode() << std::endl;
//    root_node->move(root_node->getX(), root_node->getY(), false);  // Update
    // all children of the new root node
    root_node->propagateCosts();  // Propagate costs through the tree
    root_queue.clear();
    root_queue.push_back(root_node); // Add the new root node to the queue
//    std::cout << "moving root " << std::endl;
}


void RRT::moveRoot(float x, float y, float angle) {
    if(x != root_node->getX() || y != root_node->getY() || angle != root_node->getAngle()) {
        done = false;
        std::lock_guard<std::mutex> lock(mutex);
        rewired_set.clear();
        // Create new node as a clone of the old root node.
        Node *new_node = new Node(root_node->getX(), root_node->getY(), color,
                                  this, rrt_nodes.size());
        //    new_node->setAngle(root_node->getAngle());
        rrt_nodes.push_back(new_node);
        auto old_neighbors = root_node->neighboring_nodes;
        new_node->neighboring_nodes = old_neighbors;
//        new_node->neighboring_nodes.push_back(root_node);
        for (auto &node: old_neighbors) {
            node->neighboring_nodes.push_back(new_node);
        }
        // Add the new node to the kdtree
        kdtree->addPoints(rrt_nodes.size() - 1, rrt_nodes.size() - 1);
        root_node->move(x, y);  // Move the node
        root_node->setAngle(angle);
        root_node->neighboring_nodes = nearNodes(x, y);  // Update the
        // neighboring nodes
        auto children = root_node->getChildren();
        new_node->setParent(root_node);
        for (auto &child: children) {
            child->setParent(new_node);
            new_node->addChild(child);
        }
        root_node->clearChilden();
        root_node->addChild(new_node);
        root_queue.clear();
        root_queue.push_back(root_node);
        root_node->propagateCosts();
    }
}

bool RRT::rewireDone() const {
    return done;
}

void RRT::buildPathWithCorridor(Node *goal_node) {
    segments.clear();
    auto node = goal_node;
    while(node){
        auto parent = node->getParent();
        if(parent){
            segments.emplace_back(node->getX(), node->getY(), parent->getX(), parent->getY());
        }
        node = parent;
    }
}

const std::vector<Segment> *RRT::getSegments() const {
    return &segments;
}

float RRT::localInfoNodePair(Node *node_parent, Node *node_child) const {
    // TODO Reuse distance calculation from localCostNodePair
    auto dist = distanceNodePair(node_parent, node_child);
    auto entropy_0 = field->getUncertainty(node_parent->getX(),
                                            node_parent->getY())
                        + params->info_p.weight_source
                        * field->getSourceLikelihood(node_parent->getX(),
                                                    node_parent->getY());
    auto entropy_1 = field->getUncertainty(node_child->getX(),
                                            node_child->getY())
                        + params->info_p.weight_source
                        * field->getSourceLikelihood(node_child->getX(),
                                                     node_child->getY());
    auto entropy = (entropy_0 + entropy_1)/2.0;
    return dist * entropy * params->info_p.discount_info;
}

float RRT::cumulativeCost(Node *node_parent, Node *node_child) const {
    float local_cost = localCostNodePair(node_parent, node_child);
    float local_info = localInfoFunction(node_parent, node_child);
    float cost = node_parent->getCostToNode() + local_cost;
    float info = node_parent->getInfoToNode() +
                 pow(params->info_p
                             .gamma, node_parent->getPathLength() + 1)*local_info;
//    return powf(cost, params->info_p.alpha)/info;
    auto return_value = node_parent->getUtilToNode() +
                        pow(local_cost, params->info_p.alpha)/(1.0 +
                        local_info);
    return return_value;
}


void RRT::updateTree() {
    done = false;
    update_flag = true;
}

void RRT::updateFromRoot() {
    update_once = true;
}

float RRT::getSizeX() const {
    return params->size_x;
}

float RRT::getSizeY() const {
    return params->size_y;
}

void RRT::updateFromRoot(float goal_angle) {
//    if(root_node->getAngle() != goal_angle) {
//        std::lock_guard<std::mutex> lock(mutex);
//        root_node->setAngle(goal_angle);
//    }
    update_once = true;
}

float RRT::localInfoNodePairUD(Node *node_parent, Node *node_child) const {
    // TODO Reuse distance calculation from localCostNodePair
    // Check if the current node is part of the parent's path
    Node* node = node_parent->getParent();
    auto neighboring_nodes = node_child->neighboring_nodes;
    while(node){
        if (std::find(neighboring_nodes.begin(), neighboring_nodes.end(),
                      node) !=
            neighboring_nodes.end
                    ()) {
            return 0.01;
        }
        node = node->getParent();
    }
    auto dist = distanceNodePair(node_parent, node_child);
    auto entropy_0 = (field->getUncertainty(node_parent->getX(),
                                            node_parent->getY()));
    auto entropy_1 = (field->getUncertainty(node_child->getX(),
                                            node_child->getY()));
    auto entropy = (entropy_0 + entropy_1)/2.0;
    return dist * entropy + 0.01;
}

float RRT::cumulativeUtil(Node *node_parent, Node *node_child) const {
    float local_cost = localCostNodePair(node_parent, node_child);
    float local_info = localInfoFunction(node_parent, node_child);
    float cost = node_parent->getCostToNode() + local_cost;
    float info = node_parent->getInfoToNode() +
                 pow(params->info_p
                             .gamma, node_parent->getPathLength() + 1)*local_info;
//    return powf(cost, params->info_p.alpha)/info;
    auto return_value =
            -info/pow(cost, params->info_p.alpha);
    return return_value;
}

Node *RRT::getBestPath() {
    return bestPathFunction();
}

void RRT::offsetRoot(int steps) {
    std::lock_guard<std::mutex> lock(mutex);
    while(!root_offset_queue.empty()){
        root_offset_queue.pop();
    }
    Node* node = goal_node;
    while(node->getPathLength() > steps){
        node = node->getParent();
    }
    while(node->getParent()){
        root_offset_queue.push(node);
        node = node->getParent();
    }
}

void RRT::exportTree(std::string path) {
//    if(exists(path)){
//        throw std::runtime_error("Logging into existing file! Clear first!");
//    }
    std::fstream file;
    file.open(path, std::ios::out);
    // node_id, x, y, parent_id, cost_to_node, info_to_node
    file << "node_id,x,y,parent_id,cost_to_node,info_to_node" << std::endl;
    for(auto node: rrt_nodes){
        file << node->exportNodeAsString() << std::endl;
    }
}


RewireQueue::RewireQueue() {

}

void RewireQueue::push_back(Node *node) {
    if(!set.contains(node)){
        queue.push_back(node);
        set.insert(node);
    }
}

void RewireQueue::push_front(Node *node) {
    if(!set.contains(node)){
        queue.push_front(node);
        set.insert(node);
    }
}

Node *RewireQueue::pop_front() {
    Node* node = queue.front();
    queue.pop_front();
    set.erase(node);
    return node;
}

Node *RewireQueue::pop_back() {
    Node* node = queue.back();
    queue.pop_back();
    set.erase(node);
    return node;
}

void RewireQueue::clear() {
    set.clear();
    queue.clear();
}

Segment::Segment(float x0, float y0, float x1, float y1){
    this->x0 = x0;
    this->y0 = y0;
    this->x1 = x1;
    this->y1 = y1;

}

void Segment::draw(sf::RenderTarget &target, sf::RenderStates states) const {
    float angle = atan2f(y1 - y0, x1 - x0) + M_PI/2.0f;

    sf::Vertex line[] =
            {
                    sf::Vertex(sf::Vector2f(x0, y0), sf::Color::Green),
                    sf::Vertex(sf::Vector2f(x1, y1), sf::Color::Green)
            };
    sf::Vertex line2[] =
            {
                    sf::Vertex(sf::Vector2f(x0 + cosf(angle) * upper_bound,
                                            y0 + sinf(angle) * upper_bound),
                               sf::Color::Green),
                    sf::Vertex(sf::Vector2f(x1 + cosf(angle) * upper_bound,
                                            y1 + sinf(angle) * upper_bound),
                               sf::Color::Green)
            };
    sf::Vertex line3[] =
            {
                    sf::Vertex(sf::Vector2f(x0 + cosf(angle) * lower_bound,
                                            y0 + sinf(angle) * lower_bound),
                               sf::Color::Green),
                    sf::Vertex(sf::Vector2f(x1 + cosf(angle) * lower_bound,
                                            y1 + sinf(angle) * lower_bound),
                               sf::Color::Green)
            };
    target.draw(line, 2, sf::Lines);
//    target.draw(line2, 2, sf::Lines);
//    target.draw(line3, 2, sf::Lines);
}



float SineTestField::getUncertainty(float x, float y) const {
//    return (sin(y)*cos(x)+ 1.0)*2.0;
//    return (sin(y) + 1.0)*(sin(x) + 1.0);
    return (sin(y)*sin(x - 5.5) +1.0);
}

SineTestField::SineTestField() = default;
