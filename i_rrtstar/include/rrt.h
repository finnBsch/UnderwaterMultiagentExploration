//
// Created by finn on 6/6/23.
//

#ifndef I_RRTSTAR_RRT_H
#define I_RRTSTAR_RRT_H
#include <vector>
#include "util.h"
#include <SFML/Graphics.hpp>
#include <thread>
#include <atomic>
#include <nanoflann.hpp>
#include <deque>
#include <set>
#include <stack>
#include "rrt_configs.h"
class NodeWrapper;

// Typedefs

using my_kd_tree_t = nanoflann::KDTreeSingleIndexDynamicAdaptor<
        nanoflann::L2_Simple_Adaptor<float, NodeWrapper>,
        NodeWrapper, 2 /* dim */>;
// Forward declaring for convenience
class Node;
class RRT;
class RRTObstacle;
class Viz;

class RRTObstacle: public sf::Drawable{
protected:
    sf::VertexArray line;
public:
    float x0;
    float y0;
    float x1;
    float y1;
    float length;
    float angle;
    RRTObstacle(float x0, float y0, float x1, float y1);
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    double getDistance(float x, float y) const;

};

class RRTField{
private:
public:
    virtual float getUncertainty(float x, float y) const = 0;
    virtual float getSourceLikelihood(float x, float y) const = 0;
};

class SineTestField : public RRTField{
private:
public:
    SineTestField();
    float getUncertainty(float x, float y) const;
    float getSourceLikelihood(float x, float y) const {
        return 0;
    }
};

class RewireQueue {
private:
    std::deque<Node*> queue;
    std::set<Node*> set;
public:
    RewireQueue();
    void push_back(Node* node);
    void push_front(Node* node);
    void clear();
    Node* pop_front();
    Node* pop_back();
    size_t size(){return queue.size();}
    bool empty(){return queue.empty();}
};

inline float atan2_approximation1(float y, float x)
{
//http://pubs.opengroup.org/onlinepubs/009695399/functions/atan2.html
//Volkan SALMA

    const float ONEQTR_PI = M_PI / 4.0;
    const float THRQTR_PI = 3.0 * M_PI / 4.0;
    float r, angle;
    float abs_y = fabs(y) + 1e-10f;      // kludge to prevent 0/0 condition
    if ( x < 0.0f )
    {
        r = (x + abs_y) / (abs_y - x);
        angle = THRQTR_PI;
    }
    else
    {
        r = (x - abs_y) / (x + abs_y);
        angle = ONEQTR_PI;
    }
    angle += (0.1963f * r * r - 0.9817f) * r;
    if ( y < 0.0f )
        return( -angle );     // negate if in quad III or IV
    else
        return( angle );
}

class Node: public sf::Drawable {
private:
    int node_id;
    bool parent_changed = false;
    RRT* rrt;
    float x;
    float y;
    sf::Color color;
    float node_time = 0.0f;
    float node_cost = 0.0f;
    float cost_to_node = 0.0f;
    float info_to_node = 0.0f;
    int path_length = 0;
    float util_to_node = 0.0f;
    float node_info = 0.0f;
    Node* parent_node = nullptr;
    std::vector<Node*> children;
    sf::VertexArray line;
    float angle_to_parent = 0.0f;
public:
    std::vector<Node*> neighboring_nodes;
    Node(float x, float y, RRT *rrt, int node_id);
    Node(float x, float y, const sf::Color &color, RRT *rrt,int node_id);
    void setParent(Node* parent_node);
    void setParent(Node* parent_node, bool propagate);
    void updateLocalCost();
    void addChild(Node* child_node);
    void setLocalCost(float cost);
    void propagateCosts();
    void removeChild(Node* node);
    void clearChilden();
    int getPathLength() const;
    float getCostToNode() const;
    float getInfoToNode() const;
    float getUtilToNode() const;
    float getLocalCost() const;
    float getLocalInfo() const;
    float getX() const;
    float getY() const;
    float getAngle() const;
    Node* getParent() const;
    int getNodeId() const;
    const std::vector<Node*>& getChildren() const;
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    void move(float x, float y);
    void move(float x, float y, bool propagate);
    void setAngle(float angle);
    void setAngle(float angle, bool propagate);
    const std::set<Node*>& getNodesOnPath() const;
    bool isOnPath(Node* node) const;
    void makeRoot();
    std::string exportNodeAsString();
};


class NodeWrapper: public std::vector<Node*> {
private:
    RRTParams* params;
public:
    NodeWrapper(RRTParams* params){this->params = params;};
    size_t kdtree_get_point_count() const{
        return size();
    }
    inline float kdtree_get_pt(const size_t idx, int dim) const{
        if(dim == 0) {
            return this->operator[](idx)->getX();
        }
        else{
            return this->operator[](idx)->getY();
        }
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX &bb) const
    {
        bb[0].low = 0; bb[0].high = params->size_x;  // 0th dimension limits
        bb[1].low = 0; bb[1].high = params->size_y;  // 1st dimension limits
        return true;
    }

};


class Segment : public sf::Drawable {
private:
public:

    float x0;
    float y0;
    float x1;
    float y1;
    float upper_bound;
    float lower_bound;
    Segment(float x0, float y0, float x1, float y1);
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
};

class RRT: public sf::Drawable {
    friend class Node;
private:
    sf::CircleShape root_circle;
    void makeRootNode(Node* node);
    RRTField* field;
    // Class Members
    RRTParams* params{};
    sf::Color color;
    sf::CircleShape* circ = nullptr;

    // RRT Members
    std::vector<Segment> segments;
    std::atomic<bool> kill_thread{false};
    std::atomic<bool> done{false};
    std::atomic<bool> update_once{false};
    std::atomic<bool> update_flag{false};
    std::set<Node*> rewired_set;
    RewireQueue root_queue;
    RewireQueue random_queue;
    Node* root_node{};
    Node* goal_node{};
    NodeWrapper rrt_nodes;  // Store all nodes. While all nodes are
    // implicitly stored by owning the root node, it is beneficial to have
    // them in a vector to allow for faster iteration, e.g. for 1NN-search

    my_kd_tree_t* kdtree;
    std::vector<RRTObstacle*>* obstacles_ptr;

    float goal_x;
    float goal_y;
    // RRT-specific util functions
    float localCostNodePair(Node* node_parent, Node* node_child) const;
    float localCostNodePair(Node* node_parent, Node* node_child, float angle) const;
    float localInfoNodePair(Node* node_parent, Node* node_child) const;
    float localInfoNodePairUD(Node* node_parent, Node* node_child) const;
    float cumulativeCost(Node* node_parent, Node* node_child) const;  // for
    // cost-driven
    float cumulativeUtil(Node* node_parent, Node* node_child) const; // for
    // util-driven
    float distanceNodePair(Node* node_0, Node* node_1) const;

    Node* getBestPathN();  // Cost driven
    Node* getBestPathUD();
    std::function<float(Node*, Node*)> cumulativeFunction;
    std::function<float(Node*, Node*)> localInfoFunction;
    std::function<Node*()> bestPathFunction;


    Node* nearestNode(float& x, float & y) const;
    std::vector<Node*> nearNodes(float& x, float &y) const;
    bool isPointInBounds(float x, float y) const;
    bool obstacleFree(float x0, float y0, float x1, float y1) const;

    bool isGoalNode(Node* node) const;
    bool isGoalPoint(float x, float y) const;

    // Private runtime functions
    std::stack<Node*> root_offset_queue;
    int root_move_counter = 0;
    void expand_();  // Expands the tree by one node
    void rewireFromRoot();  // Rewires the tree from the root node, using the
                    // root Queue
    void rewireRandom();  // Rewires the tree from a random node, using the
                    // random Queue
    std::mutex mutex;
    std::thread thread;

public:
    void killThread();
    RRT(RRTParams* params, std::vector<RRTObstacle*>* obstacles_ptr, RRTField*
    field, float x0, float y0);
    ~RRT() override;
    void draw(sf::RenderTarget &target, sf::RenderStates states) const override;
    void exportTree(std::string path);

    // RRT Runtime functions
    void buildTree();  // Builds the tree with max nodes
    Node* getPathTo(float x, float y);  // Returns a path to the given point
    void moveRoot(float x, float y, float angle);  // Moves the root node to
    Node* getBestPath();
    // the given point
    void updateTree();
    void offsetRoot(int steps);
    void buildPathWithCorridor(Node* goal_node);
    const std::vector<Segment>* getSegments() const;
    float getBestCost() const;
    const Node* getGoalNode() const;
    const Node* getRootNode() const;

    // Temporary
    void updateFromRoot();
    void updateFromRoot(float goal_angle);
    bool rewireDone() const;

    float getSizeX() const;
    float getSizeY() const;

};
#endif //I_RRTSTAR_RRT_H
