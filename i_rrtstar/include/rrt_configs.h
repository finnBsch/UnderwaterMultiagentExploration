//
// Created by finn on 6/15/23.
//

#ifndef MASTERMAIN_RRT_CONFIGS_H
#define MASTERMAIN_RRT_CONFIGS_H
#include <cmath>
#include <string>

struct PathCostParams{
    float weight_distance;
    float weight_steering;
};

struct RRTInformationParams{
    float gamma;
    float alpha;
    float weight_mean;
    float weight_source;
    float discount_info;
};

struct RRTParams{
    int batch_size = 1;
    bool static_tree = true;
    bool draw_all = false;
    int num_nodes = 20000;
    float expand_radius = 0.2;
    float rewire_radius = .6f;
    float rewire_radius_sq;
    float size_x = 10.0f;
    float size_y = 10.0f;
    int offset_replan;
    int min_tail_length;
    bool maximize_util = false;

    RRTInformationParams info_p;
    PathCostParams cost_p;
    RRTParams();
    RRTParams(std::string path);
    void exportConfig(std::string path);
};

#endif //MASTERMAIN_RRT_CONFIGS_H
