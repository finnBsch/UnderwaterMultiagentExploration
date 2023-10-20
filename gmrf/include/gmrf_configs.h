//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//

#ifndef GMRF_CONFIGS_H
#define GMRF_CONFIGS_H
#include <cmath>
#include <string>

struct GmrfParams{
    // General params
    int N_X = 10;
    int N_Y = 10;
    double size_x = 10.0f;
    double size_y = 10.0f;

    // Real field params
    double prior_mean = 0.0f;
    double prior_sigma = 100.0f;
    // double sigma_p = 1.0f;
    double sigma_p_per_m = 1.0f;
    double sigma_zeta = 0.0f*.01f;
    double sigma_s = 1.0f;
    double sigma_flow = 1.0f;

    // Flow field params
    double ff_prior_mean_x = 0.0f;
    double ff_prior_mean_y = 0.0f;
    double ff_prior_sigma = 100.0f;
    double ff_sigma_zeta = 0.0f*.01f;
    // double ff_sigma_p = 1.0f;
    double ff_sigma_p_per_m = 1.0f;
    double ff_sigma_s = 1.0f;
    double ff_sigma_c = 1.0f;
    double ff_sigma_o = 1.0f;

    // Energy function config
    bool field_neighbors = true;
    bool flow_neighbors = true;
    bool flow_conservation = true;
    bool obstacle_flow = true;
    bool flow_field = true;

    bool all_neighbors = false;
    GmrfParams();
    GmrfParams(std::string path);
    void exportConfig(std::string path);
};


/**
 * Parameters describing the viz.
 * @param res_x: Window size in x
 * @param res_y: Window size in y. If 0, will be automatically determined to match x.
 * @param margin: proportion of width as measure for margin
 * @param s_x: size of simulation in x
 * @param s_y: size of simulation in y
 * @param infobox: proportion of width used for infobox.
 */
struct Viz_Params{
    bool draw_map = true;
    bool draw_uncertainty = true;
    bool draw_gt = true;
    int res_x = 1920;
    int res_y = 0;
    int offset_x = 0;
    int offset_y = 0;
    double margin = 0.01;
    double infobox = 0.1;
    Viz_Params();
};


#endif //GMRF_CONFIGS_H
