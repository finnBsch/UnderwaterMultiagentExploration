//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//
#include <cmath>
#include "color_mixer.h"
#include <iostream>

color_int mix_color(color_int& col0, color_int&col1, double mix){
    color col0_ = {(double)col0[0]/255.0f,(double)col0[1]/255.0f, (double)col0[2]/255.0f};
    color col1_ = {(double)col1[0]/255.0f,(double)col1[1]/255.0f, (double)col1[2]/255.0f};
    return mix_color(col0_, col1_, mix);
}
color_int mix_color(color& col0, color& col1, double mix){
    sRGBInvCompanding(col0);
    sRGBInvCompanding(col1);
    double r = lin_interp(col0[0], col1[0], mix);
    double g = lin_interp(col0[1], col1[1], mix);
    double b = lin_interp(col0[2], col1[2], mix);
    double gamma = 0.43;
    double b0 = compute_brightness(col0, gamma);
    double b1 = compute_brightness(col1, gamma);
    double brightness = lin_interp(b0, b1, mix);
    double intensity = pow(brightness, 1.0f/gamma);
    if((r+g+b) != 0){
        double fac = (intensity/(r + g + b));
        r *= fac;
        g *= fac;
        b *= fac;
    }
    color ret_col = {r, g, b};
    sRGBCompanding(ret_col);
    color_int col = {(int)(ret_col[0]*255), (int)(ret_col[1]*255), (int)(ret_col[2]*255)};
    return col;
}
void sRGBCompanding(color& col){
    for(int i = 0; i < 3; i++){
        if(col[i] <= 0.003130){
            col[i] *= 12.92;
        }
        else{
            col[i] = 1.055f * (pow(col[i], (1/2.4))) - 0.055f;
        }
    }
}

double lin_interp(double x0, double x1, double mix){
    return x0 + (x1 - x0)*mix;
}

double compute_brightness(color& col, double gamma){
    return pow(col[0] + col[1] + col[2], gamma);
}

void sRGBInvCompanding(color& col){
    for(int i = 0; i < 3; i++){
        if(col[i] <= 0.04045){
            col[i] /= 12.92;
        }
        else{
            col[i] = pow((col[i] + 0.055f)/1.055f, 2.4);
        }
    }
}