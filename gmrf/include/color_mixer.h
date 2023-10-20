//
// Created by Finn Lukas Busch
// finn.lukas.busch@gmail.com
//

#ifndef GMRF_COLOR_MIXER_H
#define GMRF_COLOR_MIXER_H
#include <array>

using color = std::array<double, 3>;
using color_int = std::array<int, 3>;
void sRGBInvCompanding(color& col);
void sRGBCompanding(color& col);
double lin_interp(double x0, double x1, double mix);
double compute_brightness(color& col, double gamma);
color_int mix_color(color& col0, color& col1, double mix);
color_int mix_color(color_int& col0, color_int& col1, double mix);


#endif //GMRF_COLOR_MIXER_H
