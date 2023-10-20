//
// Created by finn on 6/10/23.
//
#include "info_field.h"

std::array<double, 2> TestField::getGradient(double x, double y) {
    auto g = grad_fun(DMDict {{"i0",x},{"i1", y}});
    auto gum = g["o0"];

    return std::array<double, 2>({gum(0).scalar(), gum(1).scalar()});
}

double TestField::getValue(double x, double y) {
    auto v = val_fun(DMDict {{"i0",x},{"i1", y}});

    return v["o0"].scalar();
}

TestField::TestField() {
    MX x = MX::sym("x");
    MX y = MX::sym("y");
    MX f = sin(x * 2) * cos(y * 2);
    val_fun = Function("val_fun", {x, y}, {f});
    auto vars = vertcat(x, y);
    grad_fun = Function("grad_fun", {x, y}, {jacobian(f, vars)});
}