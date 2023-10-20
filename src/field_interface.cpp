//
// Created by finn on 6/10/23.
//
#include "field_interface.h"

FieldInterface::FieldInterface(GMRF *gmrf) {
    this->gmrf = gmrf;
    update();
    // Get dimensions
}

void FieldInterface::update() {
    // Get from estimates map. This copies the current estimates to allow for flexible updates.
    this->uncertainty = this->gmrf->getCovariances();
    this->source_likelihoods = this->gmrf->getSourceLikelihoods();
}

std::array<double, 2> FieldInterface::getGradient(double x, double y) {
    return gmrf->uncertaintyGradient(x, y, this->uncertainty);
}

double FieldInterface::getValue(double x, double y) {
    return gmrf->interpolateValue(x, y, this->uncertainty);
}

float FieldInterface::getUncertainty(float x, float y) const {
    return gmrf->interpolateValue(x, y, this->uncertainty);
}

double FieldInterface::totalUncertainty() {
    return uncertainty.sum();
}

float FieldInterface::getSourceLikelihood(float x, float y) const {
    return gmrf->interpolateValue(x, y, this->source_likelihoods);
}
