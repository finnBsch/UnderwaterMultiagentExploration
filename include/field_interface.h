//
// Created by finn on 6/10/23.
//

#ifndef MASTERMAIN_FIELD_INTERFACE_H
#define MASTERMAIN_FIELD_INTERFACE_H
#include "gmrf/gmrf.h"
#include "info_field.h"

#include "rrt.h"
#include <eigen3/Eigen/Dense>

class FieldInterface : public InfoField, public RRTField {
private:
    GMRF* gmrf;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> uncertainty;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> source_likelihoods;
public:
    void update();

    double getValue(double x, double y) override;

    float getUncertainty(float x, float y) const override;

    std::array<double, 2> getGradient(double x, double y) override;

    explicit FieldInterface(GMRF* gmrf);

    double totalUncertainty();

    float getSourceLikelihood(float x, float y) const override;
};


#endif //MASTERMAIN_FIELD_INTERFACE_H
