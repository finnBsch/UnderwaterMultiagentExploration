//
// Created by finn on 5/29/23.
//

#ifndef SQP_TEST_FIELD_H
#define SQP_TEST_FIELD_H
#include <array>
#include <casadi/casadi.hpp>

using namespace casadi;

class InfoField{
private:
public:
    virtual std::array<double, 2> getGradient(double x, double y) = 0;
    virtual double getValue(double x, double y) = 0;
};


class TestField : public InfoField{
private:
    Function val_fun;
    Function grad_fun;
public:
    TestField();
    std::array<double, 2> getGradient(double x, double y) override;
    double getValue(double x, double y) override;
};



// f(x, y) = sin(x) * cos(y) + x^2
template<typename T> T field_val(T x, T y){
    return sin(x * 2) * cos(y * 2);//+ pow(x, 2) + pow(y, 2) * 0.1;
//    return  -x-y;
}


#endif //SQP_TEST_FIELD_H
