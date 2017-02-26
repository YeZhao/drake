#ifndef COSTFUNCTIONCARTPOLE_H
#define COSTFUNCTIONCARTPOLE_H

#include "config.h"
#include "costfunction.h"
#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

class CostFunctionCartPole : public CostFunction
{
public:
    CostFunctionCartPole();
private:
    double dt;
protected:
    // attributes //
public:
private:

protected:
    // methods //
public:
    void computeAllCostDeriv(const stateVec_t& X, const commandVec_t& U);
    void computeFinalCostDeriv(const stateVec_t& X);
private:
protected:
    // accessors //
public:

};

#endif // COSTFUNCTIONCARTPOLE_H