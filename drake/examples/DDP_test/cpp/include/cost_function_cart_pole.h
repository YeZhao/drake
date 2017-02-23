#ifndef COSTFUNCTIONCARTPOLE_H
#define COSTFUNCTIONCARTPOLE_H

#include "config.h"

#include "costfunction.h"

#include <Eigen/Dense>

using namespace Eigen;

class CostFunctionCartPole : public CostFunction
{
public:
    CostFunctionCartPole();
private:
    stateMat_t Q;
    stateMat_t Qf;
    commandMat_t R;
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
    stateMat_t getQ();
    stateMat_t getQf();
    commandMat_t getR();

private:
protected:
    // accessors //
public:

};

#endif // COSTFUNCTIONCARTPOLE_H
