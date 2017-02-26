#ifndef COSTFUNCTION_H
#define COSTFUNCTION_H

#include "config.h"

class CostFunction
{
public:
    CostFunction();
private:
protected:
    // attributes //
public:
private:

protected:
    double dt;
    stateVec_t lx;
    stateMat_t lxx;
    commandVec_t lu;
    commandMat_t luu;
    commandR_stateC_t lux;
    stateR_commandC_t lxu;
    stateMat_t Q;
    stateMat_t Qf;
    commandMat_t R;
    // methods //
public:
    virtual void computeAllCostDeriv(const stateVec_t& X, const commandVec_t& U)=0;
    virtual void computeFinalCostDeriv(const stateVec_t& X)=0;
private:
protected:
    // accessors //
public:
    stateVec_t& getlx();
    stateMat_t& getlxx();
    commandVec_t& getlu();
    commandMat_t& getluu();
    commandR_stateC_t& getlux();
    stateR_commandC_t& getlxu();
    stateMat_t& getQ();
    stateMat_t& getQf();
    commandMat_t& getR();
};

#endif // COSTFUNCTION_H
