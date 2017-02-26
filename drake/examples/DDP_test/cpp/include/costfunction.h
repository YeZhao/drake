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

    stateVecTab_t cx_new; 
    commandVecTab_t cu_new; 
    stateMatTab_t cxx_new; 
    commandR_stateC_tab_t cux_new; 
    commandMatTab_t cuu_new;
    double c_new;

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
    stateVecTab_t& getcx();
    commandVecTab_t& getcu();
    stateMatTab_t& getcxx();
    commandR_stateC_tab_t& getcux();
    commandMatTab_t& getcuu();
    double& getc();
    unsigned int T;
};

#endif // COSTFUNCTION_H
