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
    stateMat_t Q;
    stateMat_t Qf;
    commandMat_t R;

    stateVecTab_t cx_new; 
    commandVecTab_t cu_new; 
    stateMatTab_t cxx_new; 
    commandR_stateC_tab_t cux_new; 
    commandMatTab_t cuu_new;
    double c_new;
public:
private:
protected:
    // accessors //
public:
    stateMat_t& getQ();
    stateMat_t& getQf();
    commandMat_t& getR();
    stateVecTab_t& getcx();
    commandVecTab_t& getcu();
    stateMatTab_t& getcxx();
    commandR_stateC_tab_t& getcux();
    commandMatTab_t& getcuu();
    double& getc();
    unsigned int N;
};

#endif // COSTFUNCTION_H
