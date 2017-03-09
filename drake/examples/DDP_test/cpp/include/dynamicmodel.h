#ifndef DYNAMICMODEL_H
#define DYNAMICMODEL_H

#include "config.h"
#include "costfunction.h"
#include "cost_function_cart_pole.h"

#include <Eigen/Dense>

using namespace Eigen;

class DynamicModel
{
// constructors //
public:
    //DynamicModel();

// attributes //
public:
protected:
    unsigned int stateNb;
    unsigned int commandNb;
    double dt;
    commandVec_t lowerCommandBounds;
    commandVec_t upperCommandBounds;

    stateMat_t fx;
    stateTens_t fxx;
    stateR_commandC_t fu;
    stateR_commandC_commandD_t fuu;
    stateR_stateC_commandD_t fxu;
    stateR_commandC_stateD_t fux;

    stateMatTab_t fxList;
    stateR_commandC_tab_t fuList;
    stateTensTab_t fxxList;
    stateTensTab_t fxuList;
    stateR_commandC_Tens_t fuuList;
public:

protected:
// methods 
public:
    virtual stateVec_t cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U)=0;
    virtual void cart_pole_dyn_cst_ilqr(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunction*& costFunction)=0;
    virtual void cart_pole_dyn_cst_min_output(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr, stateVec_t& xList_next, CostFunction*& costFunction)=0;
    virtual void cart_pole_dyn_cst_v3(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList, CostFunction*& costFunction)=0;
    virtual void cart_pole_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunction*& costFunction)=0;
    virtual stateVec_t update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B)=0;
    virtual void grad(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B)=0;
    virtual void hessian(const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu)=0;
    
private:
protected:
    // accessors //
public:
    unsigned int getStateNb();
    unsigned int getCommandNb();
    commandVec_t& getLowerCommandBounds();
    commandVec_t& getUpperCommandBounds();
    stateMatTab_t& getfxList();
    stateR_commandC_tab_t& getfuList();
};

#endif // DYNAMICMODEL_H