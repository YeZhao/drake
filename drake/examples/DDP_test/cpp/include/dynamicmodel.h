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
    //unsigned int T;
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
    stateTens_t fxx_new;
    stateR_commandC_commandD_t fuu_new;
    stateR_stateC_commandD_t fxu_new;
    stateR_commandC_stateD_t fux_new;
    
public:


protected:


// methods //
public:
    virtual stateVec_t cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U)=0;
    virtual void cart_pole_dyn_cst(const int& nargout, const double& dt, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, commandR_stateC_tab_t& cux, commandMatTab_t& cuu, double& c)=0;
    virtual stateVec_t update(const int& nargout, const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B)=0;
    virtual void grad(const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B)=0;
    virtual void hessian(const double& dt, const stateVec_t& X, const commandVec_t& U)=0;

    virtual stateVec_t computeNextState(double& dt, const stateVec_t& X,const stateVec_t& Xdes,const commandVec_t& U)=0;
    virtual void computeAllModelDeriv(double& dt, const stateVec_t& X,const stateVec_t& Xdes,const commandVec_t& U)=0;
    virtual stateMat_t computeTensorContxx(const stateVec_t& nextVx)=0;
    virtual commandMat_t computeTensorContuu(const stateVec_t& nextVx)=0;
    virtual commandR_stateC_t computeTensorContux(const stateVec_t& nextVx)=0;
private:
protected:
    // accessors //
public:
    unsigned int getStateNb();
    unsigned int getCommandNb();
    commandVec_t& getLowerCommandBounds();
    commandVec_t& getUpperCommandBounds();
    stateMat_t& getfx();
    stateTens_t& getfxx();
    stateR_commandC_t &getfu();
    stateR_commandC_commandD_t& getfuu();
    stateR_stateC_commandD_t& getfxu();
    stateR_commandC_stateD_t& getfux();
    stateMatTab_t& getfxList();
    stateR_commandC_tab_t& getfuList();
};

#endif // DYNAMICMODEL_H
