#ifndef CARTPOLE_H
#define CARTPOLE_H

#include "config.h"
#include "matrixUtil.h"

#include "cost_function_cart_pole.h"

#include <iostream>
#include "dynamicmodel.h"
#include <Eigen/Dense>

#ifndef DEBUG_CART_POLE
#define DEBUG_CART_POLE 1
#else
    #if PREFIX1(DEBUG_CART_POLE)==1
    #define DEBUG_CART_POLE 1
    #endif
#endif

#define TRACE_CART_POLE(x) do { if (DEBUG_CART_POLE) printf(x);} while (0)

using namespace Eigen;

class CartPole : public DynamicModel
{
public:
    CartPole(double& mydt, unsigned int& myT);
private:
protected:

    // attributes //
public:
private:
    double dt;
    unsigned int T;
    //static const unsigned int stateNb=4;
    //static const unsigned int commandNb=1;
public:
    static const double mc;
    static const double mp;
    static const double l;
    static const double g;

    static const double k;
    static const double R;
    static const double Jm;
    static const double Jl;
    static const double fvm;
    static const double Cf0;
    static const double a;

    //commandVec_t lowerCommandBounds;
    //commandVec_t upperCommandBounds;
    CostFunctionCartPole costFunction_cart_pole;
private:
    stateVec_t Xreal;
    stateMat_t Id;
    stateMat_t A;
    stateMat_t Ad;
    stateR_commandC_t B;
    stateR_commandC_t Bd;
    double A13atan;
    double A33atan;
    stateMat_t fxBase;
    stateR_commandC_t fuBase;

    stateMat_t QxxCont;
    commandMat_t QuuCont;
    commandR_stateC_t QuxCont;

    stateMat_half_t H;
    stateMat_half_t C;
    stateVec_half_t G;
    stateR_half_commandC_t Bu;
    
    stateVec_t Xdot1;
    stateVec_t Xdot2;
    stateVec_t Xdot3;
    stateVec_t Xdot4;
    
    stateMat_t A1;
    stateMat_t A2;
    stateMat_t A3;
    stateMat_t A4;
    stateMat_t IdentityMat;

    stateR_commandC_t B1;
    stateR_commandC_t B2;
    stateR_commandC_t B3;
    stateR_commandC_t B4;

    stateVec_t Xp;
    stateVec_t Xp1;
    stateVec_t Xp2;
    stateVec_t Xp3;
    stateVec_t Xp4;

    stateVec_t Xm;
    stateVec_t Xm1;
    stateVec_t Xm2;
    stateVec_t Xm3;
    stateVec_t Xm4;
    bool debugging_print;
protected:
    // methods
public:
    stateVec_t cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U);
    void cart_pole_dyn_cst(const int& nargout, const double& dt, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, commandR_stateC_tab_t& cux, commandMatTab_t& cuu, double& c);
    void cart_pole_dyn_cst_short(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr, const stateVec_t& xgoal, stateVec_t& xList_next, double& c);
    stateVec_t update(const int& nargout, const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B);
    void grad(const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B);
    void hessian(const double& dt, const stateVec_t& X, const commandVec_t& U);
    
private:
protected:
        // accessors //
public:

};

#endif // CARTPOLE_H
