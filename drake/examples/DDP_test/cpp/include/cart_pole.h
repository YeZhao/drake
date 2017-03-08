#ifndef CARTPOLE_H
#define CARTPOLE_H

#include "config.h"
#include "matrixUtil.h"
#include "cost_function_cart_pole.h"
#include "dynamicmodel.h"

#include <cstdio>
#include <iostream>
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
    CartPole(double& mydt, unsigned int& myN);
private:
protected:
    // attributes
public:
private:
    double dt;
    unsigned int N;
public:
    static const double mc, mp, l, g;

    CostFunctionCartPole costFunction_cart_pole;
    CostFunctionCartPole *costFunctionCartPole;
private:
    
    stateMat_half_t H, C;
    stateVec_half_t G;
    stateR_half_commandC_t Bu;
    
    stateVec_t Xdot1, Xdot2, Xdot3, Xdot4;
    stateMat_t A1, A2, A3, A4, IdentityMat;
    stateR_commandC_t B1, B2, B3, B4;
    stateVec_t Xp, Xp1, Xp2, Xp3, Xp4, Xm, Xm1, Xm2, Xm3, Xm4;
    bool debugging_print;
protected:
    // methods
public:
    stateVec_t cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U);
    void cart_pole_dyn_cst(const int& nargout, const double& dt, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, commandR_stateC_tab_t& cux, commandMatTab_t& cuu, double& c);
    void cart_pole_dyn_cst_short(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr, const stateVec_t& xgoal, stateVec_t& xList_next, double& c);
    void cart_pole_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, commandR_stateC_tab_t& cux, commandMatTab_t& cuu, double& c);
    void cart_pole_dyn_cst_v3(const int& nargout, const double& dt, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList, stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, commandR_stateC_tab_t& cux, commandMatTab_t& cuu, double& c);
    stateVec_t update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B);
    void grad(const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B);
    void hessian(const double& dt, const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu);
    
private:
protected:
        // accessors //
public:
};

#endif // CARTPOLE_H