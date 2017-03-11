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
using namespace std;

class CartPole : public DynamicModel
{
public:
    CartPole(double& mydt, unsigned int& myN, stateVec_t& myxgoal);
private:
protected:
    // attributes
public:
private:
    double dt;
    unsigned int N;
public:
    static const double mc, mp, l, g;
private:
    
    stateMat_half_t H, C;
    stateVec_half_t G;
    stateR_half_commandC_t Bu;
    stateVec_half_t velocity;
    stateVec_half_t accel;
    stateVec_t X_new;
    
    stateVec_t Xdot1, Xdot2, Xdot3, Xdot4;
    stateMat_t A1, A2, A3, A4, IdentityMat;
    stateR_commandC_t B1, B2, B3, B4;
    stateVec_t Xp, Xp1, Xp2, Xp3, Xp4, Xm, Xm1, Xm2, Xm3, Xm4;
    stateVec_t xgoal;
    bool debugging_print;
    stateMat_t AA;
    stateVec_t BB;
    stateMatTab_t A_temp;//dummy matrix
    stateR_commandC_tab_t B_temp;//dummy matrix
    
protected:
    // methods
public:
    stateVec_t cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U);
    void cart_pole_dyn_cst_ilqr(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunction*& costFunction);
    void cart_pole_dyn_cst_min_output(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr, stateVec_t& xList_next, CostFunction*& costFunction);
    void cart_pole_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunction*& costFunction);
    void cart_pole_dyn_cst_v3(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList, CostFunction*& costFunction);
    stateVec_t update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B);
    void grad(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B);
    void hessian(const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu);    
private:
protected:
        // accessors //
public:
};

#endif // CARTPOLE_H