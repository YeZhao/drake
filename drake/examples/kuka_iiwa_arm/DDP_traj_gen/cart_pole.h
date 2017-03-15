#pragma once

#ifndef CARTPOLE_H
#define CARTPOLE_H

#include "config.h"
#include "cost_function_cart_pole.h"

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/util/drakeGeometryUtil.h"

#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

#define pi 3.141592653

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

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

class CartPole
{
public:
    CartPole();
    explicit CartPole(double& mydt, unsigned int& myN, stateVec_t& myxgoal)
    {
        stateNb = 4;
        commandNb = 1;
        dt = mydt;
        N = myN;
        xgoal = myxgoal;
        fxList.resize(N);
        fuList.resize(N);

        fxxList.resize(stateSize);
        for(unsigned int i=0;i<stateSize;i++)
            fxxList[i].resize(N);
        fxuList.resize(commandSize);
        fuuList.resize(commandSize);
        for(unsigned int i=0;i<commandSize;i++){
            fxuList[i].resize(N);
            fuuList[i].resize(N);
        }

        fxx[0].setZero();
        fxx[1].setZero();
        fxx[2].setZero();
        fxx[3].setZero();
        fuu[0].setZero();
        fux[0].setZero();
        fxu[0].setZero();

        lowerCommandBounds << -50.0;
        upperCommandBounds << 50.0;

        H.setZero();
        C.setZero();
        G.setZero();
        Bu.setZero();
        velocity.setZero();
        accel.setZero();
        Xdot_new.setZero();

        A1.setZero();
        A2.setZero();
        A3.setZero();
        A4.setZero();
        B1.setZero();
        B2.setZero();
        B3.setZero();
        B4.setZero();
        IdentityMat.setIdentity();

        Xp1.setZero();
        Xp2.setZero();
        Xp3.setZero();
        Xp4.setZero();

        Xm1.setZero();
        Xm2.setZero();
        Xm3.setZero();
        Xm4.setZero();

        AA.setZero();
        BB.setZero();
        A_temp.resize(N);
        B_temp.resize(N);
        
        debugging_print = 0;
        
        initial_phase_flag_ = 1;
        q.resize(stateSize/2);
        qd.resize(stateSize/2);
    }

    ~CartPole(){};
private:
protected:
    // attributes
    unsigned int stateNb;
    unsigned int commandNb;
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
private:
    double dt;
    unsigned int N;
    bool initial_phase_flag_;
public:
    static const double mc, mp, l, g;

    struct traj_test
    {
        stateVecTab_t xList;
        commandVecTab_t uList;
        unsigned int iter;
        double finalCost;
        double finalGrad;
        double finalLambda;
        Eigen::VectorXd time_forward, time_backward, time_derivative, time_range1, time_range2;
    };
    stateVec_t xgoal;
private:
    
    stateMat_half_t H, C;
    stateVec_half_t G;
    stateR_half_commandC_t Bu;
    stateVec_half_t velocity;
    stateVec_half_t accel;
    stateVec_t Xdot_new;
    stateVec_half_t vd;

    stateVec_t Xdot1, Xdot2, Xdot3, Xdot4;
    stateMat_t A1, A2, A3, A4, IdentityMat;
    stateR_commandC_t B1, B2, B3, B4;
    stateVec_t Xp, Xp1, Xp2, Xp3, Xp4, Xm, Xm1, Xm2, Xm3, Xm4;
    
    bool debugging_print;
    stateMat_t AA;
    stateR_commandC_t BB;
    stateMatTab_t A_temp;
    stateR_commandC_tab_t B_temp;
    
    struct traj_test lastTraj;
    std::unique_ptr<RigidBodyTree<double>> tree_{nullptr};
    Eigen::VectorXd q;
    Eigen::VectorXd qd;
protected:
    // methods
public:
    stateVec_t cart_pole_dynamics(const stateVec_t& X, const commandVec_t& tau);
    void cart_pole_dyn_cst_ilqr(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunctionCartPole*& costFunction);
    void cart_pole_dyn_cst_min_output(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr,  const bool& isUNan, stateVec_t& xList_next, CostFunctionCartPole*& costFunction);
    void cart_pole_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunctionCartPole*& costFunction);
    void cart_pole_dyn_cst_v3(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList, CostFunctionCartPole*& costFunction);
    stateVec_t update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    void grad(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    void hessian(const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu);    

    unsigned int getStateNb();
    unsigned int getCommandNb();
    commandVec_t& getLowerCommandBounds();
    commandVec_t& getUpperCommandBounds();
    stateMatTab_t& getfxList();
    stateR_commandC_tab_t& getfuList();
private:
protected:
        // accessors //
public:
};

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

#endif // CARTPOLE_H