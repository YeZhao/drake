#pragma once

#ifndef KUKAARM_H
#define KUKAARM_H

#include "config.h"
#include "cost_function_kuka_arm.h"

#include "drake/common/drake_path.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"

#include <cstdio>
#include <iostream>
#include <Eigen/Dense>
#include <math.h>

// #include <mutex>
// std::mutex mtx;

#define pi 3.141592653

#ifndef DEBUG_KUKA_ARM
#define DEBUG_KUKA_ARM 1
#else
    #if PREFIX1(DEBUG_KUKA_ARM)==1
    #define DEBUG_KUKA_ARM 1
    #endif
#endif

#define TRACE_KUKA_ARM(x) do { if (DEBUG_KUKA_ARM) printf(x);} while (0)

using namespace Eigen;
using namespace std;

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

class KukaArm
{
public:
    KukaArm();
    KukaArm(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal);
    ~KukaArm(){};
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
    struct timeprofile
    {
        double time_period1, time_period2, time_period3, time_period4;
        unsigned int counter0_, counter1_, counter2_;
    };

private:
    double dt;
    unsigned int N;
    bool initial_phase_flag_;
    struct timeprofile finalTimeProfile;
    struct timeval tbegin_period, tend_period, tbegin_period2, tend_period2, tbegin_period3, tend_period3, tbegin_period4, tend_period4;

public:
    static const double mc, mp, l, g;

    stateVec_t xgoal;
private:
    
    stateMat_half_t H, C;
    stateVec_half_t G;
    stateR_half_commandC_t Bu;
    stateVec_half_t velocity;
    stateVec_half_t accel;
    stateVec_t Xdot_new;
    stateVec_half_t vd;
    stateVecTab_half_t vd_thread;
    stateVecTab_t Xdot_new_thread;

    stateVec_t Xdot1, Xdot2, Xdot3, Xdot4;
    stateMat_t A1, A2, A3, A4, IdentityMat;
    stateR_commandC_t B1, B2, B3, B4;
    stateVec_t Xp, Xp1, Xp2, Xp3, Xp4, Xm, Xm1, Xm2, Xm3, Xm4;
    
    bool debugging_print;
    stateMat_t AA;
    stateR_commandC_t BB;
    stateMatTab_t A_temp;
    stateR_commandC_tab_t B_temp;
    
    std::unique_ptr<RigidBodyTree<double>> robot_thread_{nullptr};

    Eigen::VectorXd q;
    Eigen::VectorXd qd;
    std::vector<Eigen::VectorXd> q_thread, qd_thread;
protected:
    // methods
public:
    stateVec_t kuka_arm_dynamics(const stateVec_t& X, const commandVec_t& tau);
    stateVec_t kuka_arm_dynamicsThread1(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread2(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread3(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread4(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread5(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread6(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread7(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread8(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread9(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread10(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread11(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread12(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread13(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread14(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread15(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread16(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread17(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread18(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread19(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread20(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread21(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);
    stateVec_t kuka_arm_dynamicsThread22(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index);

    void kuka_arm_dyn_cst_ilqr(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunctionKukaArm*& costFunction);
    void kuka_arm_dyn_cst_min_output(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr,  const bool& isUNan, stateVec_t& xList_next, CostFunctionKukaArm*& costFunction);
    void kuka_arm_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, CostFunctionKukaArm*& costFunction);
    void kuka_arm_dyn_cst_v3(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList, CostFunctionKukaArm*& costFunction);
    stateVec_t update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    void grad(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B);
    void hessian(const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu);    
    struct timeprofile getFinalTimeProfile();

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

#endif // KUKAARM_H
