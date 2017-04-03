#pragma once

#ifndef UDPSOLVER_H
#define UDPSOLVER_H

#include "config.h"
#include "kuka_arm.h"
#include "cost_function_kuka_arm.h"
#include <string>

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>
//#include <qpOASES.hpp>
//#include <qpOASES/QProblemB.hpp>

#ifndef DEBUG_UDP
#define DEBUG_UDP 1
#else
    #if PREFIX1(DEBUG_UDP)==1
    #define DEBUG_UDP 1
    #endif
#endif

#define TRACE_UDP(x) do { if (DEBUG_UDP) printf(x);} while (0)

using namespace Eigen;
//USING_NAMESPACE_QPOASES

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {

class UDPSolver
{
public:
    struct traj
    {
        stateVecTab_t xList;
        commandVecTab_t uList;
        unsigned int iter;
        double finalCost;
        double finalGrad;
        double finalLambda;
        Eigen::VectorXd time_forward, time_backward, time_derivative, time_range1, time_range2, time_range3, time_iteration;
    };

    struct tOptSet {
        int n_hor;
        int debug_level;
        stateVec_t xInit;
        double new_cost, cost, dcost, lambda, dlambda, g_norm, expected;
        double **p;
        const double *alpha;
        int n_alpha;
        double lambdaMax;
        double lambdaMin;
        double lambdaInit;
        double dlambdaInit;
        double lambdaFactor;
        unsigned int max_iter;
        double tolGrad;
        double tolFun;
        double tolConstraint;
        double zMin;
        int regType;
        int iterations;
        int *log_linesearch;
        double *log_z;
        double *log_cost;
        double dV[2];
        
        double w_pen_l;
        double w_pen_f;
        double w_pen_max_l;
        double w_pen_max_f;
        double w_pen_init_l;
        double w_pen_init_f;
        double w_pen_fact1;
        double w_pen_fact2;
        
        int print;
        double print_head; // print headings every print_head lines
        double last_head;
        Eigen::VectorXd time_backward, time_forward, time_derivative, time_range1, time_range2, time_range3, time_iteration;
        Eigen::VectorXd alphaList;
        // traj_t *nominal;
        // traj_t *candidates[NUMBER_OF_THREADS]; 
        // traj_t trajectories[NUMBER_OF_THREADS+1];
        // multipliers_t multipliers;
    };

public:
    UDPSolver(KukaArm& iiwaDynamicModel, CostFunctionKukaArm& iiwaCostFunction,bool fullDDP=0,bool QPBox=0);
private:
protected:
    // attributes //
public:
private:
    KukaArm* dynamicModel;
    CostFunctionKukaArm* costFunction;
    unsigned int stateNb;
    unsigned int commandNb;
    stateVec_t xInit;
    stateVec_t xgoal;
    unsigned int N;
    unsigned int iter;
    double dt;
    double scale;

    stateVecTab_t xList;
    commandVecTab_t uList;
    commandVecTab_t uListFull;
    commandVec_t u_NAN;
    stateVecTab_t updatedxList;
    commandVecTab_t updateduList;
    stateVecTab_t FList;
    costVecTab_t costList;
    costVecTab_t costListNew;
    struct traj lastTraj;
    struct timeval tbegin_time_fwd, tend_time_fwd, tbegin_time_bwd, tend_time_bwd, tbegin_time_deriv, tend_time_deriv;
    struct timeval tbegin_test, tend_test, tbegin_test2, tend_test2, tbegin_test3, tend_test3, tbegin_iteration, tend_iteration;

    stateVecTab_t Vx;
    stateMatTab_t Vxx;

    stateVec_t Qx;
    stateMat_t Qxx;
    commandVec_t Qu;
    commandMat_t Quu;
    commandMat_t QuuF;
    commandMat_t QuuInv;
    commandR_stateC_t Qux;
    commandVec_t k;
    commandR_stateC_t K;
    commandVecTab_t kList;
    commandR_stateC_tab_t KList;
    double alpha;

    stateMat_t lambdaEye;
    unsigned int backPassDone, fwdPassDone, initFwdPassDone;
    unsigned int diverge;

    /* QP variables */
    //QProblemB* qp;
    bool enableQPBox, enableFullDDP;
    commandMat_t H;
    commandVec_t g;
    commandVec_t lowerCommandBounds, upperCommandBounds;
    commandVec_t lb, ub;
    int nWSR;
    //real_t* xOpt;

    tOptSet Op;
    Eigen::Vector2d dV;
    bool debugging_print;    
    int newDeriv;

    stateVecTab_t Xdot1, Xdot2, Xdot3, Xdot4, XnextThread;
    double g_norm_i, g_norm_max, g_norm_sum;

    /* matrix in doBackwardPass() */
    Eigen::MatrixXd augMatrix, Sig, augState, G, D, df, M, HH;
    stateAug_t QxQu, mu;
    commandR_stateC_t ZeroLowerLeftMatrix;
    stateR_commandC_t ZeroUpperRightMatrix;
    stateMat_t Vxx_next_inverse;
    commandMat_t cuu_inverse;
    stateVec_t X_new;
    bool isUNan;
    #if MULTI_THREAD
        std::vector<std::thread> thread;
    #endif
    bool enable_rk4_;
    bool enable_euler_;
protected:
    // methods
public:
    void firstInitSolver(stateVec_t& iiwaxInit, stateVec_t& iiwaxDes, unsigned int& iiwaN,
                    double& iiwadt, double& iiwascale, unsigned int& iiwamax_iter, double& iiwatolFun, double& iiwatolGrad);
    void solveTrajectory();
    void initializeTraj();
    void standardizeParameters(tOptSet *o);
    struct traj getLastSolvedTrajectory();
    void doBackwardPass();
    void doForwardPass();
    bool isPositiveDefinite(const commandMat_t & Quu);
    stateVec_t rungeKuttaStepBackward(stateAug_t augX, double& dt, unsigned int i);
    stateVec_t eulerStepBackward(stateAug_t augX, double& dt, unsigned int i);
    stateVec_t rungeKutta3StepBackward(stateAug_t augX, commandVec_t U_previous, double& dt, unsigned int i);
    void rungeKuttaStepBackwardThread(stateAug_t augX, double dt, unsigned int i);
    void rungeKuttaStepBackwardTwoSigmaPointsThread1(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread2(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread3(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread4(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread5(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread6(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread7(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread8(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread9(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread10(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread11(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    void rungeKuttaStepBackwardTwoSigmaPointsThread12(stateAug_t augXThread, stateAug_t augXThreadNext, stateAug_t augXThreadNextNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread13(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread14(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread15(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread16(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread17(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread18(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread19(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
    // void rungeKuttaStepBackwardTwoSigmaPointsThread20(stateAug_t augXThread, stateAug_t augXThreadNext, double dt, unsigned int iThread);
protected:
};

}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

#endif // UDPSOLVER_H
