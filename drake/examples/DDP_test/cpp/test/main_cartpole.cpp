#include <iostream>
#include <fstream>

#include "config.h"
#include "matrixUtil.h"

#include "ilqrsolver.h"
#include "cart_pole.h"
#include "cost_function_cart_pole.h"

#include <time.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;

#define pi M_PI

int main()
{
    struct timeval tbegin,tend;
    double texec = 0.0;
    stateVec_t xinit,xgoal;

    xinit << 0.0,0.0,0.0,0.0;
    xgoal << 0.0,pi,0.0,0.0;

    double T = 5;
    double dt = 0.1;
    unsigned int N = (int)T/dt;
    double stopCrit = 1e-5;
    double tolFun = 1e-10;
    double tolGrad = 1e-10;
    unsigned int iterMax = 150;

    stateVecTab_t xList;
    commandVecTab_t uList;
    ILQRSolver::traj lastTraj;
    
    CartPole cartPoleModel(dt, N);
    CostFunctionCartPole costCartPole;
    ILQRSolver testSolverCartPole(cartPoleModel,costCartPole,ENABLE_FULLDDP,ENABLE_QPBOX);

    testSolverCartPole.FirstInitSolver(xinit,xgoal,N,dt,iterMax,stopCrit, tolFun, tolGrad);

    // run multiple times and then average
    int Num_run = 1;
    gettimeofday(&tbegin,NULL);
    for(int i=0;i<Num_run;i++) testSolverCartPole.solveTrajectory();
    gettimeofday(&tend,NULL);

    lastTraj = testSolverCartPole.getLastSolvedTrajectory();
    xList = lastTraj.xList;
    uList = lastTraj.uList;
    unsigned int iter = lastTraj.iter;

    double finalCost = lastTraj.finalCost;
    double finalGrad = lastTraj.finalGrad;
    double finalLambda = lastTraj.finalLambda;
    Eigen::VectorXd time_backward, time_forward, time_derivative;
    time_backward = lastTraj.time_backward;
    //cout << "time_backward.size: " << time_backward.size() << endl;
    double time_backward_sum = time_backward.sum();
    time_forward = lastTraj.time_forward;
    double time_forward_sum = time_forward.sum();
    time_derivative = lastTraj.time_derivative;
    double time_derivative_sum = time_derivative.sum();

    texec=((double)(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
    texec /= Num_run;

    cout << endl;
    cout << "Final cost: " << finalCost << endl;
    cout << "Final gradient: " << finalGrad << endl;
    cout << "Final lambda: " << finalLambda << endl;

    cout << "Number of iterations: " << iter << endl;
    cout << "Execution time by time step (second): ";
    cout << texec/N << endl;
    cout << "Execution time per iteration (second): ";
    cout << texec/iter << endl;
    cout << "Total execution time of the solver (second): ";
    cout << texec << endl;
    cout << "\tTime of derivative (second): " << time_derivative_sum << " (" << 100.0*time_derivative_sum/texec << "%)" << endl;
    cout << "\tTime of backward pass (second): " << time_backward_sum << " (" << 100.0*time_backward_sum/texec << "%)" << endl;
    cout << "\tTime of forward pass (second): " << time_forward_sum << " (" << 100.0*time_forward_sum/texec << "%)" << endl;
    
    ofstream file("results.csv",ios::out | ios::trunc);

    if(file)
    {
        file << "tau,tauDot,q,qDot,u" << endl;
        for(int i=0;i<N;i++) file << xList[i](0,0) << "," << xList[i](1,0) << "," << xList[i](2,0) << "," << xList[i](3,0) << "," << uList[i](0,0) << endl;
        file << xList[N](0,0) << "," << xList[N](1,0) << "," << xList[N](2,0) << "," << xList[N](3,0) << "," << 0.0 << endl;
        file.close();
    }
    else
        cerr << "error in open file" << endl;

    return 0;

}
