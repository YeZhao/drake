#include <iostream>
#include <fstream>

#include "config.h"
#include "matrixUtil.h"

#include "ilqrsolver.h"
#include "udpsolver.h"
#include "cart_pole.h"
#include "cost_function_cart_pole.h"

#include <time.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;

#define pi M_PI
#define useILQRSolver 0
#define useUDPSolver 1

int main()
{
    struct timeval tbegin,tend;
    double texec = 0.0;
    stateVec_t xinit,xgoal;

    xinit << 0.0,0.0,0.0,0.0;
    xgoal << 0.0,pi,0.0,0.0;

    double T = TimeHorizon;
    double dt = TimeStep;
    unsigned int N = (int)(T/dt);
    double tolFun = 1e-5;//relaxing default value: 1e-10;
    double tolGrad = 1e-5;//relaxing default value: 1e-10;
    unsigned int iterMax = 150;
    #if useILQRSolver
        ILQRSolver::traj lastTraj;
        CartPole cartPoleModel(dt, N);
        CostFunctionCartPole costCartPole;
        ILQRSolver testSolverCartPole(cartPoleModel,costCartPole,ENABLE_FULLDDP,ENABLE_QPBOX);
        testSolverCartPole.firstInitSolver(xinit, xgoal, N, dt, iterMax, tolFun, tolGrad);    
    #endif
    #if useUDPSolver    
        double scale = 0.01;
        UDPSolver::traj lastTraj;
        CartPole cartPoleModel(dt, N);
        CostFunctionCartPole costCartPole;
        UDPSolver testSolverCartPole(cartPoleModel,costCartPole,ENABLE_FULLDDP,ENABLE_QPBOX);
        testSolverCartPole.firstInitSolver(xinit, xgoal, N, dt, scale, iterMax, tolFun, tolGrad);    
    #endif

    // run one or multiple times and then average
    int Num_run = 1;
    gettimeofday(&tbegin,NULL);
    for(int i=0;i<Num_run;i++) testSolverCartPole.solveTrajectory();
    gettimeofday(&tend,NULL);

    lastTraj = testSolverCartPole.getLastSolvedTrajectory();

    texec=((double)(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
    texec /= Num_run;

    cout << endl;
    cout << "Number of iterations: " << lastTraj.iter << endl;
    cout << "Final cost: " << lastTraj.finalCost << endl;
    cout << "Final gradient: " << lastTraj.finalGrad << endl;
    cout << "Final lambda: " << lastTraj.finalLambda << endl;
    cout << "Execution time by time step (second): " << texec/N << endl;
    cout << "Execution time per iteration (second): " << texec/lastTraj.iter << endl;
    cout << "Total execution time of the solver (second): " << texec << endl;
    cout << "\tTime of derivative (second): " << lastTraj.time_derivative.sum() << " (" << 100.0*lastTraj.time_derivative.sum()/texec << "%)" << endl;
    cout << "\tTime of backward pass (second): " << lastTraj.time_backward.sum() << " (" << 100.0*lastTraj.time_backward.sum()/texec << "%)" << endl;
    cout << "\tTime of forward pass (second): " << lastTraj.time_forward.sum() << " (" << 100.0*lastTraj.time_forward.sum()/texec << "%)" << endl;
    
    ofstream file("results.csv",ios::out | ios::trunc);

    if(file)
    {
        file << "x,theta,xDot,thetaDot,u" << endl;
        for(int i=0;i<N;i++) file << lastTraj.xList[i](0,0) << "," << lastTraj.xList[i](1,0) << "," << lastTraj.xList[i](2,0) << "," << lastTraj.xList[i](3,0) << "," << lastTraj.uList[i](0,0) << endl;
        file << lastTraj.xList[N](0,0) << "," << lastTraj.xList[N](1,0) << "," << lastTraj.xList[N](2,0) << "," << lastTraj.xList[N](3,0) << "," << 0.0 << endl;
        file.close();
    }
    else
        cerr << "error in open file" << endl;

    return 0;
}