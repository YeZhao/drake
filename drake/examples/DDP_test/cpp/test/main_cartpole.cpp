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

    xinit << -3.0,0.0,0.0,0.0;
    xgoal << 0.0,0.0,0.0,0.0;

    unsigned int T = 50;
    double dt=1e-4;

/*    xinit << 0.0,0.0,0.0,0.0;
    xgoal << 0.0,pi,0.0,0.0;

    unsigned int T = 50;
    double dt = 0.1;*/
    double stopCrit = 1e-5;
    double tolFun = 1e-10;
    double tolGrad = 1e-10;
    unsigned int iterMax = 150;

    stateVecTab_t xList;
    commandVecTab_t uList;
    ILQRSolver::traj lastTraj;
    
    CartPole cartPoleModel(dt);
    CostFunctionCartPole costCartPole;
    ILQRSolver testSolverCartPole(cartPoleModel,costCartPole,ENABLE_FULLDDP,ENABLE_QPBOX);

    testSolverCartPole.FirstInitSolver(xinit,xgoal,T,dt,iterMax,stopCrit, tolFun, tolGrad);

    // run multiple times and then average
    int N = 1;//100
    gettimeofday(&tbegin,NULL);
    for(int i=0;i<N;i++) testSolverCartPole.solveTrajectory();
    gettimeofday(&tend,NULL);

    lastTraj = testSolverCartPole.getLastSolvedTrajectory();
    xList = lastTraj.xList;
    uList = lastTraj.uList;
    unsigned int iter = lastTraj.iter;

    texec=((double)(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
    texec /= N;

    cout << endl;
    cout << "Total execution time of the solver ";
    cout << texec << endl;
    cout << "Execution time by time step ";
    cout << texec/T << endl;
    cout << "Number of iterations : " << iter << endl;

    ofstream file("results.csv",ios::out | ios::trunc);

    if(file)
    {
        file << "tau,tauDot,q,qDot,u" << endl;
        for(int i=0;i<T;i++) file << xList[i](0,0) << "," << xList[i](1,0) << "," << xList[i](2,0) << "," << xList[i](3,0) << "," << uList[i](0,0) << endl;
        file << xList[T](0,0) << "," << xList[T](1,0) << "," << xList[T](2,0) << "," << xList[T](3,0) << "," << 0.0 << endl;
        file.close();
    }
    else
        cerr << "error in open file" << endl;

    return 0;

}
