#include "cost_function_cart_pole.h"

CostFunctionCartPole::CostFunctionCartPole()
{
    // Q << 100.0,0.0,0.0,0.0,
    //             0.0,0.0,0.0,0.0,
    //             0.0,0.0,0.0,0.0,
    //             0.0,0.0,0.0,0.0;
    // R << 0.1;

    // Qf = Q;

    Q << .1,0.0,0.0,0.0,
         0.0,0.1,0.0,0.0,
         0.0,0.0,0.1,0.0,
         0.0,0.0,0.0,0.1;
    Qf << 1000.0,0.0,0.0,0.0,
          0.0,1000,0.0,0.0,
          0.0,0.0,1000,0.0,
          0.0,0.0,0.0,1000;
    R << 0.01;

    lxx = Q;
    luu = R;
    lux << 0.0,0.0,0.0,0.0;
    lxu << 0.0,0.0,0.0,0.0;
    lx.setZero();
}

void CostFunctionCartPole::computeAllCostDeriv(const stateVec_t& X, const commandVec_t& U)
{
//    lx = Q*(X-Xdes);
    lx = Q*X;
    lu = R*U;
}

void CostFunctionCartPole::computeFinalCostDeriv(const stateVec_t& X)
{
    lx = Qf*X;
}