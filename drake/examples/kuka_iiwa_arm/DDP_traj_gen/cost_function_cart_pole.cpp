#include "cost_function_cart_pole.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {
	
CostFunctionCartPole::CostFunctionCartPole()
{
    Q << .1,0.0,0.0,0.0,
         0.0,0.1,0.0,0.0,
         0.0,0.0,0.1,0.0,
         0.0,0.0,0.0,0.1;
    Qf << 1000.0,0.0,0.0,0.0,
          0.0,1000,0.0,0.0,
          0.0,0.0,1000,0.0,
          0.0,0.0,0.0,1000;
    R << 0.01;

    N = TimeHorizon/TimeStep;
    cx_new.resize(N+1);
    cu_new.resize(N+1);
    cxx_new.resize(N+1);
    cux_new.resize(N+1);
    cuu_new.resize(N+1);
}

stateMat_t& CostFunctionCartPole::getQ()
{
    return Q;
}

stateMat_t& CostFunctionCartPole::getQf()
{
    return Qf;
}

commandMat_t& CostFunctionCartPole::getR()
{
    return R;
}

stateVecTab_t& CostFunctionCartPole::getcx()
{
    return cx_new;
}

commandVecTab_t& CostFunctionCartPole::getcu()
{
    return cu_new;
}

stateMatTab_t& CostFunctionCartPole::getcxx()
{
    return cxx_new;
}

commandR_stateC_tab_t& CostFunctionCartPole::getcux()
{
    return cux_new;
}

commandMatTab_t& CostFunctionCartPole::getcuu()
{
    return cuu_new;
}

double& CostFunctionCartPole::getc()
{
    return c_new;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake