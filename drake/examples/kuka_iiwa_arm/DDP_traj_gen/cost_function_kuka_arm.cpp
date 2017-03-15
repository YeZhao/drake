#include "cost_function_kuka_arm.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {
	
CostFunctionKukaArm::CostFunctionKukaArm()
{
    
    // Q << .1,0.0,0.0,0.0,
    //      0.0,0.1,0.0,0.0,
    //      0.0,0.0,0.1,0.0,
    //      0.0,0.0,0.0,0.1;
    // Qf << 1000.0,0.0,0.0,0.0,
    //       0.0,1000,0.0,0.0,
    //       0.0,0.0,1000,0.0,
    //       0.0,0.0,0.0,1000;
    // R << 0.01;
    
    QDiagElementVec << 100, 100, 100, 100, 100, 100, 100,  10, 10, 10, 10, 10, 10, 10;
    QfDiagElementVec << 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0,    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0;
    RDiagElementVec << 0.005, 0.005, 0.007, 0.007, 0.02, 0.02, 0.05;

    Q = QDiagElementVec.asDiagonal();
    Qf = QfDiagElementVec.asDiagonal();
    R = RDiagElementVec.asDiagonal();

    N = TimeHorizon/TimeStep;
    cx_new.resize(N+1);
    cu_new.resize(N+1);
    cxx_new.resize(N+1);
    cux_new.resize(N+1);
    cuu_new.resize(N+1);
}

stateMat_t& CostFunctionKukaArm::getQ()
{
    return Q;
}

stateMat_t& CostFunctionKukaArm::getQf()
{
    return Qf;
}

commandMat_t& CostFunctionKukaArm::getR()
{
    return R;
}

stateVecTab_t& CostFunctionKukaArm::getcx()
{
    return cx_new;
}

commandVecTab_t& CostFunctionKukaArm::getcu()
{
    return cu_new;
}

stateMatTab_t& CostFunctionKukaArm::getcxx()
{
    return cxx_new;
}

commandR_stateC_tab_t& CostFunctionKukaArm::getcux()
{
    return cux_new;
}

commandMatTab_t& CostFunctionKukaArm::getcuu()
{
    return cuu_new;
}

double& CostFunctionKukaArm::getc()
{
    return c_new;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake
