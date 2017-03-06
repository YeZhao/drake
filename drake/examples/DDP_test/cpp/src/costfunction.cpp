#include "costfunction.h"

CostFunction::CostFunction()
{
	N = 50;
	cx_new.resize(N+1);
    cu_new.resize(N+1);
    cxx_new.resize(N+1);
    cux_new.resize(N+1);
    cuu_new.resize(N+1);
}

stateMat_t& CostFunction::getQ()
{
    return Q;
}

stateMat_t& CostFunction::getQf()
{
    return Qf;
}

commandMat_t& CostFunction::getR()
{
    return R;
}

stateVecTab_t& CostFunction::getcx()
{
    return cx_new;
}

commandVecTab_t& CostFunction::getcu()
{
    return cu_new;
}

stateMatTab_t& CostFunction::getcxx()
{
    return cxx_new;
}

commandR_stateC_tab_t& CostFunction::getcux()
{
    return cux_new;
}

commandMatTab_t& CostFunction::getcuu()
{
    return cuu_new;
}

double& CostFunction::getc()
{
    return c_new;
}