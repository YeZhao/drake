#include "costfunction.h"

CostFunction::CostFunction()
{
	T = 50;
	cx_new.resize(T+1);
    cu_new.resize(T+1);
    cxx_new.resize(T+1);
    cux_new.resize(T+1);
    cuu_new.resize(T+1);
}

stateVec_t& CostFunction::getlx()
{
    return lx;
}

stateMat_t& CostFunction::getlxx()
{
    return lxx;
}

commandVec_t& CostFunction::getlu()
{
    return lu;
}

commandMat_t& CostFunction::getluu()
{
    return luu;
}

commandR_stateC_t& CostFunction::getlux()
{
    return lux;
}

stateR_commandC_t& CostFunction::getlxu()
{
    return lxu;
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