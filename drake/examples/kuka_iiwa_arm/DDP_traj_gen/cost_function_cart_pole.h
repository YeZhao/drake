#pragma once

#ifndef COSTFUNCTIONCARTPOLE_H
#define COSTFUNCTIONCARTPOLE_H

#include "config.h"
#include <iostream>

#include <Eigen/Dense>

using namespace Eigen;

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

class CostFunctionCartPole
{
public:
    CostFunctionCartPole();
private:
protected:
	stateMat_t Q;
	stateMat_t Qf;
	commandMat_t R;

	stateVec_t QDiagElementVec;
	stateVec_t QfDiagElementVec;
	commandVec_t RDiagElementVec;

	stateVecTab_t cx_new;
	commandVecTab_t cu_new; 
	stateMatTab_t cxx_new; 
	commandR_stateC_tab_t cux_new; 
	commandMatTab_t cuu_new;
	double c_new;
    // attributes
public:
	stateMat_t& getQ();
	stateMat_t& getQf();
	commandMat_t& getR();
	stateVecTab_t& getcx();
	commandVecTab_t& getcu();
	stateMatTab_t& getcxx();
	commandR_stateC_tab_t& getcux();
	commandMatTab_t& getcuu();
	double& getc();

	unsigned int N;
private:

protected:
    // methods
public:
private:
protected:
    // accessors
public:

};

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

#endif // COSTFUNCTIONCARTPOLE_H