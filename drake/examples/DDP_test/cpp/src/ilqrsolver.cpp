#include "ilqrsolver.h"

/* Debug */
#include <iostream>
using namespace std;
/* */

using namespace Eigen;
using Eigen::VectorXd;

ILQRSolver::ILQRSolver(DynamicModel& myDynamicModel, CostFunction& myCostFunction, bool fullDDP, bool QPBox)
{
    //TRACE("initialize dynamic model and cost function\n");
    dynamicModel = &myDynamicModel;
    costFunction = &myCostFunction;
    stateNb = myDynamicModel.getStateNb();
    commandNb = myDynamicModel.getCommandNb();
    enableQPBox = QPBox;
    enableFullDDP = fullDDP;

    if(enableQPBox) cout << "Box QP is enabled" << endl;
    else cout << "Box QP is disabled" << endl;

    if(enableFullDDP) cout << "Full DDP is enabled" << endl;
    else cout << "Full DDP is disabled" << endl;

    if(QPBox)
    {
        qp = new QProblemB(commandNb);
        Options myOptions;
        myOptions.printLevel = PL_LOW;
        myOptions.enableRegularisation = BT_TRUE;
        myOptions.initialStatusBounds = ST_INACTIVE;
        myOptions.numRefinementSteps = 1;
        myOptions.enableCholeskyRefactorisation = 1;
        qp->setOptions(myOptions);

        xOpt = new real_t[commandNb];
        lowerCommandBounds = myDynamicModel.getLowerCommandBounds();
        upperCommandBounds = myDynamicModel.getUpperCommandBounds();
    }

    //tOptSet Op = INIT_OPTSET;
}

void ILQRSolver::FirstInitSolver(stateVec_t& myxInit, stateVec_t& myxgoal, unsigned int& myT,
                       double& mydt, unsigned int& mymax_iter, double& mystopCrit, double& mytolFun, double& mytolGrad)
{
    // TODO: double check opt params
    xInit = myxInit; // remove myxgoal. Double check whether this makes sense.
    xgoal = myxgoal;
    T = myT;
    dt = mydt;
    stopCrit = mystopCrit;
    
    // Eigen::VectorXd default_alpha;
    // default_alpha.setZero();
    // default_alpha << 1.0, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    //TRACE("initialize option parameters\n");
    //Op = INIT_OPTSET;
    standard_parameters(&Op);
    Op.xInit = myxInit;
    Op.n_hor = T-1; // TODO: to be checked
    Op.tolFun = mytolFun;
    Op.tolGrad = mytolGrad;
    Op.max_iter = mymax_iter;

    xList.resize(T+1);
    uList.resize(T);
    updatedxList.resize(T+1);
    updateduList.resize(T);
    tmpxPtr.resize(T+1);
    tmpuPtr.resize(T);
    k.setZero();
    K.setZero();
    kList.resize(T);
    KList.resize(T);
    
    FList.resize(T+1);
    
    Vx.resize(T);
    Vxx.resize(T);
    dV.setZero();

    /*xList = new stateVec_t[myT+1];
    uList = new commandVec_t[myT];
    updatedxList = new stateVec_t[myT+1];
    updateduList = new commandVec_t[myT];
    k.setZero();
    K.setZero();
    kList = new commandVec_t[myT];
    KList = new commandR_stateC_t[myT];*/

    // parameters for line search [TODO: to be optimized]
    alphaList[0] = 1.0;
    alphaList[1] = 0.5012;
    alphaList[2] = 0.2512;
    alphaList[3] = 0.1259;
    alphaList[4] = 0.0631;
    alphaList[5] = 0.0316;
    alphaList[6] = 0.0158;
    alphaList[7] = 0.0079;
    alphaList[8] = 0.0040;
    alphaList[9] = 0.0020;
    alphaList[10] = 0.0020;
    alphaList[10] = 0.0010;

    alpha = 1.0;
}

void ILQRSolver::solveTrajectory()
{
    initializeTraj(&Op);

    int diverge, backPassDone, fwdPassDone, newDeriv;
    double dlambda= Op.dlambdaInit;

    Op.lambda= Op.lambdaInit;
    Op.w_pen_l= Op.w_pen_init_l;
    Op.w_pen_f= Op.w_pen_init_f;
    newDeriv = 1; // i.e., flgChange

    // TODO: update multipliers
    //update_multipliers(Op, 1);

    for(iter=0;iter<Op.max_iter;iter++)
    {
        //TRACE("STEP 1: differentiate dynamics and cost along new trajectory\n");
        if(newDeriv){
            int nargout = 7;//fx,fu,cx,cu,cxx,cxu,cuu
            dynamicModel->cart_pole_dyn_cst(nargout, dt, xList, uList, xgoal, FList, costFunction->getcx(), costFunction->getcu(), costFunction->getcxx(), costFunction->getcux(), costFunction->getcuu(), costFunction->getc());
            newDeriv = 0;
        }
        //TRACE("Finish STEP 1\n");

        // ====== STEP 2: backward pass, compute optimal control law and cost-to-go
        // backPassDone = 0;
        // TRACE(("Back pass:\n"));
        // while(!backPassDone) {
        //     if(back_pass(o)) {
        //         if(Op.debug_level>=1)
        //             TRACE(("Back pass failed.\n"));

        //         dlambda= max(dlambda * Op.lambdaFactor, Op.lambdaFactor);
        //         Op.lambda= max(Op.lambda * dlambda, Op.lambdaMin);
        //         if(Op.lambda > Op.lambdaMax)
        //             break;
        //     } else {
        //         backPassDone= 1;
        //         TRACE(("...done\n"));
        //     }
        // }
        
        // // check for termination due to small gradient
        // // TODO: add constraint tolerance check
        // if(o->g_norm < o->tolGrad && o->lambda < 1e-5) {
        //     dlambda= min(dlambda / o->lambdaFactor, 1.0/o->lambdaFactor);
        //     o->lambda= o->lambda * dlambda * (o->lambda > o->lambdaMin);
        //     if(o->debug_level>=1)
        //         TRACE(("\nSUCCESS: gradient norm < tolGrad\n"));
        //     break;
        // }


        backwardLoop();
        forwardLoop(&Op);
        cout << "iteration:  " << iter << endl;
        if(changeAmount<stopCrit)
        {
            break;
        }
        tmpxPtr = xList;
        tmpuPtr = uList;
        xList = updatedxList;
        updatedxList = tmpxPtr;
        uList = updateduList;
        updateduList = tmpuPtr;
    }
}

void ILQRSolver::initializeTraj(tOptSet *Op)
{
    xList[0] = Op->xInit;
    commandVec_t zeroCommand;
    zeroCommand.setZero();
    verbosity = Op->print;
    // (low priority) TODO: implement control limit selection
    // (low priority) TODO: initialize trace data structure

    for(unsigned int i=0;i<T;i++)
    {
        uList[i] = zeroCommand;
        xList[i+1] = dynamicModel->computeNextState(dt,xList[i],xgoal,zeroCommand);
    }
}

void ILQRSolver::standard_parameters(tOptSet *o) {
    //o->alpha= default_alpha;
    o->n_alpha= 11;
    o->tolFun= 1e-4;
    o->tolConstraint= 1e-7; // TODO: to be modified
    o->tolGrad= 1e-4;
    o->max_iter= 500;
    o->lambdaInit= 1;
    o->dlambdaInit= 1;
    o->lambdaFactor= 1.6;
    o->lambdaMax= 1e10;
    o->lambdaMin= 1e-6;
    o->regType= 1;
    o->zMin= 0.0;
    o->debug_level= 2; // TODO: to be modified
    o->w_pen_init_l= 1.0; // TODO: to be modified
    o->w_pen_init_f= 1.0; // TODO: to be modified
    o->w_pen_max_l= 1e100;//set to INF originally
    o->w_pen_max_f= 1e100;//set to INF originally
    o->w_pen_fact1= 4.0; // 4...10 Bertsekas p. 123
    o->w_pen_fact2= 1.0; 
    o->print = 2;
}

void ILQRSolver::initTrajectory(tOptSet *Op)
{
    xList[0] = Op->xInit;
    commandVec_t zeroCommand;
    zeroCommand.setZero();
    for(unsigned int i=0;i<T;i++)
    {
        uList[i] = zeroCommand;
        xList[i+1] = dynamicModel->computeNextState(dt,xList[i],xgoal,zeroCommand);
    }
}

void ILQRSolver::backwardLoop()
{
    costFunction->computeFinalCostDeriv(xList[T]);
    nextVx = costFunction->getlx();
    nextVxx = costFunction->getlxx();

    lambda = 0.0;
    completeBackwardFlag = 0;

    while(!completeBackwardFlag)
    {
        completeBackwardFlag = 1;
        if(Op.regType == 1)
            lambdaEye = lambda*stateMat_t::Zero();//[to be checked, change it to One()]
        else
            lambdaEye = lambda*stateMat_t::Identity();
        //std::cout << "lambdaEye: " << lambdaEye << std::endl;

        diverge = 0;
        
        Vx[T] = costFunction->getcx()[T];
        Vxx[T] = costFunction->getcxx()[T];

        for(int i=T-1;i>=0;i--)
        {
            x = xList[i];
            u = uList[i];

            dynamicModel->computeAllModelDeriv(dt,x,xgoal,u);
            costFunction->computeAllCostDeriv(x,u);

            Qx = costFunction->getcx()[i] + dynamicModel->getfx().transpose() * nextVx;
            Qu = costFunction->getcu()[i] + dynamicModel->getfu().transpose() * nextVx;
            Qxx = costFunction->getcxx()[i] + dynamicModel->getfx().transpose() * (nextVxx+lambdaEye) * dynamicModel->getfx();
            Quu = costFunction->getcuu()[i] + dynamicModel->getfu().transpose() * (nextVxx+lambdaEye) * dynamicModel->getfu();
            Qux = costFunction->getcux()[i] + dynamicModel->getfu().transpose() * (nextVxx+lambdaEye) * dynamicModel->getfx();

            Qx = costFunction->getlx() + dynamicModel->getfx().transpose() * nextVx;
            Qu = costFunction->getlu() + dynamicModel->getfu().transpose() * nextVx;
            Qxx = costFunction->getlxx() + dynamicModel->getfx().transpose() * (nextVxx+lambdaEye) * dynamicModel->getfx();
            Quu = costFunction->getluu() + dynamicModel->getfu().transpose() * (nextVxx+lambdaEye) * dynamicModel->getfu();
            Qux = costFunction->getlux() + dynamicModel->getfu().transpose() * (nextVxx+lambdaEye) * dynamicModel->getfx();

            //for debugging
            if (i==T-1){
                std::cout << "--------- new iteration ---------" << std::endl;
                std::cout << "costFunction->getlx(): " << costFunction->getlx() << std::endl;
                std::cout << "costFunction->getcx(): " << costFunction->getcx()[i] << std::endl;
                
                std::cout << "costFunction->getlu(): " << costFunction->getlu() <<std::endl;
                std::cout << "costFunction->getcu(): " << costFunction->getcu()[i] << std::endl;
                
                std::cout << "costFunction->getlxx(): " << costFunction->getlxx() << std::endl;
                std::cout << "costFunction->getcxx(): " << costFunction->getcxx()[i] << std::endl;
                
                std::cout << "costFunction->getluu(): " << costFunction->getluu() <<std::endl;
                std::cout << "costFunction->getcuu(): " << costFunction->getcuu()[i] << std::endl;
                
                std::cout << "costFunction->getlux(): " << costFunction->getlux() <<std::endl;
                std::cout << "costFunction->getcux(): " << costFunction->getcux()[i] << std::endl;
            }

            // Qx = costFunction->getcx()[i] + dynamicModel->getfxList()[i].transpose()*Vx[i+1];
            // Qu = costFunction->getcu()[i] + dynamicModel->getfuList()[i].transpose()*Vx[i+1];
            // Qxx = costFunction->getcxx()[i] + dynamicModel->getfxList()[i].transpose()*(Vxx[i+1]+lambdaEye)*dynamicModel->getfxList()[i];
            // Quu = costFunction->getcuu()[i] + dynamicModel->getfuList()[i].transpose()*(Vxx[i+1]+lambdaEye)*dynamicModel->getfuList()[i];
            // Qux = costFunction->getcux()[i] + dynamicModel->getfuList()[i].transpose()*(Vxx[i+1]+lambdaEye)*dynamicModel->getfxList()[i];

            if(Op.regType == 1)
                QuuF = Quu + lambda*commandMat_t::Identity();
            else
                QuuF = Quu;
            
            if(enableFullDDP)
            {
                //Qxx += dynamicModel->computeTensorContxx(nextVx);
                //Qux += dynamicModel->computeTensorContux(nextVx);
                //Quu += dynamicModel->computeTensorContuu(nextVx);

                Qxx += dynamicModel->computeTensorContxx(Vx[i+1]);
                Qux += dynamicModel->computeTensorContux(Vx[i+1]);
                Quu += dynamicModel->computeTensorContuu(Vx[i+1]);
                QuuF += dynamicModel->computeTensorContuu(Vx[i+1]);
            }

            QuuInv = Quu.inverse();

            if(!isQuudefinitePositive(Quu))
            {
                /*
                  To be Implemented : Regularization (is Quu definite positive ?)
                */
                if(lambda==0.0) lambda += 1e-4;
                else lambda *= 10;
                completeBackwardFlag = 0;
                break;
            }

            if(enableQPBox)
            {
                //TRACE("Use Box QP");
                nWSR = 10; //[to be checked]
                H = Quu;
                g = Qu;
                lb = lowerCommandBounds - u;
                ub = upperCommandBounds - u;
                qp->init(H.data(),g.data(),lb.data(),ub.data(),nWSR);
                qp->getPrimalSolution(xOpt);
                k = Map<commandVec_t>(xOpt);
                K = -QuuInv*Qux;
                for(unsigned int i_cmd=0;i_cmd<commandNb;i_cmd++)
                {
                    if((k[i_cmd] == lowerCommandBounds[i_cmd]) | (k[i_cmd] == upperCommandBounds[i_cmd]))
                    {
                        K.row(i_cmd).setZero();
                    }
                }
            }
            else
            {
                // Cholesky decomposition by using upper triangular matrix
                //TRACE("Use Cholesky decomposition");
                Eigen::LLT<MatrixXd> lltOfQuuF(QuuF);
                Eigen::MatrixXd L = lltOfQuuF.matrixU(); 
                //assume QuuF is positive definite
                //TODO: check the non-PD case

                k = - L.inverse()*L.transpose().inverse()*Qu;
                K = - L.inverse()*L.transpose().inverse()*Qux;

                //k = -QuuInv*Qu;
                //K = -QuuInv*Qux;
            }   

            /*nextVx = Qx - K.transpose()*Quu*k;
            nextVxx = Qxx - K.transpose()*Quu*K;*/
            nextVx = Qx + K.transpose()*Quu*k + K.transpose()*Qu + Qux.transpose()*k;
            nextVxx = Qxx + K.transpose()*Quu*K+ K.transpose()*Qux + Qux.transpose()*K;
            nextVxx = 0.5*(nextVxx + nextVxx.transpose());

            kList[i] = k;
            KList[i] = K;
        }
    }
}

void ILQRSolver::forwardLoop(tOptSet *Op)
{
    changeAmount = 0.0;
    updatedxList[0] = Op->xInit;
    //cout << "Op->lambdaInit: " << Op->lambdaInit << endl;
    // TODO: Line search to be implemented
    alpha = 1.0;
    for(unsigned int i=0;i<T;i++)
    {
        updateduList[i] = uList[i] + alpha*kList[i] + KList[i]*(updatedxList[i]-xList[i]);
        updatedxList[i+1] = dynamicModel->computeNextState(dt,updatedxList[i],xgoal,updateduList[i]);
        for(unsigned int j=0;j<commandNb;j++)
        {
            changeAmount += abs(uList[i](j,0) - updateduList[i](j,0));
        }
    }
}

ILQRSolver::traj ILQRSolver::getLastSolvedTrajectory()
{
    lastTraj.xList = updatedxList;
    for(int i=0;i<T+1;i++)lastTraj.xList[i] += xgoal;
    lastTraj.uList = updateduList;
    lastTraj.iter = iter;
    return lastTraj;
}

bool ILQRSolver::isQuudefinitePositive(const commandMat_t & Quu)
{
    /*
      Todo : check if Quu is definite positive
    */
    //Eigen::JacobiSVD<commandMat_t> svd_Quu (Quu, ComputeThinU | ComputeThinV);
    Eigen::VectorXcd singular_values = Quu.eigenvalues();

    for(long i = 0; i < Quu.cols(); ++i)
    {
        if (singular_values[i].real() < 0.)
        {
            std::cout << "not sdp" << std::endl;
            return false;
        }
    }
    return true;
}
