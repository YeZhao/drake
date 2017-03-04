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
    Op.n_hor = T; // TODO: to be checked, it was set to T-1
    Op.tolFun = mytolFun;
    Op.tolGrad = mytolGrad;
    Op.max_iter = mymax_iter;

    xList.resize(T+1);
    uList.resize(T);
    updatedxList.resize(T+1);
    updateduList.resize(T);
    costList.resize(T+1);
    costListNew.resize(T+1);
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
    alphaList.resize(11);
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
    alphaList[10] = 0.0010;

    debugging_print = 0;
}

void ILQRSolver::solveTrajectory()
{
    initializeTraj();

    int diverge, backPassDone, fwdPassDone, newDeriv;
    
    //dlambda= Op.dlambdaInit;
    //lambda   = Op.lambda;
    
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

        //====== STEP 2: backward pass, compute optimal control law and cost-to-go
        backPassDone = 0;
        while(!backPassDone){
            backwardLoop();
            // if(diverge){//[Never entered, but entered here]
            //     if(verbosity > 2) printf("Cholesky failed at timestep %d.\n",diverge);
            //     Op.dlambda   = max(Op.dlambda * Op.lambdaFactor, Op.lambdaFactor);
            //     Op.lambda    = max(Op.lambda * Op.dlambda, Op.lambdaMin);
            //     if(Op.lambda > Op.lambdaMax) break;

            //         continue;
            // }
            backPassDone = 1;
        }

        // check for termination due to small gradient
        // TODO: add constraint tolerance check
        if(Op.g_norm < Op.tolGrad && Op.lambda < 1e-5){
            Op.dlambda= min(Op.dlambda / Op.lambdaFactor, 1.0/Op.lambdaFactor);
            Op.lambda= Op.lambda * Op.dlambda * (Op.lambda > Op.lambdaMin);
            if(Op.debug_level>=1){
                TRACE(("\nSUCCESS: gradient norm < tolGrad\n"));
            }
            break;
        }


        //====== STEP 3: line-search to find new control sequence, trajectory, cost
        fwdPassDone = 0;
        if(backPassDone){
            //only implement serial backtracking line-search
            for(int alpha_index = 0; alpha_index < alphaList.size(); alpha_index++){
                alpha = alphaList[alpha_index];
                forwardLoop();
                Op.dcost = accumulate(costListNew.begin(), costListNew.end(), 0.0) - accumulate(costList.begin(), costList.end(), 0.0);
                Op.expected = -alpha*(dV(0) + alpha*dV(1));
                double z;
                if(Op.expected > 0) {
                    z = Op.dcost/Op.expected;
                }else {
                    z = (double)(-signbit(Op.dcost));//[TODO:doublecheck]
                    TRACE("non-positive expected reduction: should not occur");//warning
                }
                if(z > Op.zMin){
                    fwdPassDone = 1;
                    break;
                }
            }
            if(!fwdPassDone) alpha = sqrt(-1.0);    
        }
        

        //====== STEP 4: accept step (or not), draw graphics, print status
        //cout << "iteration:  " << iter << endl;

        if (verbosity > 1 && Op.last_head == Op.print_head){
            Op.last_head = 0;
            TRACE("iteration,\t cost, \t reduction, \t expected, \t gradient, \t log10(lambda) \n");
            //printf("%-14d%-12.6g%-12.3g%-15.3g%-19.3g%-17.1f\n", iter+1, accumulate(costList.begin(), costList.end(), 0.0), Op.dcost, Op.expected, Op.g_norm, log10(Op.lambda));
            //printf("iter: %-3d  REJECTED    expected: %-11.3g    actual: %-11.3g    log10lam: %3.1f w_pen_l: %-9.3g w_pen_l: %-9.3g\n", iter+1, Op.expected , Op.dcost, Op.lambda, Op.w_pen_l, Op.w_pen_f);
        }
        
        if(fwdPassDone){
            // print status
            if (verbosity > 1){
                printf("%-14d%-12.6g%-15.3g%-15.3g%-19.3g%-17.1f\n", iter+1, accumulate(costList.begin(), costList.end(), 0.0), Op.dcost, Op.expected, Op.g_norm, log10(Op.lambda));
                Op.last_head = Op.last_head+1;
            }

            // decrease lambda
            Op.dlambda = min(Op.dlambda / Op.lambdaFactor, 1/Op.lambdaFactor);
            Op.lambda = Op.lambda * Op.dlambda * (Op.lambda > Op.lambdaMin);
            
            // accept changes
            tmpxPtr = xList;
            tmpuPtr = uList;
            xList = updatedxList;
            //updatedxList = tmpxPtr;
            uList = updateduList;
            //updateduList = tmpuPtr;
            costList = costListNew;
            newDeriv = 1;

            // terminate ?
            // TODO: add constraint tolerance check
            if(Op.dcost < Op.tolFun) {
                if(verbosity >= 1)
                    TRACE(("\nSUCCESS: cost change < tolFun\n"));
            
                break;
            }
        }else { // no cost improvement
            // increase lambda
            Op.dlambda= max(Op.dlambda * Op.lambdaFactor, Op.lambdaFactor);
            Op.lambda= max(Op.lambda * Op.dlambda, Op.lambdaMin);

            // if(o->w_pen_fact2>1.0) {
            //     o->w_pen_l= min(o->w_pen_max_l, o->w_pen_l*o->w_pen_fact2);
            //     o->w_pen_f= min(o->w_pen_max_f, o->w_pen_f*o->w_pen_fact2);
            //     forward_pass(o->nominal, o, 0.0, &o->cost, 1);
            // }
            
            // print status
            if(verbosity >= 1){
                printf("%-14d%-12.6s%-15.3g%-15.3g%-19.3g%-17.1f\n", iter+1, "No STEP", Op.dcost, Op.expected, Op.g_norm, log10(Op.lambda));
                Op.last_head = Op.last_head+1;
            }
                
            // terminate ?
            if(Op.lambda > Op.lambdaMax) {
                if(verbosity >= 1)
                    TRACE(("\nEXIT: lambda > lambdaMax\n"));
                break;
            }
        }

        // if(changeAmount<stopCrit)
        // {
        //     break;
        // }
        // tmpxPtr = xList;
        // tmpuPtr = uList;
        // xList = updatedxList;
        // updatedxList = tmpxPtr;
        // uList = updateduList;
        // updateduList = tmpuPtr;
    }

    Op.iterations = iter;

    if(!backPassDone) {
        if(verbosity >= 1)
            TRACE(("\nEXIT: no descent direction found.\n"));
        
        return;    
    } else if(iter >= Op.max_iter) {
        if(verbosity >= 1)
            TRACE(("\nEXIT: Maximum iterations reached.\n"));
        
        return;
    }
}

void ILQRSolver::initializeTraj()
{
    xList[0] = Op.xInit;
    commandVec_t zeroCommand;
    zeroCommand.setZero();
    verbosity = Op.print;
    // (low priority) TODO: implement control limit selection
    // (low priority) TODO: initialize trace data structure
    Op.print_head = 6;
    Op.last_head = Op.print_head;

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

void ILQRSolver::initTrajectory()
{
    xList[0] = Op.xInit;
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

    //lambda = 0.0;
    
    if(Op.regType == 1)
        lambdaEye = Op.lambda*stateMat_t::Zero();//[to be checked, change it to One()]
    else
        lambdaEye = Op.lambda*stateMat_t::Identity();
    //std::cout << "lambdaEye: " << lambdaEye << std::endl;

    diverge = 0;
    double g_norm_i, g_norm_max, g_norm_sum;
    
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
        if(debugging_print){
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
        }
        
        // Qx = costFunction->getcx()[i] + dynamicModel->getfxList()[i].transpose()*Vx[i+1];
        // Qu = costFunction->getcu()[i] + dynamicModel->getfuList()[i].transpose()*Vx[i+1];
        // Qxx = costFunction->getcxx()[i] + dynamicModel->getfxList()[i].transpose()*(Vxx[i+1]+lambdaEye)*dynamicModel->getfxList()[i];
        // Quu = costFunction->getcuu()[i] + dynamicModel->getfuList()[i].transpose()*(Vxx[i+1]+lambdaEye)*dynamicModel->getfuList()[i];
        // Qux = costFunction->getcux()[i] + dynamicModel->getfuList()[i].transpose()*(Vxx[i+1]+lambdaEye)*dynamicModel->getfxList()[i];

        if(Op.regType == 1)
            QuuF = Quu + Op.lambda*commandMat_t::Identity();
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
            if(Op.lambda==0.0) Op.lambda += 1e-4;
            else Op.lambda *= 10;
            backPassDone = 0;
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

            //A temporary solution: check the non-PD case
            if(lltOfQuuF.info() == Eigen::NumericalIssue)
                {
                    diverge = i;
                    TRACE("Possibly non semi-positive definitie matrix!");
                    return;
                }

            k = - L.inverse()*L.transpose().inverse()*Qu;
            K = - L.inverse()*L.transpose().inverse()*Qux;

            // old version
            //k = -QuuInv*Qu;
            //K = -QuuInv*Qux;
        }   

        // old version
        /*nextVx = Qx - K.transpose()*Quu*k;
        nextVxx = Qxx - K.transpose()*Quu*K;*/
        nextVx = Qx + K.transpose()*Quu*k + K.transpose()*Qu + Qux.transpose()*k;
        nextVxx = Qxx + K.transpose()*Quu*K+ K.transpose()*Qux + Qux.transpose()*K;
        nextVxx = 0.5*(nextVxx + nextVxx.transpose());

        //update cost-to-go approximation
        dV(0) += k.transpose()*Qu;
        commandMat_t c_mat_to_scalar;
        c_mat_to_scalar = 0.5*k.transpose()*Quu*k;
        dV(1) += c_mat_to_scalar(0,0);
        Vx[i] = Qx + K.transpose()*Quu*k + K.transpose()*Qu + Qux.transpose()*k;
        Vxx[i] = Qxx + K.transpose()*Quu*K+ K.transpose()*Qux + Qux.transpose()*K;
        Vxx[i] = 0.5*(Vxx[i] + Vxx[i].transpose());

        kList[i] = k;
        KList[i] = K;

        g_norm_max= 0.0;
        for(unsigned int j= 0; j<commandSize; j++) {
            g_norm_i= fabs(kList[i](j,0)) / (fabs(uList[i](j,0))+1.0);
            if(g_norm_i>g_norm_max) g_norm_max= g_norm_i;
        }
        g_norm_sum+= g_norm_max;
    }
    Op.g_norm= g_norm_sum/((double)(Op.n_hor));

    //TODO: handle Cholesky failure case
    // while(!backPassDone) {
    //     if(back_pass(o)) {
    //         if(Op.debug_level>=1)
    //             TRACE(("Back pass failed.\n"));

    //         Op.dlambda= max(Op.dlambda * Op.lambdaFactor, Op.lambdaFactor);
    //         Op.lambda= max(Op.lambda * Op.dlambda, Op.lambdaMin);
    //         if(Op.lambda > Op.lambdaMax)
    //             break;
    //     } else {
    //         backPassDone= 1;
    //         TRACE(("...done\n"));
    //     }
    // }
}

void ILQRSolver::forwardLoop()
{
    changeAmount = 0.0;
    updatedxList[0] = Op.xInit;
    //cout << "Op->lambdaInit: " << Op->lambdaInit << endl;
    // TODO: Line search to be implemented
    int nargout = 2;

    for(unsigned int i=0;i<T;i++)
    {
        updateduList[i] = uList[i] + alpha*kList[i] + KList[i]*(updatedxList[i]-xList[i]);
        //updatedxList[i+1] = dynamicModel->computeNextState(dt,updatedxList[i],xgoal,updateduList[i]); // [ToBeCommented]
        dynamicModel->cart_pole_dyn_cst_short(nargout, dt, updatedxList[i], updateduList[i], xgoal, updatedxList[i+1], costListNew[i]);// [ToBeUnCommented]
        //dynamicModel->cart_pole_dyn_cst(nargout, dt, updatedxList, updateduList, xgoal, FList, costFunction->getcx(), costFunction->getcu(), costFunction->getcxx(), costFunction->getcux(), costFunction->getcuu(), costFunction->getc());
        // [ToBeUnCommented]
        for(unsigned int j=0;j<commandNb;j++)
        {
            changeAmount += abs(uList[i](j,0) - updateduList[i](j,0));
        }
    }
    // [ToBeUnCommented]
    stateVec_t x_unused;
    x_unused.setZero();
    commandVec_t u_NAN;
    u_NAN << sqrt(-1.0);
    dynamicModel->cart_pole_dyn_cst_short(nargout, dt, updatedxList[T], u_NAN, xgoal, x_unused, costListNew[T]);
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
