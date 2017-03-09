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

    if(enableQPBox) TRACE("Box QP is enabled");
    else TRACE("Box QP is disabled");

    if(enableFullDDP) TRACE("Full DDP is enabled");
    else TRACE("Full DDP is disabled");

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

void ILQRSolver::firstInitSolver(stateVec_t& myxInit, stateVec_t& myxgoal, unsigned int& myN,
                       double& mydt, unsigned int& mymax_iter, double& mytolFun, double& mytolGrad)
{
    // TODO: double check opt params
    xInit = myxInit; // removed myxgoal. Double check whether this makes sense.
    xgoal = myxgoal;
    N = myN;
    dt = mydt;
    
    // Eigen::VectorXd default_alpha;
    // default_alpha.setZero();
    // default_alpha << 1.0, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    //TRACE("initialize option parameters\n");
    //Op = INIT_OPTSET;
    standardizeParameters(&Op);
    Op.xInit = myxInit;
    Op.n_hor = N;
    Op.tolFun = mytolFun;
    Op.tolGrad = mytolGrad;
    Op.max_iter = mymax_iter;
    Op.time_backward.resize(Op.max_iter);
    Op.time_backward.setZero();
    Op.time_forward.resize(Op.max_iter);
    Op.time_forward.setZero();
    Op.time_derivative.resize(Op.max_iter);
    Op.time_derivative.setZero();

    xList.resize(N+1);
    uList.resize(N);
    updatedxList.resize(N+1);
    updateduList.resize(N);
    costList.resize(N+1);
    costListNew.resize(N+1);
    k.setZero();
    K.setZero();
    kList.resize(N);
    KList.resize(N);
    
    FList.resize(N+1);
    
    Vx.resize(N);
    Vxx.resize(N);
    dV.setZero();

    // parameters for line search
    Op.alphaList.resize(11);
    Op.alphaList << 1.0, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    debugging_print = 0;
}

void ILQRSolver::solveTrajectory()
{
    initializeTraj();
    
    Op.lambda = Op.lambdaInit;
    Op.dlambda = Op.dlambdaInit;
    
    // TODO: update multipliers
    //update_multipliers(Op, 1);

    for(iter=0;iter<Op.max_iter;iter++)
    {
        //TRACE("STEP 1: differentiate dynamics and cost along new trajectory\n");
        if(newDeriv){
            int nargout = 7;//fx,fu,cx,cu,cxx,cxu,cuu
            commandVecTab_t uListFull;
            uListFull.resize(N+1);
            commandVec_t u_NAN;
            for(unsigned int i=0;i<u_NAN.size();i++)
                u_NAN(i,0) = sqrt(-1.0);
            
            for(unsigned int i=0;i<uList.size();i++)
                uListFull[i] = uList[i];
            uListFull[uList.size()] = u_NAN;

            // cout << "xList[0]: " << xList[0] << endl;
            // cout << "xList[10]: " << xList[10] << endl;
            // cout << "uListFull[0]: " << uListFull[0] << endl;
            // cout << "uListFull[10]: " << uListFull[10] << endl;

            gettimeofday(&tbegin_time_deriv,NULL);
            dynamicModel->cart_pole_dyn_cst_ilqr(nargout, xList, uListFull, FList, costFunction);
            //dynamicModel->cart_pole_dyn_cst(nargout, dt, xList, uListFull, xgoal, FList, costFunction->getcx(), costFunction->getcu(), costFunction->getcxx(), costFunction->getcux(), costFunction->getcuu(), costFunction->getc());
            gettimeofday(&tend_time_deriv,NULL);
            Op.time_derivative(iter) = ((double)(1000*(tend_time_deriv.tv_sec-tbegin_time_deriv.tv_sec)+((tend_time_deriv.tv_usec-tbegin_time_deriv.tv_usec)/1000)))/1000.0;
            //cout << "Op.time_derivative(iter): " << Op.time_derivative(iter) << endl;

            newDeriv = 0;
            // cout << "(initial position0) costFunction->getcx()[N]: " << costFunction->getcx()[N] << endl;
            // cout << "(initial position0) costFunction->getcx()[N-1]: " << costFunction->getcx()[N-1] << endl;
        }
        //TRACE("Finish STEP 1\n");

        //====== STEP 2: backward pass, compute optimal control law and cost-to-go
        backPassDone = 0;
        while(!backPassDone){

            gettimeofday(&tbegin_time_bwd,NULL);
            doBackwardPass();
            gettimeofday(&tend_time_bwd,NULL);
            Op.time_backward(iter) = ((double)(1000*(tend_time_bwd.tv_sec-tbegin_time_bwd.tv_sec)+((tend_time_bwd.tv_usec-tbegin_time_bwd.tv_usec)/1000)))/1000.0;
            //cout << "Op.time_backward(iter): " << Op.time_backward(iter) << endl;

            //TRACE("handle Cholesky failure case");
            if(diverge){
                if(Op.debug_level > 2) printf("Cholesky failed at timestep %d.\n",diverge);
                Op.dlambda   = max(Op.dlambda * Op.lambdaFactor, Op.lambdaFactor);
                Op.lambda    = max(Op.lambda * Op.dlambda, Op.lambdaMin);
                if(Op.lambda > Op.lambdaMax) break;

                    continue;
            }
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
            gettimeofday(&tbegin_time_fwd,NULL);
            //only implement serial backtracking line-search
            for(int alpha_index = 0; alpha_index < Op.alphaList.size(); alpha_index++){
                alpha = Op.alphaList[alpha_index];
                doForwardPass();
                Op.dcost = accumulate(costList.begin(), costList.end(), 0.0) - accumulate(costListNew.begin(), costListNew.end(), 0.0);
                Op.expected = -alpha*(dV(0) + alpha*dV(1));
                // std::cout << "costListNew[2]: " << costListNew[2] << std::endl;
                // std::cout << "costListNew[20]: " << costListNew[20] << std::endl;
                // std::cout << "costListNew[50]: " << costListNew[50] << std::endl;
                // std::cout << "costList[2]: " << costList[2] << std::endl;
                // std::cout << "costList[20]: " << costList[20] << std::endl;
                // std::cout << "costList[50]: " << costList[50] << std::endl;
                //std::cout << "accumulate(costListNew): " << accumulate(costListNew.begin(), costListNew.end(), 0.0) << std::endl;
                //std::cout << "accumulate(costList): " << accumulate(costList.begin(), costList.end(), 0.0) << std::endl;
                //std::cout << "Op.dcost: " << Op.dcost << std::endl;
                //std::cout << "Op.expected: " << Op.expected << std::endl;
                //std::cout << "alpha: " << alpha << std::endl;
                //std::cout << "dV: " << dV << std::endl;
                double z;
                if(Op.expected > 0) {
                    z = Op.dcost/Op.expected;
                }else {
                    z = (double)(-signbit(Op.dcost));//[TODO:doublecheck]
                    TRACE("non-positive expected reduction: should not occur \n");//warning
                }
                if(z > Op.zMin){
                    fwdPassDone = 1;
                    break;
                }
            }
            if(!fwdPassDone) alpha = sqrt(-1.0);
            gettimeofday(&tend_time_fwd,NULL);
            Op.time_forward(iter) = ((double)(1000*(tend_time_fwd.tv_sec-tbegin_time_fwd.tv_sec)+((tend_time_fwd.tv_usec-tbegin_time_fwd.tv_usec)/1000)))/1000.0;
    
        }
        
        //====== STEP 4: accept step (or not), draw graphics, print status
        if (Op.debug_level > 1 && Op.last_head == Op.print_head){
            Op.last_head = 0;
            TRACE("iteration,\t cost, \t reduction, \t expected, \t gradient, \t log10(lambda) \n");
        }
        
        if(fwdPassDone){
            // print status
            if (Op.debug_level > 1){
                if(!debugging_print) printf("%-14d%-12.6g%-15.3g%-15.3g%-19.3g%-17.1f\n", iter+1, accumulate(costList.begin(), costList.end(), 0.0), Op.dcost, Op.expected, Op.g_norm, log10(Op.lambda));
                Op.last_head = Op.last_head+1;
            }

            Op.dlambda = min(Op.dlambda / Op.lambdaFactor, 1.0/Op.lambdaFactor);
            Op.lambda = Op.lambda * Op.dlambda * (Op.lambda > Op.lambdaMin);

            // accept changes
            xList = updatedxList;
            uList = updateduList;
            costList = costListNew;
            newDeriv = 1;

            // terminate ?
            // TODO: add constraint tolerance check
            if(Op.dcost < Op.tolFun) {
                if(Op.debug_level >= 1)
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
            if(Op.debug_level >= 1){
                if(!debugging_print) printf("%-14d%-12.9s%-15.3g%-15.3g%-19.3g%-17.1f\n", iter+1, "No STEP", Op.dcost, Op.expected, Op.g_norm, log10(Op.lambda));
                Op.last_head = Op.last_head+1;
            }

            // terminate ?
            if(Op.lambda > Op.lambdaMax) {
                if(Op.debug_level >= 1)
                    TRACE(("\nEXIT: lambda > lambdaMax\n"));
                break;
            }
        }
    }

    Op.iterations = iter;

    if(!backPassDone) {
        if(Op.debug_level >= 1)
            TRACE(("\nEXIT: no descent direction found.\n"));
        
        return;    
    } else if(iter >= Op.max_iter) {
        if(Op.debug_level >= 1)
            TRACE(("\nEXIT: Maximum iterations reached.\n"));
        
        return;
    }
}

void ILQRSolver::initializeTraj()
{
    xList[0] = Op.xInit;
    commandVec_t zeroCommand;
    zeroCommand.setZero();
    // (low priority) TODO: implement control limit selection
    // (low priority) TODO: initialize trace data structure
    
    initFwdPassDone = 0;
    diverge = 1;
    for(int alpha_index = 0; alpha_index < Op.alphaList.size(); alpha_index++){
        alpha = Op.alphaList[alpha_index];    
        for(unsigned int i=0;i<N;i++)
        {
            uList[i] = zeroCommand;
        }
        doForwardPass();
        //simplistic divergence test
        int diverge_element_flag = 0;
        for(unsigned int i = 0; i < xList.size(); i++){
            for(unsigned int j = 0; j < xList[i].size(); j++){
                if(fabs(xList[i](j,0)) > 1e8)
                    diverge_element_flag = 1;
            }
        }
        if(!diverge_element_flag){
            diverge = 0;
            break;
        }
    }
    
    initFwdPassDone = 1;

    // std::cout << "updatedxList: " << std::endl;
    // for(unsigned int i = 0; i < updatedxList.size(); i++)
    //     std::cout <<  updatedxList[i] << std::endl;
    // std::cout << "updateduList: " << std::endl;
    // for(unsigned int i = 0; i < updatedxList.size(); i++)
    //     std::cout <<  updateduList[i] << std::endl;
    // std::cout << "costListNew: " << std::endl;
    // for(unsigned int i = 0; i < updatedxList.size(); i++)
    //     std::cout <<  costListNew[i] << std::endl;
    
    //constants, timers, counters
    newDeriv = 1; //flgChange
    Op.lambda= Op.lambdaInit;
    Op.w_pen_l= Op.w_pen_init_l;
    Op.w_pen_f= Op.w_pen_init_f;
    Op.dcost = 0;
    Op.expected = 0;
    Op.print_head = 6;
    Op.last_head = Op.print_head;
    if(Op.debug_level > 0) TRACE("\n=========== begin iLQG ===========\n");
}

void ILQRSolver::standardizeParameters(tOptSet *o) {
    //o->alpha= default_alpha;
    o->n_alpha = 11;
    o->tolFun = 1e-4;
    o->tolConstraint = 1e-7; // TODO: to be modified
    o->tolGrad = 1e-4;
    o->max_iter = 500;
    o->lambdaInit = 1;
    o->dlambdaInit = 1;
    o->lambdaFactor = 1.6;
    o->lambdaMax = 1e10;
    o->lambdaMin = 1e-6;
    o->regType = 1;
    o->zMin = 0.0;
    o->debug_level = 2; // == verbosity in matlab code
    o->w_pen_init_l = 1.0; // TODO: to be modified
    o->w_pen_init_f = 1.0; // TODO: to be modified
    o->w_pen_max_l = 1e100;//set to INF originally
    o->w_pen_max_f = 1e100;//set to INF originally
    o->w_pen_fact1 = 4.0; // 4...10 Bertsekas p. 123
    o->w_pen_fact2 = 1.0; 
    o->print = 2;
}

void ILQRSolver::doBackwardPass()
{    
    if(Op.regType == 1)
        lambdaEye = Op.lambda*stateMat_t::Identity();
    else
        lambdaEye = Op.lambda*stateMat_t::Zero();

    diverge = 0;
    double g_norm_i, g_norm_max, g_norm_sum;
    
    //cout << "(initial position) costFunction->getcx()[N]: " << costFunction->getcx()[N] << endl;
    Vx[N] = costFunction->getcx()[N];
    Vxx[N] = costFunction->getcxx()[N];
    dV.setZero();

    for(int i=N-1;i>=0;i--)
    {
        //for debugging
        if(debugging_print){
            if (i==N-1){
                std::cout << "--------- new iteration ---------" << std::endl;
                std::cout << "costFunction->getcx(): " << costFunction->getcx()[i] << std::endl;                
                std::cout << "costFunction->getcu(): " << costFunction->getcu()[i] << std::endl;                
                std::cout << "costFunction->getcxx(): " << costFunction->getcxx()[i] << std::endl;                
                std::cout << "costFunction->getcuu(): " << costFunction->getcuu()[i] << std::endl;                
                std::cout << "costFunction->getcux(): " << costFunction->getcux()[i] << std::endl;
            }    
        }
        
        Qx = costFunction->getcx()[i] + dynamicModel->getfxList()[i].transpose()*Vx[i+1];
        Qu = costFunction->getcu()[i] + dynamicModel->getfuList()[i].transpose()*Vx[i+1];
        Qxx = costFunction->getcxx()[i] + dynamicModel->getfxList()[i].transpose()*(Vxx[i+1])*dynamicModel->getfxList()[i];
        Quu = costFunction->getcuu()[i] + dynamicModel->getfuList()[i].transpose()*(Vxx[i+1])*dynamicModel->getfuList()[i];
        Qux = costFunction->getcux()[i] + dynamicModel->getfuList()[i].transpose()*(Vxx[i+1])*dynamicModel->getfxList()[i];

        // if(i > N-4){
        //     cout << "index i: " << i << endl;
        //     cout << "costFunction->getcx()[i]: " << costFunction->getcx()[i] << endl;
        //     cout << "costFunction->getcu()[i]: " << costFunction->getcu()[i] << endl;
        //     cout << "dynamicModel->getfxList()[i]: " << dynamicModel->getfxList()[i] << endl;
        //     cout << "Vx[i+1]: " << Vx[i+1] << endl;
        // }

        if(Op.regType == 1)
            QuuF = Quu + Op.lambda*commandMat_t::Identity();
        else
            QuuF = Quu;
        
        if(enableFullDDP)
        {
            //Qxx += dynamicModel->computeTensorContxx(nextVx);
            //Qux += dynamicModel->computeTensorContux(nextVx);
            //Quu += dynamicModel->computeTensorContuu(nextVx);

            // Qxx += dynamicModel->computeTensorContxx(Vx[i+1]);
            // Qux += dynamicModel->computeTensorContux(Vx[i+1]);
            // Quu += dynamicModel->computeTensorContuu(Vx[i+1]);
            // QuuF += dynamicModel->computeTensorContuu(Vx[i+1]);
        }

        // if(i > N -4){
            //cout << "QuuF: " << QuuF << endl;
            //cout << "Qu: " << Qu << endl;
            //cout << "Qux: " << Qux << endl;
        // }
        
        QuuInv = QuuF.inverse();

        if(!isPositiveDefinite(Quu))
        {
            
            //To be Implemented : Regularization (is Quu definite positive ?)
            TRACE("Quu is not positive definite");
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
            lb = lowerCommandBounds - uList[i];
            ub = upperCommandBounds - uList[i];
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

            //cout << "L: " << L << endl;
            // if(i > N-4){
            //     //cout << "index i: " << i << endl;
            //     cout << "k: " << k << endl;
            //     cout << "K: " << K << endl;
            //     cout << "Qx: " << Qx << endl;
            //     cout << "Qu: " << Qu << endl;
            //     cout << "Quu: " << Quu << endl;
            //     cout << "Qux: " << Qux << endl;
            //     cout << "QuuF: " << QuuF << endl;
            //     cout << "Op.lambda: " << Op.lambda << endl;
            //     cout << "L: " << L << endl;
            // }
        }   

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
        for(unsigned int j=0; j<commandSize; j++) {
            g_norm_i = fabs(kList[i](j,0)) / (fabs(uList[i](j,0))+1.0);
            if(g_norm_i > g_norm_max) g_norm_max = g_norm_i;
        }
        g_norm_sum += g_norm_max;
    }
    Op.g_norm = g_norm_sum/((double)(Op.n_hor));
}

void ILQRSolver::doForwardPass()
{
    updatedxList[0] = Op.xInit;
    int nargout = 2;

    stateVec_t x_unused;
    x_unused.setZero();
    commandVec_t u_NAN;
    u_NAN << sqrt(-1.0);

    //[TODO: to be optimized]
    if(!initFwdPassDone){
        //TRACE("initial forward pass");
        for(unsigned int i=0;i<N;i++)
        {
            updateduList[i] = uList[i];
            dynamicModel->cart_pole_dyn_cst_min_output(nargout, dt, updatedxList[i], updateduList[i], updatedxList[i+1], costFunction);
            costList[i] = costFunction->getc();
        }
        dynamicModel->cart_pole_dyn_cst_min_output(nargout, dt, updatedxList[N], u_NAN, x_unused, costFunction);
        costList[N] = costFunction->getc();
    }else{
        // cout << "kList[5]: " << kList[5] << endl;
        // cout << "kList[10]: " << kList[10] << endl;
        // cout << "kList[20]: " << kList[20] << endl;
        // cout << "kList[40]: " << kList[40] << endl;
        // cout << "KList[5]: " << KList[5] << endl;
        // cout << "KList[10]: " << KList[10] << endl;
        // cout << "KList[20]: " << KList[20] << endl;
        // cout << "KList[40]: " << KList[40] << endl;

        for(unsigned int i=0;i<N;i++)
        {
            updateduList[i] = uList[i] + alpha*kList[i] + KList[i]*(updatedxList[i]-xList[i]);
            // if(i<3){
            //     cout << "index i: " << i << endl;
            //     cout << "uList[i]: " << uList[i] << endl;
            //     cout << "kList[i]: " << kList[i] << endl;
            //     cout << "KList[i]: " << KList[i] << endl;
            //     cout << "updatedxList[i]: " << updatedxList[i] << endl;
            //     cout << "xList[i]: " << xList[i] << endl;
            //     cout << "updateduList[i]: " << updateduList[i] << endl;
            //     cout << "updatedxList[i]: " << updatedxList[i] << endl;
            // }
            dynamicModel->cart_pole_dyn_cst_min_output(nargout, dt, updatedxList[i], updateduList[i], updatedxList[i+1], costFunction);
            costListNew[i] = costFunction->getc();
        }
        dynamicModel->cart_pole_dyn_cst_min_output(nargout, dt, updatedxList[N], u_NAN, x_unused, costFunction);
        costListNew[N] = costFunction->getc();
    }
}

ILQRSolver::traj ILQRSolver::getLastSolvedTrajectory()
{
    lastTraj.xList = updatedxList;
    for(int i=0;i<N+1;i++)lastTraj.xList[i] += xgoal;//retrieve original state with xgoal
    lastTraj.uList = updateduList;
    lastTraj.iter = iter;
    lastTraj.finalCost = accumulate(costList.begin(), costList.end(), 0.0);
    lastTraj.finalGrad = Op.g_norm;
    lastTraj.finalLambda = log10(Op.lambda);
    lastTraj.time_forward = Op.time_forward;
    lastTraj.time_backward = Op.time_backward;
    lastTraj.time_derivative = Op.time_derivative;
    return lastTraj;
}

bool ILQRSolver::isPositiveDefinite(const commandMat_t & Quu)
{
    //Eigen::JacobiSVD<commandMat_t> svd_Quu (Quu, ComputeThinU | ComputeThinV);
    Eigen::VectorXcd singular_values = Quu.eigenvalues();

    for(long i = 0; i < Quu.cols(); ++i)
    {
        if (singular_values[i].real() < 0.)
        {
            TRACE("Matrix is not SDP");
            return false;
        }
    }
    return true;
}
