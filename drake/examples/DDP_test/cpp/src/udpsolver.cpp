#include "udpsolver.h"

/* Debug */
#include <iostream>
using namespace std;
/* */

using namespace Eigen;
using Eigen::VectorXd;

UDPSolver::UDPSolver(DynamicModel& myDynamicModel, CostFunction& myCostFunction, bool fullDDP, bool QPBox)
{
    //TRACE_UDP("initialize dynamic model and cost function\n");
    dynamicModel = &myDynamicModel;
    costFunction = &myCostFunction;
    stateNb = myDynamicModel.getStateNb();
    commandNb = myDynamicModel.getCommandNb();
    enableQPBox = QPBox;
    enableFullDDP = fullDDP;

    if(enableQPBox) TRACE_UDP("Box QP is enabled");
    else TRACE_UDP("Box QP is disabled");

    if(enableFullDDP) TRACE_UDP("Full DDP is enabled");
    else TRACE_UDP("Full DDP is disabled");

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

void UDPSolver::firstInitSolver(stateVec_t& myxInit, stateVec_t& myxgoal, unsigned int& myN,
                       double& mydt, double& myscale, unsigned int& mymax_iter, double& mytolFun, double& mytolGrad)
{
    // TODO: double check opt params
    xInit = myxInit; // removed myxgoal. Double check whether this makes sense.
    xgoal = myxgoal;
    N = myN;
    dt = mydt;
    scale = myscale;

    // Eigen::VectorXd default_alpha;
    // default_alpha.setZero();
    // default_alpha << 1.0, 0.5012, 0.2512, 0.1259, 0.0631, 0.0316, 0.0158, 0.0079, 0.0040, 0.0020, 0.0010;

    //TRACE_UDP("initialize option parameters\n");
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
    uListFull.resize(N+1);
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

    // initialization in doBackwardPass
    augMatrix.resize(fullstatecommandSize, fullstatecommandSize);
    Sig.resize(fullstatecommandSize, 2*fullstatecommandSize);
    augState.resize(fullstatecommandSize, 1);
    G.resize(2*fullstatecommandSize, 1);
    D.resize(fullstatecommandSize, fullstatecommandSize);
    df.resize(fullstatecommandSize, 1);
    M.resize(fullstatecommandSize, fullstatecommandSize);
    HH.resize(fullstatecommandSize, fullstatecommandSize);
}

void UDPSolver::solveTrajectory()
{
    initializeTraj();
    
    Op.lambda = Op.lambdaInit;
    Op.dlambda = Op.dlambdaInit;
    
    // TODO: update multipliers
    //update_multipliers(Op, 1);

    for(iter=0;iter<Op.max_iter;iter++)
    {
        //TRACE_UDP("STEP 1: differentiate cost along new trajectory\n");
        if(newDeriv){
            int nargout = 7;//fx,fu,cx,cu,cxx,cxu,cuu
            for(unsigned int i=0;i<u_NAN.size();i++)
                u_NAN(i,0) = sqrt(-1.0);
            for(unsigned int i=0;i<uList.size();i++)
                uListFull[i] = uList[i];

            uListFull[uList.size()] = u_NAN;

            gettimeofday(&tbegin_time_deriv,NULL);
            dynamicModel->cart_pole_dyn_cst_udp(nargout, xList, uListFull, FList, costFunction);
            gettimeofday(&tend_time_deriv,NULL);
            Op.time_derivative(iter) = ((double)(1000*(tend_time_deriv.tv_sec-tbegin_time_deriv.tv_sec)+((tend_time_deriv.tv_usec-tbegin_time_deriv.tv_usec)/1000)))/1000.0;
            newDeriv = 0;
        }
        //TRACE_UDP("Finish STEP 1\n");

        //====== STEP 2: backward pass, compute optimal control law and cost-to-go
        backPassDone = 0;
        while(!backPassDone){
            gettimeofday(&tbegin_time_bwd,NULL);
            doBackwardPass();
            gettimeofday(&tend_time_bwd,NULL);
            Op.time_backward(iter) = ((double)(1000*(tend_time_bwd.tv_sec-tbegin_time_bwd.tv_sec)+((tend_time_bwd.tv_usec-tbegin_time_bwd.tv_usec)/1000)))/1000.0;

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
                TRACE_UDP(("\nSUCCESS: gradient norm < tolGrad\n"));
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
                double z;
                if(Op.expected > 0) {
                    z = Op.dcost/Op.expected;
                }else {
                    z = (double)(-signbit(Op.dcost));//[TODO:doublecheck]
                    TRACE_UDP("non-positive expected reduction: should not occur \n");//warning
                }
                if(z > Op.zMin){
                    fwdPassDone = 1;
                    break;
                }
            }
            if(!fwdPassDone) alpha = sqrt(-1.0);
            gettimeofday(&tend_time_fwd,NULL);
            Op.time_forward(iter) = ((double)(1000*(tend_time_fwd.tv_sec-tbegin_time_fwd.tv_sec)+((tend_time_fwd.tv_usec-tbegin_time_fwd.tv_usec)/1000)))/1000.0;
            cout << "Op.time_forward(iter): " << Op.time_forward(iter) << endl;
        }
        
        //====== STEP 4: accept step (or not), draw graphics, print status
        if (Op.debug_level > 1 && Op.last_head == Op.print_head){
            Op.last_head = 0;
            TRACE_UDP("iteration,\t cost, \t reduction, \t expected, \t gradient, \t log10(lambda) \n");
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
                    TRACE_UDP(("\nSUCCESS: cost change < tolFun\n"));

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
                if(Op.debug_level >= 0)
                    TRACE_UDP(("\nEXIT: lambda > lambdaMax\n"));
                break;
            }
        }
        //cout << "final alpha value after one iteration: " << alpha << endl;
    }

    Op.iterations = iter;

    if(!backPassDone) {
        if(Op.debug_level >= 1)
            TRACE_UDP(("\nEXIT: no descent direction found.\n"));
        
        return;    
    } else if(iter >= Op.max_iter) {
        if(Op.debug_level >= 0)
            TRACE_UDP(("\nEXIT: Maximum iterations reached.\n"));
        
        return;
    }
}

void UDPSolver::initializeTraj()
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
    
    //constants, timers, counters
    newDeriv = 1; //flgChange
    Op.lambda= Op.lambdaInit;
    Op.w_pen_l= Op.w_pen_init_l;
    Op.w_pen_f= Op.w_pen_init_f;
    Op.dcost = 0;
    Op.expected = 0;
    Op.print_head = 6;
    Op.last_head = Op.print_head;
    if(Op.debug_level > 0) TRACE_UDP("\n=========== begin UDP ===========\n");
}

void UDPSolver::standardizeParameters(tOptSet *o) {
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

void UDPSolver::doBackwardPass()
{    
    //Perform the Ricatti-Mayne backward pass with unscented transform
    if(Op.regType == 1)
        lambdaEye = Op.lambda*stateMat_t::Identity();
    else
        lambdaEye = Op.lambda*stateMat_t::Zero();

    diverge = 0;
    g_norm_sum = 0;
    Vx[N] = costFunction->getcx()[N];
    Vxx[N] = costFunction->getcxx()[N];
    dV.setZero();

    for(int i=N-1;i>=0;i--){
        //Generate sigma points from Vxx(i+1) and cuu(i+1)
        commandR_stateC_t ZeroLowerLeftMatrix;
        stateR_commandC_t ZeroUpperRightMatrix;
        ZeroLowerLeftMatrix.setZero();
        ZeroUpperRightMatrix.setZero();
        augMatrix << Vxx[i+1].inverse(), ZeroUpperRightMatrix, 
                ZeroLowerLeftMatrix, costFunction->getcuu()[i].inverse();

        gettimeofday(&tbegin_test,NULL);        
        Eigen::LLT<MatrixXd> lltOfaugMatrix(augMatrix);
        Eigen::MatrixXd S = lltOfaugMatrix.matrixL(); 
        //assume augMatrix is positive definite
        gettimeofday(&tend_test,NULL);

        //A temporary solution: check the non-PD case
        if(lltOfaugMatrix.info() == Eigen::NumericalIssue)
        {
            diverge = i;
            TRACE_UDP("Possibly non semi-positive definitie matrix!");
            return;
        }
        S = scale*S;
        Sig << S, -S;

        for(unsigned int j=0;j<2*fullstatecommandSize;j++){
            augState << xList[i+1], uList[i];
            Sig.col(j) += augState;
        }

        // Project Vx(i+1) onto sigma points
        for(unsigned int j=0;j<fullstatecommandSize;j++){
            G(j) = Vx[i+1].transpose()*S.col(j).head(stateSize);//TODO: double make sure its correctness
            G(j+fullstatecommandSize) = -G(j);
        }

        //Propagate sigma points through backwards dynamics
        for(unsigned int j=0;j<2*fullstatecommandSize;j++)
            Sig.col(j).head(stateSize) = rungeKuttaStepBackward(Sig.col(j), dt);

        //Calculate [Qu; Qx] from sigma points
        for(unsigned int j=0;j<fullstatecommandSize;j++){
            D.row(j) =  Sig.col(j).transpose() - Sig.col(j+fullstatecommandSize).transpose();
            df(j) = G(j) - G(fullstatecommandSize+j);
        }

        QxQu = D.inverse()*df;
        Qx = QxQu.head(stateSize) + costFunction->getcx()[i]; //add on one-step cost
        Qu = QxQu.tail(commandSize) + costFunction->getcu()[i]; //add on one-step cost
        
        mu.setZero();
        //Calculate Hessian w.r.t. [xList[i]; uList[i]] from sigma points
        for(unsigned int j=0;j<2*fullstatecommandSize;j++)
            mu += 1.0/(2.0*fullstatecommandSize)*Sig.col(j);

        M.setZero();
        for(unsigned int j=0;j<2*fullstatecommandSize;j++)
            M += (0.5/pow(scale, 2.0))*(Sig.col(j) - mu)*(Sig.col(j).transpose() - mu.transpose());

        HH = M.inverse();
        HH.block(0,0,stateSize,stateSize) += costFunction->getcxx()[i]; //add in one-step state cost for this timestep
        
        Qxx = HH.block(0,0,stateSize,stateSize);
        Quu = HH.block(stateSize,stateSize,commandSize,commandSize);
        Qux = HH.block(stateSize,0,commandSize,stateSize);

        if(Op.regType == 1)
            QuuF = Quu + Op.lambda*commandMat_t::Identity();

        QuuInv = QuuF.inverse();

        if(!isPositiveDefinite(Quu))
        {
            TRACE_UDP("Quu is not positive definite");
            if(Op.lambda==0.0) Op.lambda += 1e-4;
            else Op.lambda *= 10;
            backPassDone = 0;
            break;
        }

        if(enableQPBox)
        {
            //TRACE_UDP("Use Box QP");
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
            gettimeofday(&tbegin_test2,NULL);
            // Cholesky decomposition by using upper triangular matrix
            //TRACE_UDP("Use Cholesky decomposition");
            Eigen::LLT<MatrixXd> lltOfQuuF(QuuF);
            Eigen::MatrixXd L = lltOfQuuF.matrixU(); 
            //assume QuuF is positive definite
            gettimeofday(&tend_test2,NULL);

            //A temporary solution: check the non-PD case
            if(lltOfQuuF.info() == Eigen::NumericalIssue)
                {
                    diverge = i;
                    TRACE_UDP("Possibly non semi-positive definitie matrix!");
                    return;
                }

            k = - L.inverse()*L.transpose().inverse()*Qu;
            K = - L.inverse()*L.transpose().inverse()*Qux;

            // if(i > N-4){
            //     cout << "index i: " << i << endl;
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

void UDPSolver::doForwardPass()
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
        for(unsigned int i=0;i<N;i++){
            updateduList[i] = uList[i] + alpha*kList[i] + KList[i]*(updatedxList[i]-xList[i]);
            dynamicModel->cart_pole_dyn_cst_min_output(nargout, dt, updatedxList[i], updateduList[i], updatedxList[i+1], costFunction);
            costListNew[i] = costFunction->getc();
        }
        dynamicModel->cart_pole_dyn_cst_min_output(nargout, dt, updatedxList[N], u_NAN, x_unused, costFunction);
        costListNew[N] = costFunction->getc();
    }
}

UDPSolver::traj UDPSolver::getLastSolvedTrajectory()
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

bool UDPSolver::isPositiveDefinite(const commandMat_t & Quu)
{
    //Eigen::JacobiSVD<commandMat_t> svd_Quu (Quu, ComputeThinU | ComputeThinV);
    Eigen::VectorXcd singular_values = Quu.eigenvalues();

    for(long i = 0; i < Quu.cols(); ++i)
    {
        if (singular_values[i].real() < 0.)
        {
            TRACE_UDP("Matrix is not SDP");
            return false;
        }
    }
    return true;
}

stateVec_t UDPSolver::rungeKuttaStepBackward(stateAug_t augX, double& dt){
    // Backwards 4^th order Runge-Kutta step from X_{k+1} to X_k
    Xdot1 = dynamicModel->cart_pole_dynamics(augX.head(stateSize), augX.tail(commandSize));
    Xdot2 = dynamicModel->cart_pole_dynamics(augX.head(stateSize) - 0.5*dt*Xdot1, augX.tail(commandSize));
    Xdot3 = dynamicModel->cart_pole_dynamics(augX.head(stateSize) - 0.5*dt*Xdot2, augX.tail(commandSize));
    Xdot4 = dynamicModel->cart_pole_dynamics(augX.head(stateSize) - dt*Xdot3, augX.tail(commandSize));
    stateVec_t X_new;
    X_new = augX.head(stateSize) - (dt/6)*(Xdot1 + 2*Xdot2 + 2*Xdot3 + Xdot4);
    return X_new;
}