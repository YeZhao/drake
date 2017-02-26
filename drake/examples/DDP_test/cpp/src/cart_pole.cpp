#include "cart_pole.h"
#include <math.h>

#define pi M_PI

const double CartPole::mc=10;
const double CartPole::mp=1;
const double CartPole::l=0.5;
const double CartPole::g=9.81;

const double CartPole::k=1000.0;
const double CartPole::R=200.0;
const double CartPole::Jm=138*1e-7;
const double CartPole::Jl=0.1;
const double CartPole::fvm=0.01;
const double CartPole::Cf0=0.1;
const double CartPole::a=10.0;

CartPole::CartPole(double& mydt, unsigned int& myT)
{
    stateNb=4;
    commandNb=1;
    dt = mydt;
    T = myT;
    Id.setIdentity();

    A.setZero();
    A(0,1) = 1.0;
    A(2,3) = 1.0;
    A(1,0) = -((k/Jl)+(k/(Jm*R*R)));
    A(1,1) = -(fvm/Jm);
    A(1,3) = -((fvm*k)/Jm);
    A(3,0) = 1.0/Jl;
    Ad = (A*dt+Id);

    A13atan = dt*(2.0*Jm*R/(pi*Jl))*Cf0;
    A33atan = dt*(2.0/(pi*Jl))*Cf0;

    B <<  0.0,
          k/(R*Jm),
          0.0,
          0.0;
    Bd = dt*B;

    fxBase <<   1.0,      dt,      0.0,      0.0,
                dt*(-(k/Jl)-(k/(Jm*R*R))),     1 - dt*(fvm/Jm),      0.0,      -dt*((fvm*k)/Jm),
                0.0,      0.0,      1.0,      dt,
                dt/Jl,      0.0,      0.0,      1.0;

    fxx[0].setZero();
    fxx[1].setZero();
    fxx[2].setZero();
    fxx[3].setZero();

    fxu[0].setZero();
    fxu[0].setZero();
    fuBase << 0.0,
              k/(R*Jm),
              0.0,
              0.0;
    fu = dt* fuBase;
    fuu[0].setZero();
    fux[0].setZero();
    fxu[0].setZero();

    QxxCont.setZero();
    QuuCont.setZero();
    QuxCont.setZero();

    lowerCommandBounds << -50.0;
    upperCommandBounds << 50.0;

    H.setZero();
    C.setZero();
    G.setZero();
    Bu.setZero();
    A1.setZero();
    A2.setZero();
    A3.setZero();
    A4.setZero();
    B1.setZero();
    B2.setZero();
    B3.setZero();
    B4.setZero();
    IdentityMat.setIdentity();

    Xp1.setZero();
    Xp2.setZero();
    Xp3.setZero();
    Xp4.setZero();

    Xm1.setZero();
    Xm2.setZero();
    Xm3.setZero();
    Xm4.setZero();
}

stateVec_t CartPole::cart_pole_dynamics(const stateVec_t& X, const commandVec_t& U)
{
    H << mc + mp, mp*l*cos(X(1,0)),
         mp*l*cos(X(1,0)), mp*pow(l,2.0);
    C << 0, -mp*X(3,0)*l*sin(X(1,0)),
         0, 0;
    G << 0,
         mp*g*l*sin(X(1,0));
    Bu << 1,
         0;     
    stateVec_half_t velocity;
    velocity << X(2),
                X(3);
    stateVec_half_t accel = - H.inverse() * (C*velocity + G - Bu*U);

    stateVec_t X_new;
    X_new << velocity(0),
             velocity(1),
             accel(0),
             accel(1);

    return X_new;
}

void CartPole::cart_pole_dyn_cst(const int& nargout, const double& dt, const stateVecTab_t& xList, const commandVecTab_t& uList, const stateVec_t& xgoal, stateVecTab_t& FList, 
                                stateMatTab_t& fx, stateR_commandC_tab_t& fu, stateVecTab_t& cx, commandVecTab_t& cu, stateMatTab_t& cxx, stateR_commandC_tab_t& cxu, commandMatTab_t& cuu, double& c){
    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    // costFunctionCartPole
    TRACE_CART_POLE("initialize dimensions\n");
    int N = xList.size();//[TODO: to be checked]
    int n = xList[0].rows();
    int m = uList[0].rows();
    CostFunctionCartPole *costFunctionCartPole = &costFunction_cart_pole;

    c = 0;
    double R = 0.1; //later, use costFunction->getR()
    stateMat_t AA;
    stateVec_t BB;
    AA.setZero();
    BB.setZero();

    TRACE_CART_POLE("compute cost function\n");

    commandMat_t c_mat_to_scalar;

    stateMatTab_t A_temp;
    stateR_commandC_tab_t B_temp;
    A_temp.resize(T);
    B_temp.resize(T);
    
    if(nargout == 2){
        const int nargout_update1 = 1;        
        for(unsigned int k=0;k<N;k++){
            if(isNan(uList[k])){ //[double check the type of c, scalar or 1x1 matirx]
                TRACE_CART_POLE("before the update1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose()) * costFunctionCartPole->getQf() * (xList[k] - xgoal);
                c += c_mat_to_scalar(0,0);
                TRACE_CART_POLE("after the update1\n");
            }else{
                TRACE_CART_POLE("before the update2\n");
                FList[k] = update(nargout_update1, dt, xList[k], uList[k], AA, BB);
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunctionCartPole->getQf()*(xList[k] - xgoal);
                TRACE_CART_POLE("after the update2\n");
                c_mat_to_scalar += 0.5*uList[k].transpose()*R*uList[k];
                c += c_mat_to_scalar(0,0);
            }
        }

    }else{
        const int nargout_update2 = 1;
        for(unsigned int k=0;k<N;k++){
            if(isNan(uList[k])){ //[double check the type of c, scalar or 1x1 matirx]
                TRACE_CART_POLE("before the update3\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunctionCartPole->getQf()*(xList[k] - xgoal);
                c += c_mat_to_scalar(0,0);
                TRACE_CART_POLE("after the update3\n");
            }else{
                TRACE_CART_POLE("before the update4\n");
                TRACE_CART_POLE("pass here\n");
                FList[k] = update(nargout_update2, dt, xList[k], uList[k], AA, BB);
                TRACE_CART_POLE("before the update4-1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunctionCartPole->getQf()*(xList[k] - xgoal);
                TRACE_CART_POLE("after the update4\n");

                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunctionCartPole->getR()*uList[k];
                c += c_mat_to_scalar(0,0); // TODO: to be checked
                TRACE_CART_POLE("after the update5\n");
                
                A_temp[k] = AA;
                B_temp[k] = BB;
            }
        }

        stateVec_t cx_temp;
        
        TRACE_CART_POLE("compute dynamics and cost derivative\n");

        for(unsigned int k=0;k<N-1;k++){
            fx[k] = A_temp[k];
            fu[k] = B_temp[k];
            cx_temp << xList[k](0,0)-xgoal(0), xList[k](1,0)-xgoal(1), xList[k](2,0)-xgoal(2), xList[k](3,0)-xgoal(3);
            cx[k] = costFunctionCartPole->getQ()*cx_temp;
            cu[k] = costFunctionCartPole->getR()*uList[k];
            cxx[k] = costFunctionCartPole->getQ();
            cxu[k].setZero();
            cuu[k] = costFunctionCartPole->getR();
        }
        TRACE_CART_POLE("update the final value of cost derivative \n");

        cx[N-1] = costFunctionCartPole->getQf()*(xList[N-1]-xgoal);
        cu[N-1] = costFunctionCartPole->getR()*uList[N-1];
        cxx[N-1] = costFunctionCartPole->getQf();
        cxu[N-1].setZero();
        cuu[N-1] = costFunctionCartPole->getR();

        TRACE_CART_POLE("set unused matrices to zero \n");

        // the following matrices and scalars are set to Zero instead of empty, not supported by Eigen.
        //fxx, fxu, fuu are not defined since never used
        for(unsigned int k=0;k<N;k++){
            FList[k].setZero();
        }    
        c = 0;
    }
    TRACE_CART_POLE("finish cart_pole_dyn_cst\n");
}

stateVec_t CartPole::update(const int& nargout, const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B){
    // 4th-order Runge-Kutta step
    TRACE_CART_POLE("update: 4th-order Runge-Kutta step\n");
    Xdot1 = cart_pole_dynamics(X, U);
    Xdot2 = cart_pole_dynamics(X + 0.5*dt*Xdot1, U);
    Xdot3 = cart_pole_dynamics(X + 0.5*dt*Xdot2, U);
    Xdot4 = cart_pole_dynamics(X + 0.5*dt*Xdot3, U);
    stateVec_t X_new;
    X_new = X + (dt/6)*(Xdot1 + Xdot2 + Xdot3 + Xdot4);
    
    TRACE_CART_POLE("update: X_new\n");

    if(nargout > 1){
        int n = X.size();
        int m = U.size();

        double delta = 1e-7;
        stateMat_t Dx;
        commandMat_t Du;
        Dx.setIdentity();
        Dx = delta*Dx;
        Du.setIdentity();
        Du = delta*Du;

        for(unsigned int i=0;i<n;i++){
            Xp1 = cart_pole_dynamics(X+Dx.col(i),U);
            Xm1 = cart_pole_dynamics(X-Dx.col(i),U);
            A1.col(i) = (Xp1 - Xm1)/(2*delta);

            Xp2 = cart_pole_dynamics(X+0.5*dt*Xdot1+Dx.col(i),U);
            Xm2 = cart_pole_dynamics(X+0.5*dt*Xdot1-Dx.col(i),U);
            A2.col(i) = (Xp2 - Xm2)/(2*delta);

            Xp3 = cart_pole_dynamics(X+0.5*dt*Xdot2+Dx.col(i),U);
            Xm3 = cart_pole_dynamics(X+0.5*dt*Xdot2-Dx.col(i),U);
            A3.col(i) = (Xp3 - Xm3)/(2*delta);

            Xp4 = cart_pole_dynamics(X+0.5*dt*Xdot3+Dx.col(i),U);
            Xm4 = cart_pole_dynamics(X+0.5*dt*Xdot3-Dx.col(i),U);
            A4.col(i) = (Xp4 - Xm4)/(2*delta);
        }

        for(unsigned int i=0;i<m;i++){
            Xp1 = cart_pole_dynamics(X,U+Du.col(i));
            Xm1 = cart_pole_dynamics(X,U-Du.col(i));
            B1.col(i) = (Xp1 - Xm1)/(2*delta);

            Xp2 = cart_pole_dynamics(X+0.5*dt*Xdot1,U+Du.col(i));
            Xm2 = cart_pole_dynamics(X+0.5*dt*Xdot1,U-Du.col(i));
            B2.col(i) = (Xp2 - Xm2)/(2*delta);

            Xp3 = cart_pole_dynamics(X+0.5*dt*Xdot2,U+Du.col(i));
            Xm3 = cart_pole_dynamics(X+0.5*dt*Xdot2,U-Du.col(i));
            B3.col(i) = (Xp3 - Xm3)/(2*delta);

            Xp4 = cart_pole_dynamics(X+0.5*dt*Xdot3,U+Du.col(i));
            Xm4 = cart_pole_dynamics(X+0.5*dt*Xdot3,U-Du.col(i));
            B4.col(i) = (Xp4 - Xm4)/(2*delta);
        }

        A = (IdentityMat + A4 * dt/6)*(IdentityMat + A3 * dt/3)*(IdentityMat + A2 * dt/3)*(IdentityMat + A1 * dt/6);
        B = B4 * dt/6 + (IdentityMat + A4 * dt/6) * B3 * dt/3 + (IdentityMat + A4 * dt/6)*(IdentityMat + A3 * dt/3)* B2 * dt/3 + (IdentityMat + (dt/6)*A4)*(IdentityMat + (dt/3)*A3)*(IdentityMat + (dt/3)*A2)*(dt/6)*B1;
    }
    TRACE_CART_POLE("update: X_new\n");
    return X_new;
}

void CartPole::grad(const double& dt, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateVec_t& B){
    int n = X.size();
    int m = U.size();

    double delta = 1e-7;
    stateMat_t Dx;
    commandMat_t Du;
    Dx.setIdentity();
    Dx = delta*Dx;
    Du.setIdentity();
    Du = delta*Du;

    stateMat_t AA;
    stateVec_t BB;
    AA.setZero();
    BB.setZero();

    int nargout = 1;
    for(unsigned int i=0;i<n;i++){
        Xp = update(nargout, dt, X+Dx.col(i), U, AA, BB);
        Xm = update(nargout, dt, X-Dx.col(i), U, AA, BB);
        A.col(i) = (Xp - Xm)/(2*delta);
    }

    for(unsigned int i=0;i<m;i++){
        Xp = update(nargout, dt, X, U+Du.col(i), AA, BB);
        Xm = update(nargout, dt, X, U-Du.col(i), AA, BB);
        B.col(i) = (Xp - Xm)/(2*delta);
    }
}

void CartPole::hessian(const double& dt, const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu){
    int n = X.size();
    int m = U.size();

    double delta = 1e-5;
    stateMat_t Dx;
    commandMat_t Du;
    Dx.setIdentity();
    Dx = delta*Dx;
    Du.setIdentity();
    Du = delta*Du;

    stateMat_t Ap;
    Ap.setZero();
    stateMat_t Am;
    Am.setZero();
    stateVec_t B;
    B.setZero();

    for(unsigned int i=0;i<n;i++){
        fxx[i].setZero();
        fxu[i].setZero();
        fuu[i].setZero();
    }

    for(unsigned int i=0;i<n;i++){
        grad(dt, X+Dx.col(i), U, Ap, B);
        grad(dt, X-Dx.col(i), U, Am, B);
        fxx[i] = (Ap - Am)/(2*delta);
    }

    stateVec_t Bp;
    Bp.setZero();
    stateVec_t Bm;
    Bm.setZero();

    for(unsigned int j=0;j<m;j++){
        grad(dt, X, U+Du.col(j), Ap, Bp);
        grad(dt, X, U-Du.col(j), Am, Bm);
        fxu[j] = (Ap - Am)/(2*delta);
        fuu[j] = (Bp - Bm)/(2*delta);
    }

}

stateVec_t CartPole::computeNextState(double& dt, const stateVec_t& X,const stateVec_t& Xdes,const commandVec_t& U)
{
    stateVec_t result = Ad*X + Bd*U;
    result(1,0)+=A13atan*atan(a*(X(3,0)+Xdes(3,0)));
    result(3,0)+=A33atan*atan(a*(X(3,0)+Xdes(3,0)));

    return result;
}

void CartPole::computeAllModelDeriv(double& dt, const stateVec_t& X,const stateVec_t& Xdes,const commandVec_t& U)
{
    Xreal = X + Xdes;
    fx = fxBase;
    fx(1,3) += A13atan*(a/(1+a*a*Xreal(3,0)*Xreal(3,0)));
    fx(3,3) -= A33atan*(a/(1+(a*a*Xreal(3,0)*Xreal(3,0))));
    fxx[3](1,3) = -((2*dt*Jm*R)/(pi*Jl))*Cf0*((2*a*a*a*Xreal(3,0))/((1+(a*a*Xreal(3,0)*Xreal(3,0)))*(1+(a*a*Xreal(3,0)*Xreal(3,0)))));
    fxx[3](3,3) = +((2*dt*Cf0)/(pi*Jl))*((2*a*a*a*Xreal(3,0))/((1+(a*a*Xreal(3,0)*Xreal(3,0)))*(1+(a*a*Xreal(3,0)*Xreal(3,0)))));
}

stateMat_t CartPole::computeTensorContxx(const stateVec_t& nextVx)
{
    QxxCont = nextVx[3]*fxx[3];
    return QxxCont;
}

commandMat_t CartPole::computeTensorContuu(const stateVec_t& nextVx)
{
    return QuuCont;
}

commandR_stateC_t CartPole::computeTensorContux(const stateVec_t& nextVx)
{
    return QuxCont;
}
