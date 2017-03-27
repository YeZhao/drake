#include "kuka_arm.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

KukaArm::KukaArm(){}

KukaArm::KukaArm(double& iiwa_dt, unsigned int& iiwa_N, stateVec_t& iiwa_xgoal)
{
    stateNb = 4;
    commandNb = 1;
    dt = iiwa_dt;
    N = iiwa_N;
    xgoal = iiwa_xgoal;
    fxList.resize(N);
    fuList.resize(N);

    fxxList.resize(stateSize);
    for(unsigned int i=0;i<stateSize;i++)
        fxxList[i].resize(N);
    fxuList.resize(commandSize);
    fuuList.resize(commandSize);
    for(unsigned int i=0;i<commandSize;i++){
        fxuList[i].resize(N);
        fuuList[i].resize(N);
    }

    fxx[0].setZero();
    fxx[1].setZero();
    fxx[2].setZero();
    fxx[3].setZero();
    fuu[0].setZero();
    fux[0].setZero();
    fxu[0].setZero();

    lowerCommandBounds << -50.0;
    upperCommandBounds << 50.0;

    H.setZero();
    C.setZero();
    G.setZero();
    Bu.setZero();
    velocity.setZero();
    accel.setZero();
    Xdot_new.setZero();
    Xdot_new_thread.resize(NUMBER_OF_THREAD);
    vd_thread.resize(NUMBER_OF_THREAD);

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

    AA.setZero();
    BB.setZero();
    A_temp.resize(N);
    B_temp.resize(N);
    
    debugging_print = 0;
    finalTimeProfile.counter0_ = 0;
    finalTimeProfile.counter1_ = 0;
    finalTimeProfile.counter2_ = 0;

    initial_phase_flag_ = 1;
    q.resize(stateSize/2);
    qd.resize(stateSize/2);
    q_thread.resize(NUMBER_OF_THREAD);
    qd_thread.resize(NUMBER_OF_THREAD);
    for(unsigned int i=0;i<NUMBER_OF_THREAD;i++){
        q_thread[i].resize(stateSize/2);
        qd_thread[i].resize(stateSize/2);
    }

    finalTimeProfile.time_period1 = 0;
    finalTimeProfile.time_period2 = 0;
    finalTimeProfile.time_period3 = 0;
    finalTimeProfile.time_period4 = 0;

    if(initial_phase_flag_ == 1){
        robot_thread_ = std::make_unique<RigidBodyTree<double>>();
        parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
            GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_estimated_params_fixed_gripper.urdf",
        multibody::joints::kFixed, robot_thread_.get());

        initial_phase_flag_ = 0;
    }
}

stateVec_t KukaArm::kuka_arm_dynamics(const stateVec_t& X, const commandVec_t& tau)
{
    finalTimeProfile.counter0_ += 1;

    if(finalTimeProfile.counter0_ == 10)
        gettimeofday(&tbegin_period,NULL);

    q << X.head(stateSize/2);
    qd << X.tail(stateSize/2);
    
    KinematicsCache<double> cache_ = robot_thread_->doKinematics(q, qd);
    
    //const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;
    //gettimeofday(&tbegin_period,NULL);
    MatrixX<double> M_ = robot_thread_->massMatrix(cache_); // Inertial matrix
    //gettimeofday(&tend_period,NULL);
    //finalTimeProfile.time_period2 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

    //gettimeofday(&tbegin_period,NULL);
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext;
    
    //gettimeofday(&tend_period,NULL);
    //finalTimeProfile.time_period3 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;
    
    //gettimeofday(&tbegin_period,NULL);
    VectorX<double> bias_term_ = robot_thread_->dynamicsBiasTerm(cache_, f_ext);  // Bias term: M * vd + h = tau + J^T * lambda
    //gettimeofday(&tend_period,NULL);
    //finalTimeProfile.time_period4 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

    vd = M_.inverse()*(tau - bias_term_);
    Xdot_new << qd, vd;
    
    if(finalTimeProfile.counter0_ == 10){
        gettimeofday(&tend_period,NULL);
        finalTimeProfile.time_period1 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;
    }
    
    return Xdot_new;
}

stateVec_t KukaArm::kuka_arm_dynamicsThread1(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{

    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);

    //gettimeofday(&tbegin_period,NULL);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    //gettimeofday(&tend_period,NULL);
    //finalTimeProfile.time_period2 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

    //gettimeofday(&tbegin_period,NULL);
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    //gettimeofday(&tend_period,NULL);
    //finalTimeProfile.time_period3 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

    //gettimeofday(&tbegin_period,NULL);
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda
    //gettimeofday(&tend_period,NULL);
    //finalTimeProfile.time_period4 += ((double)(1000.0*(tend_period.tv_sec-tbegin_period.tv_sec)+((tend_period.tv_usec-tbegin_period.tv_usec)/1000.0)))/1000.0;

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread2(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    finalTimeProfile.counter1_ += 1;

    if(finalTimeProfile.counter1_ == 10)
        gettimeofday(&tbegin_period2,NULL);

    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;

    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    if(finalTimeProfile.counter1_ == 10){
        gettimeofday(&tend_period2,NULL);
        finalTimeProfile.time_period2 += ((double)(1000.0*(tend_period2.tv_sec-tbegin_period2.tv_sec)+((tend_period2.tv_usec-tbegin_period2.tv_usec)/1000.0)))/1000.0;
    }

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread3(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{

    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;

    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread4(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread5(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    finalTimeProfile.counter2_ += 1;
    if(finalTimeProfile.counter2_ == 10)
        gettimeofday(&tbegin_period3,NULL);

    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    if(finalTimeProfile.counter2_ == 10){
        gettimeofday(&tend_period3,NULL);
        finalTimeProfile.time_period3 += ((double)(1000.0*(tend_period3.tv_sec-tbegin_period3.tv_sec)+((tend_period3.tv_usec-tbegin_period3.tv_usec)/1000.0)))/1000.0;
    }

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread6(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread7(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread8(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread9(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread10(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread11(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread12(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread13(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;

    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread14(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread15(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread16(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread17(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread18(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread19(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread20(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread21(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

stateVec_t KukaArm::kuka_arm_dynamicsThread22(const stateVec_t& X_thread, const commandVec_t& tau_thread_, unsigned int index)
{
    q_thread[index] << X_thread.head(stateSize/2);
    qd_thread[index] << X_thread.tail(stateSize/2);

    KinematicsCache<double> cache_thread_ = robot_thread_->doKinematics(q_thread[index], qd_thread[index]);
    MatrixX<double> M_thread_ = robot_thread_->massMatrix(cache_thread_); // Inertial matrix
    drake::eigen_aligned_std_unordered_map<RigidBody<double> const*, drake::TwistVector<double>> f_ext_thread_;
    VectorX<double> bias_term_thread_ = robot_thread_->dynamicsBiasTerm(cache_thread_, f_ext_thread_);  // Bias term: M * vd + h = tau + J^T * lambda

    vd_thread[index] = M_thread_.inverse()*(tau_thread_ - bias_term_thread_);
    Xdot_new_thread[index] << qd_thread[index], vd_thread[index];

    return Xdot_new_thread[index];
}

KukaArm::timeprofile KukaArm::getFinalTimeProfile()
{    
    return finalTimeProfile;
}

void KukaArm::kuka_arm_dyn_cst_ilqr(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, 
                                CostFunctionKukaArm*& costFunction){
    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    if(debugging_print) TRACE_KUKA_ARM("initialize dimensions\n");
    unsigned int N = xList.size();
    
    costFunction->getc() = 0;
    AA.setZero();
    BB.setZero();

    if(debugging_print) TRACE_KUKA_ARM("compute cost function\n");

    scalar_t c_mat_to_scalar;

    if(nargout == 2){
        const int nargout_update1 = 3;        
        for(unsigned int k=0;k<N;k++){
            if(k == N-1){//isNanVec(uList[k])
                if(debugging_print) TRACE_KUKA_ARM("before the update1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose()) * costFunction->getQf() * (xList[k] - xgoal);
                costFunction->getc() += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_KUKA_ARM("after the update1\n");
            }else{
                if(debugging_print) TRACE_KUKA_ARM("before the update2\n");
                FList[k] = update(nargout_update1, xList[k], uList[k], AA, BB);
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_KUKA_ARM("after the update2\n");
                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunction->getR()*uList[k];
                costFunction->getc() += c_mat_to_scalar(0,0);
            }
        }
    }else{
        const int nargout_update2 = 3;
        for(unsigned int k=0;k<N;k++){
            if(k == N-1){//isNanVec(uList[k])
                if(debugging_print) TRACE_KUKA_ARM("before the update3\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQf()*(xList[k] - xgoal);
                costFunction->getc() += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_KUKA_ARM("after the update3\n");
            }else{
                if(debugging_print) TRACE_KUKA_ARM("before the update4\n");
                FList[k] = update(nargout_update2, xList[k], uList[k], AA, BB);//assume three outputs, code needs to be optimized
                if(debugging_print) TRACE_KUKA_ARM("before the update4-1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_KUKA_ARM("after the update4\n");

                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunction->getR()*uList[k];
                costFunction->getc() += c_mat_to_scalar(0,0); // TODO: to be checked
                if(debugging_print) TRACE_KUKA_ARM("after the update5\n");
                
                A_temp[k] = AA;
                B_temp[k] = BB;
            }
        }

        stateVec_t cx_temp;
        
        if(debugging_print) TRACE_KUKA_ARM("compute dynamics and cost derivative\n");

        for(unsigned int k=0;k<N-1;k++){
            fxList[k] = A_temp[k];
            fuList[k] = B_temp[k];
            
            cx_temp << xList[k](0,0)-xgoal(0), xList[k](1,0)-xgoal(1), xList[k](2,0)-xgoal(2), xList[k](3,0)-xgoal(3);
            costFunction->getcx()[k] = costFunction->getQ()*cx_temp;
            costFunction->getcu()[k] = costFunction->getR()*uList[k];
            costFunction->getcxx()[k] = costFunction->getQ();
            costFunction->getcux()[k].setZero();
            costFunction->getcuu()[k] = costFunction->getR();            
        }
        if(debugging_print) TRACE_KUKA_ARM("update the final value of cost derivative \n");

        costFunction->getcx()[N-1] = costFunction->getQf()*(xList[N-1]-xgoal);
        costFunction->getcu()[N-1] = costFunction->getR()*uList[N-1];
        costFunction->getcxx()[N-1] = costFunction->getQf();
        costFunction->getcux()[N-1].setZero();
        costFunction->getcuu()[N-1] = costFunction->getR();

        if(debugging_print) TRACE_KUKA_ARM("set unused matrices to zero \n");

        // the following useless matrices are set to Zero.
        //fxx, fxu, fuu are not defined since never used
        for(unsigned int k=0;k<N;k++){
            FList[k].setZero();
        }
        costFunction->getc() = 0;
    }
    if(debugging_print) TRACE_KUKA_ARM("finish kuka_arm_dyn_cst\n");
}

void KukaArm::kuka_arm_dyn_cst_min_output(const int& nargout, const double& dt, const stateVec_t& xList_curr, const commandVec_t& uList_curr, const bool& isUNan, stateVec_t& xList_next, CostFunctionKukaArm*& costFunction){
    if(debugging_print) TRACE_KUKA_ARM("initialize dimensions\n");
    unsigned int N = xList_curr.cols();

    costFunction->getc() = 0;
    AA.setZero();
    BB.setZero();

    if(debugging_print) TRACE_KUKA_ARM("compute cost function\n");

    scalar_t c_mat_to_scalar;
    xList_next.setZero();

    const int nargout_update1 = 1;
    for(unsigned int k=0;k<N;k++){
        if (isUNan){ 
            if(debugging_print) TRACE_KUKA_ARM("before the update1\n");
            c_mat_to_scalar = 0.5*(xList_curr.transpose() - xgoal.transpose()) * costFunction->getQf() * (xList_curr - xgoal);
            costFunction->getc() += c_mat_to_scalar(0,0);
            if(debugging_print) TRACE_KUKA_ARM("after the update1\n");
        }else{
            if(debugging_print) TRACE_KUKA_ARM("before the update2\n");
            xList_next = update(nargout_update1, xList_curr, uList_curr, AA, BB);
            c_mat_to_scalar = 0.5*(xList_curr.transpose() - xgoal.transpose())*costFunction->getQ()*(xList_curr - xgoal);
            if(debugging_print) TRACE_KUKA_ARM("after the update2\n");
            c_mat_to_scalar += 0.5*uList_curr.transpose()*costFunction->getR()*uList_curr;
            costFunction->getc() += c_mat_to_scalar(0,0);
        }
    }
    if(debugging_print) TRACE_KUKA_ARM("finish kuka_arm_dyn_cst\n");
}

void KukaArm::kuka_arm_dyn_cst_v3(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList, stateTensTab_t& fxxList, stateTensTab_t& fxuList, stateR_commandC_Tens_t& fuuList,
                                CostFunctionKukaArm*& costFunction){
    // // for a positive-definite quadratic, no control cost (indicated by the iLQG function using nans), is equivalent to u=0
    if(debugging_print) TRACE_KUKA_ARM("initialize dimensions\n");
    unsigned int N = xList.size();//[TODO: to be checked]
    
    costFunction->getc() = 0;
    AA.setZero();
    BB.setZero();

    if(debugging_print) TRACE_KUKA_ARM("compute cost function\n");

    scalar_t c_mat_to_scalar;
    if(nargout == 2){
        const int nargout_update1 = 3;        
        for(unsigned int k=0;k<N;k++){
            if (k == N-1){//[TODO: to be double checked, original condition is if(isNanVec(uList[k])){}]
                if(debugging_print) TRACE_KUKA_ARM("before the update1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose()) * costFunction->getQf() * (xList[k] - xgoal);
                costFunction->getc() += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_KUKA_ARM("after the update1\n");
            }else{
                if(debugging_print) TRACE_KUKA_ARM("before the update2\n");
                FList[k] = update(nargout_update1, xList[k], uList[k], AA, BB);
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_KUKA_ARM("after the update2\n");
                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunction->getR()*uList[k];
                costFunction->getc() += c_mat_to_scalar(0,0);
            }
        }
    }else{
        const int nargout_update2 = 1;
        for(unsigned int k=0;k<N;k++){
            if (k == N-1){ //[TODO: to be double checked, original condition is if(isNanVec(uList[k])){}]
                if(debugging_print) TRACE_KUKA_ARM("before the update3\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQf()*(xList[k] - xgoal);
                costFunction->getc() += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_KUKA_ARM("after the update3\n");
            }else{
                if(debugging_print) TRACE_KUKA_ARM("before the update4\n");
                FList[k] = update(nargout_update2, xList[k], uList[k], AA, BB);//assume one output, code needs to be optimized
                grad(xList[k], uList[k], AA, BB);
                hessian(xList[k], uList[k],fxx,fxu,fuu);
                for(unsigned int j;j<stateSize;j++)
                    fxxList[j][k] = fxx[j];
                for(unsigned int j;j<commandSize;j++){
                    fxuList[j][k] = fxu[j];
                    fuuList[j][k] = fuu[j];
                }
                
                if(debugging_print) TRACE_KUKA_ARM("before the update4-1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_KUKA_ARM("after the update4\n");

                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunction->getR()*uList[k];
                costFunction->getc() += c_mat_to_scalar(0,0); // TODO: to be checked
                if(debugging_print) TRACE_KUKA_ARM("after the update5\n");
                
                A_temp[k] = AA;
                B_temp[k] = BB;
            }
        }

        stateVec_t cx_temp;
        
        if(debugging_print) TRACE_KUKA_ARM("compute dynamics and cost derivative\n");

        for(unsigned int k=0;k<N-1;k++){
            fxList[k] = A_temp[k];
            fuList[k] = B_temp[k];
            
            cx_temp << xList[k](0,0)-xgoal(0), xList[k](1,0)-xgoal(1), xList[k](2,0)-xgoal(2), xList[k](3,0)-xgoal(3);
            costFunction->getcx()[k] = costFunction->getQ()*cx_temp;
            costFunction->getcu()[k] = costFunction->getR()*uList[k];
            costFunction->getcxx()[k] = costFunction->getQ();
            costFunction->getcux()[k].setZero();
            costFunction->getcuu()[k] = costFunction->getR();
        }
        if(debugging_print) TRACE_KUKA_ARM("update the final value of cost derivative \n");

        costFunction->getcx()[N-1] = costFunction->getQf()*(xList[N-1]-xgoal);//[TODO: double check whether there is - xgoal]
        costFunction->getcu()[N-1] = costFunction->getR()*uList[N-1];
        costFunction->getcxx()[N-1] = costFunction->getQf();
        costFunction->getcux()[N-1].setZero();
        costFunction->getcuu()[N-1] = costFunction->getR();

        if(debugging_print) TRACE_KUKA_ARM("set unused matrices to zero \n");

        // the following useless matrices and scalars are set to Zero.
        for(unsigned int k=0;k<N;k++){
            FList[k].setZero();
        }
        costFunction->getc() = 0;
    }
    if(debugging_print) TRACE_KUKA_ARM("finish kuka_arm_dyn_cst\n");
}

void KukaArm::kuka_arm_dyn_cst_udp(const int& nargout, const stateVecTab_t& xList, const commandVecTab_t& uList, stateVecTab_t& FList,
                                CostFunctionKukaArm*& costFunction){
    if(debugging_print) TRACE_KUKA_ARM("initialize dimensions\n");
    unsigned int N = xList.size();
    
    costFunction->getc() = 0;
    AA.setZero();
    BB.setZero();
    if(debugging_print) TRACE_KUKA_ARM("compute cost function\n");

    scalar_t c_mat_to_scalar;
    c_mat_to_scalar.setZero();

    if(nargout == 2){
        const int nargout_update1 = 3;        
        for(unsigned int k=0;k<N;k++){
            if (k == N-1){
                if(debugging_print) TRACE_KUKA_ARM("before the update1\n");
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose()) * costFunction->getQf() * (xList[k] - xgoal);
                costFunction->getc() += c_mat_to_scalar(0,0);
                if(debugging_print) TRACE_KUKA_ARM("after the update1\n");
            }else{
                if(debugging_print) TRACE_KUKA_ARM("before the update2\n");
                FList[k] = update(nargout_update1, xList[k], uList[k], AA, BB);
                c_mat_to_scalar = 0.5*(xList[k].transpose() - xgoal.transpose())*costFunction->getQ()*(xList[k] - xgoal);
                if(debugging_print) TRACE_KUKA_ARM("after the update2\n");
                c_mat_to_scalar += 0.5*uList[k].transpose()*costFunction->getR()*uList[k];
                costFunction->getc() += c_mat_to_scalar(0,0);
            }
        }
    }else{
        stateVec_t cx_temp;
        if(debugging_print) TRACE_KUKA_ARM("compute cost derivative\n");
        for(unsigned int k=0;k<N-1;k++){
            cx_temp << xList[k](0,0)-xgoal(0), xList[k](1,0)-xgoal(1), xList[k](2,0)-xgoal(2), xList[k](3,0)-xgoal(3);
            costFunction->getcx()[k] = costFunction->getQ()*cx_temp;
            costFunction->getcu()[k] = costFunction->getR()*uList[k];
            costFunction->getcxx()[k] = costFunction->getQ();
            costFunction->getcux()[k].setZero();
            costFunction->getcuu()[k] = costFunction->getR();            
        }
        if(debugging_print) TRACE_KUKA_ARM("update the final value of cost derivative \n");
        costFunction->getcx()[N-1] = costFunction->getQf()*(xList[N-1]-xgoal);
        costFunction->getcu()[N-1] = costFunction->getR()*uList[N-1];
        costFunction->getcxx()[N-1] = costFunction->getQf();
        costFunction->getcux()[N-1].setZero();
        costFunction->getcuu()[N-1] = costFunction->getR();
        if(debugging_print) TRACE_KUKA_ARM("set unused matrices to zero \n");

        // the following useless matrices and scalars are set to Zero.
        for(unsigned int k=0;k<N;k++){
            FList[k].setZero();
        }
        costFunction->getc() = 0;
    }
    if(debugging_print) TRACE_KUKA_ARM("finish kuka_arm_dyn_cst\n");
}

stateVec_t KukaArm::update(const int& nargout, const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B){
    // 4th-order Runge-Kutta step
    if(debugging_print) TRACE_KUKA_ARM("update: 4th-order Runge-Kutta step\n");

    gettimeofday(&tbegin_period4,NULL);

    Xdot1 = kuka_arm_dynamics(X, U);
    Xdot2 = kuka_arm_dynamics(X + 0.5*dt*Xdot1, U);
    Xdot3 = kuka_arm_dynamics(X + 0.5*dt*Xdot2, U);
    Xdot4 = kuka_arm_dynamics(X + dt*Xdot3, U);
    stateVec_t X_new;
    X_new = X + (dt/6)*(Xdot1 + 2*Xdot2 + 2*Xdot3 + Xdot4);
    
    if(debugging_print) TRACE_KUKA_ARM("update: X_new\n");

    if(nargout > 1){
        unsigned int n = X.size();
        unsigned int m = U.size();

        double delta = 1e-7;
        stateMat_t Dx;
        commandMat_t Du;
        Dx.setIdentity();
        Dx = delta*Dx;
        Du.setIdentity();
        Du = delta*Du;

        for(unsigned int i=0;i<n;i++){
            Xp1 = kuka_arm_dynamics(X+Dx.col(i),U);
            Xm1 = kuka_arm_dynamics(X-Dx.col(i),U);
            A1.col(i) = (Xp1 - Xm1)/(2*delta);

            Xp2 = kuka_arm_dynamics(X+0.5*dt*Xdot1+Dx.col(i),U);
            Xm2 = kuka_arm_dynamics(X+0.5*dt*Xdot1-Dx.col(i),U);
            A2.col(i) = (Xp2 - Xm2)/(2*delta);

            Xp3 = kuka_arm_dynamics(X+0.5*dt*Xdot2+Dx.col(i),U);
            Xm3 = kuka_arm_dynamics(X+0.5*dt*Xdot2-Dx.col(i),U);
            A3.col(i) = (Xp3 - Xm3)/(2*delta);

            Xp4 = kuka_arm_dynamics(X+0.5*dt*Xdot3+Dx.col(i),U);
            Xm4 = kuka_arm_dynamics(X+0.5*dt*Xdot3-Dx.col(i),U);
            A4.col(i) = (Xp4 - Xm4)/(2*delta);
        }

        for(unsigned int i=0;i<m;i++){
            Xp1 = kuka_arm_dynamics(X,U+Du.col(i));
            Xm1 = kuka_arm_dynamics(X,U-Du.col(i));
            B1.col(i) = (Xp1 - Xm1)/(2*delta);

            Xp2 = kuka_arm_dynamics(X+0.5*dt*Xdot1,U+Du.col(i));
            Xm2 = kuka_arm_dynamics(X+0.5*dt*Xdot1,U-Du.col(i));
            B2.col(i) = (Xp2 - Xm2)/(2*delta);

            Xp3 = kuka_arm_dynamics(X+0.5*dt*Xdot2,U+Du.col(i));
            Xm3 = kuka_arm_dynamics(X+0.5*dt*Xdot2,U-Du.col(i));
            B3.col(i) = (Xp3 - Xm3)/(2*delta);

            Xp4 = kuka_arm_dynamics(X+0.5*dt*Xdot3,U+Du.col(i));
            Xm4 = kuka_arm_dynamics(X+0.5*dt*Xdot3,U-Du.col(i));
            B4.col(i) = (Xp4 - Xm4)/(2*delta);
        }

        A = (IdentityMat + A4 * dt/6)*(IdentityMat + A3 * dt/3)*(IdentityMat + A2 * dt/3)*(IdentityMat + A1 * dt/6);
        B = B4 * dt/6 + (IdentityMat + A4 * dt/6) * B3 * dt/3 + (IdentityMat + A4 * dt/6)*(IdentityMat + A3 * dt/3)* B2 * dt/3 + (IdentityMat + (dt/6)*A4)*(IdentityMat + (dt/3)*A3)*(IdentityMat + (dt/3)*A2)*(dt/6)*B1;
    }
    if(debugging_print) TRACE_KUKA_ARM("update: X_new\n");

    gettimeofday(&tend_period4,NULL);
    finalTimeProfile.time_period4 += ((double)(1000.0*(tend_period4.tv_sec-tbegin_period4.tv_sec)+((tend_period4.tv_usec-tbegin_period4.tv_usec)/1000.0)))/1000.0;

    return X_new;
}

void KukaArm::grad(const stateVec_t& X, const commandVec_t& U, stateMat_t& A, stateR_commandC_t& B){
    unsigned int n = X.size();
    unsigned int m = U.size();

    double delta = 1e-7;
    stateMat_t Dx;
    commandMat_t Du;
    Dx.setIdentity();
    Dx = delta*Dx;
    Du.setIdentity();
    Du = delta*Du;

    AA.setZero();
    BB.setZero();

    int nargout = 1;
    for(unsigned int i=0;i<n;i++){
        Xp = update(nargout, X+Dx.col(i), U, AA, BB);
        Xm = update(nargout, X-Dx.col(i), U, AA, BB);
        A.col(i) = (Xp - Xm)/(2*delta);
    }

    for(unsigned int i=0;i<m;i++){
        Xp = update(nargout, X, U+Du.col(i), AA, BB);
        Xm = update(nargout, X, U-Du.col(i), AA, BB);
        B.col(i) = (Xp - Xm)/(2*delta);
    }
}

void KukaArm::hessian(const stateVec_t& X, const commandVec_t& U, stateTens_t& fxx, stateR_stateC_commandD_t& fxu, stateR_commandC_commandD_t& fuu){
    unsigned int n = X.size();
    unsigned int m = U.size();

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
    stateR_commandC_t B;
    B.setZero();

    for(unsigned int i=0;i<n;i++){
        fxx[i].setZero();
        fxu[i].setZero();
        fuu[i].setZero();
    }

    for(unsigned int i=0;i<n;i++){
        grad(X+Dx.col(i), U, Ap, B);
        grad(X-Dx.col(i), U, Am, B);
        fxx[i] = (Ap - Am)/(2*delta);
    }

    stateR_commandC_t Bp;
    Bp.setZero();
    stateR_commandC_t Bm;
    Bm.setZero();

    for(unsigned int j=0;j<m;j++){
        grad(X, U+Du.col(j), Ap, Bp);
        grad(X, U-Du.col(j), Am, Bm);
        fxu[j] = (Ap - Am)/(2*delta);
        fuu[j] = (Bp - Bm)/(2*delta);
    }
}

unsigned int KukaArm::getStateNb()
{
    return stateNb;
}

unsigned int KukaArm::getCommandNb()
{
    return commandNb;
}

commandVec_t& KukaArm::getLowerCommandBounds()
{
    return lowerCommandBounds;
}

commandVec_t& KukaArm::getUpperCommandBounds()
{
    return upperCommandBounds;
}

stateMatTab_t& KukaArm::getfxList()
{
    return fxList;
}

stateR_commandC_tab_t& KukaArm::getfuList()
{
    return fuList;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake
