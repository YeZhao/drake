/// @file
///
/// kuka_plan_runner is designed to wait for LCM messages contraining
/// a robot_plan_t message, and then execute the plan on an iiwa arm
/// (also communicating via LCM using the
/// lcmt_iiwa_command/lcmt_iiwa_status messages).
///
/// When a plan is received, it will immediately begin executing that
/// plan on the arm (replacing any plan in progress).

#include <lcm/lcm-cpp.hpp>

#include "robotlocomotion/robot_plan_t.hpp"

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/common/trajectories/piecewise_polynomial.h"
#include "drake/common/trajectories/piecewise_polynomial_trajectory.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/parsers/urdf_parser.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_robot_controller_reference.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using drake::Vector1d;
using Eigen::Vector2d;
using Eigen::Vector3d;

/* DDP trajectory generation */
#include <iostream>
#include <fstream>

#include "drake/examples/kuka_iiwa_arm/DDP_traj_gen/config.h"
#include "drake/examples/kuka_iiwa_arm/DDP_traj_gen/ilqrsolver.cpp"
#include "drake/examples/kuka_iiwa_arm/DDP_traj_gen/udpsolver.cpp"
#include "drake/examples/kuka_iiwa_arm/DDP_traj_gen/kuka_arm.cpp"
#include "drake/examples/kuka_iiwa_arm/DDP_traj_gen/cost_function_kuka_arm.cpp"

#include <time.h>
#include <sys/time.h>

using namespace std;
using namespace Eigen;

#define useILQRSolver 0
#define useUDPSolver 1
/* DDP trajectory generation */

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

const char* const kLcmStatusChannel = "IIWA_STATUS";
const char* const kLcmControlRefChannel = "CONTROLLER_REFERENCE";
const char* const kLcmPlanChannel = "COMMITTED_ROBOT_PLAN";
const int kNumJoints = 7;

typedef PiecewisePolynomial<double> PPType;
typedef PPType::PolynomialType PPPoly;
typedef PPType::PolynomialMatrix PPMatrix;

class RobotPlanRunner {
 public:
  /// tree is aliased
  explicit RobotPlanRunner(const RigidBodyTree<double>& tree)
      : tree_(tree), plan_number_(0) {
    VerifyIiwaTree(tree);
    lcm_.subscribe(kLcmStatusChannel,
                    &RobotPlanRunner::HandleStatus, this);
    lcm_.subscribe(kLcmPlanChannel,
                    &RobotPlanRunner::HandlePlan, this);
  }

  void Run() {
    cout << "-------- Waiting for trajectory to be sent out! --------" << endl;

    int cur_plan_number = plan_number_;
    int64_t cur_time_us = -1;
    int64_t start_time_us = -1;

    // Initialize the timestamp to an invalid number so we can detect the first message.
    iiwa_status_.utime = cur_time_us;

    lcmt_robot_controller_reference robot_controller_reference;
    robot_controller_reference.num_joints = kNumJoints;
    robot_controller_reference.joint_position_desired.resize(kNumJoints, 0.);
    robot_controller_reference.joint_velocity_desired.resize(kNumJoints, 0.);
    robot_controller_reference.joint_accel_desired.resize(kNumJoints, 0.);
    robot_controller_reference.u_nominal.resize(kNumJoints, 0.);

    while (true) {
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }

      DRAKE_ASSERT(iiwa_status_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      if (qtraj_) {
        if (plan_number_ != cur_plan_number) {
          std::cout << "Starting new plan." << std::endl;
          start_time_us = cur_time_us;
          cur_plan_number = plan_number_;
        }
        const double cur_traj_time_s = static_cast<double>(cur_time_us - start_time_us) / 1e6;

        const auto q_ref = qtraj_->value(cur_traj_time_s);
        const auto qd_ref = qdtraj_->value(cur_traj_time_s);
        const auto qdd_ref = qddtraj_->value(cur_traj_time_s);

        robot_controller_reference.utime = iiwa_status_.utime;

        for(int joint = 0; joint < kNumJoints; joint++){
          robot_controller_reference.joint_position_desired[joint] = q_ref(joint);
          robot_controller_reference.joint_velocity_desired[joint] = qd_ref(joint);
          robot_controller_reference.joint_accel_desired[joint] = qdd_ref(joint);
        }

        // publish robot controller reference to kuka control runner
        lcm_.publish(kLcmControlRefChannel, &robot_controller_reference);
      }
    }
  }

  void RunUDP() {
    struct timeval tbegin,tend;
    double texec = 0.0;
    stateVec_t xinit,xgoal;

    xinit << 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0;
    xgoal << 0.0,pi/3,0.0,pi/2,0.0,pi/4,0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0;

    // xinit << 0.0,0.0,0.0,0.0;
    // xgoal << 0.0,pi,0.0,0.0;

    double T = TimeHorizon;
    double dt = TimeStep;
    unsigned int N = (int)(T/dt);
    double tolFun = 1e-5;//relaxing default value: 1e-10;
    double tolGrad = 1e-5;//relaxing default value: 1e-10;
    unsigned int iterMax = 150;
    #if useILQRSolver
        ILQRSolver::traj lastTraj;
        KukaArm KukaArmModel(dt, N, xgoal);
        CostFunctionKukaArm costKukaArm;
        ILQRSolver testSolverKukaArm(KukaArmModel,costKukaArm,ENABLE_FULLDDP,ENABLE_QPBOX);
        testSolverKukaArm.firstInitSolver(xinit, xgoal, N, dt, iterMax, tolFun, tolGrad);    
    #endif
    #if useUDPSolver
        double scale = 1e-4;//0.1;
        UDPSolver::traj lastTraj;
        KukaArm KukaArmModel(dt, N, xgoal);
        CostFunctionKukaArm costKukaArm;
        UDPSolver testSolverKukaArm(KukaArmModel,costKukaArm,ENABLE_FULLDDP,ENABLE_QPBOX);
        testSolverKukaArm.firstInitSolver(xinit, xgoal, N, dt, scale, iterMax, tolFun, tolGrad);    
    #endif

    // run one or multiple times and then average
    unsigned int Num_run = 1;
    gettimeofday(&tbegin,NULL);
    for(unsigned int i=0;i<Num_run;i++) testSolverKukaArm.solveTrajectory();
    gettimeofday(&tend,NULL);

    lastTraj = testSolverKukaArm.getLastSolvedTrajectory();

    texec=((double)(1000*(tend.tv_sec-tbegin.tv_sec)+((tend.tv_usec-tbegin.tv_usec)/1000)))/1000.;
    texec /= Num_run;

    cout << endl;
    cout << "Number of iterations: " << lastTraj.iter + 1 << endl;
    cout << "Final cost: " << lastTraj.finalCost << endl;
    cout << "Final gradient: " << lastTraj.finalGrad << endl;
    cout << "Final lambda: " << lastTraj.finalLambda << endl;
    cout << "Execution time by time step (second): " << texec/N << endl;
    cout << "Execution time per iteration (second): " << texec/lastTraj.iter << endl;
    cout << "Total execution time of the solver (second): " << texec << endl;
    cout << "\tTime of derivative (second): " << lastTraj.time_derivative.sum() << " (" << 100.0*lastTraj.time_derivative.sum()/texec << "%)" << endl;
    cout << "\tTime of backward pass (second): " << lastTraj.time_backward.sum() << " (" << 100.0*lastTraj.time_backward.sum()/texec << "%)" << endl;
    //cout << "\t\tTime range 1 (second): " << lastTraj.time_range1.sum() << " (" << 100.0*lastTraj.time_range2.sum()/texec << "%)" << endl;
    //cout << "\t\tTime range 2 (second): " << lastTraj.time_range2.sum() << " (" << 100.0*lastTraj.time_range2.sum()/texec << "%)" << endl;
    cout << "\tTime of forward pass (second): " << lastTraj.time_forward.sum() << " (" << 100.0*lastTraj.time_forward.sum()/texec << "%)" << endl;
    
    // debugging trajectory and control outputs
    cout << "--------- final joint state trajectory ---------" << endl;
    for(unsigned int i=0;i<=N;i++){
      cout << "lastTraj.xList[" << i << "]:" << lastTraj.xList[i].transpose() - xgoal.transpose() << endl;
    }
    cout << "--------- final joint torque trajectory ---------" << endl;
    
    for(unsigned int i=0;i<=N;i++){
      cout << "lastTraj.uList[" << i << "]:" << lastTraj.uList[i].transpose() << endl;
    }

    ofstream file("results.csv",ios::out | ios::trunc);

    if(file)
    {
        file << "x,theta,xDot,thetaDot,u" << endl;
        for(unsigned int i=0;i<N;i++) file << lastTraj.xList[i](0,0) << "," << lastTraj.xList[i](1,0) << "," << lastTraj.xList[i](2,0) << "," << lastTraj.xList[i](3,0) << "," << lastTraj.uList[i](0,0) << endl;
        file << lastTraj.xList[N](0,0) << "," << lastTraj.xList[N](1,0) << "," << lastTraj.xList[N](2,0) << "," << lastTraj.xList[N](3,0) << "," << 0.0 << endl;
        file.close();
    }
    else
        cerr << "error in open file" << endl;
    cout << "-------- DDP Trajectory Generation Finished! --------" << endl;
  }

 private:
  void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    iiwa_status_ = *status;
  }

  void HandlePlan(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                  const robotlocomotion::robot_plan_t* plan) {
    std::cout << "New plan received." << std::endl;
    Eigen::MatrixXd traj_mat(kNumJoints, plan->num_states);
    traj_mat.fill(0);

    std::map<std::string, int> name_to_idx =
        tree_.computePositionNameToIndexMap();
    for (int i = 0; i < plan->num_states; ++i) {
      const auto& state = plan->plan[i];
      for (int j = 0; j < state.num_joints; ++j) {
        if (name_to_idx.count(state.joint_name[j]) == 0) {
          continue;
        }
        traj_mat(name_to_idx[state.joint_name[j]], i) =
            state.joint_position[j];
      }
    }

    std::cout << traj_mat << std::endl;

    std::vector<double> input_time;
    for (int k = 0; k < static_cast<int>(plan->plan.size()); ++k) {
      input_time.push_back(plan->plan[k].utime / 1e6);
    }
    qtraj_.reset(new PiecewisePolynomialTrajectory(traj_mat, input_time));
    qdtraj_.reset(new PiecewisePolynomialTrajectory(qtraj_->getPP().derivative(1)));
    qddtraj_.reset(new PiecewisePolynomialTrajectory(qdtraj_->getPP().derivative(1)));
    ++plan_number_;
  }

  lcm::LCM lcm_;
  const RigidBodyTree<double>& tree_;
  int plan_number_{};
  std::unique_ptr<PiecewisePolynomialTrajectory> qtraj_;
  std::unique_ptr<PiecewisePolynomialTrajectory> qdtraj_;
  std::unique_ptr<PiecewisePolynomialTrajectory> qddtraj_;
  lcmt_iiwa_status iiwa_status_;
};

int do_main(int argc, const char* argv[]) {

  auto tree = std::make_unique<RigidBodyTree<double>>();
  parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_estimated_params_fixed_gripper.urdf",
      multibody::joints::kFixed, tree.get());

  RobotPlanRunner runner(*tree);
  runner.RunUDP();
  runner.Run();

  return 0;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::do_main(argc, argv);
}
