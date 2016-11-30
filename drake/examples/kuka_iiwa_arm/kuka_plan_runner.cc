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
#include "drake/examples/kuka_iiwa_arm/kuka_inverse_dynamics_controller.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using drake::Vector1d;
using Eigen::Vector2d;
using Eigen::Vector3d;

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

const char* const kLcmStatusChannel = "IIWA_STATUS";
const char* const kLcmCommandChannel = "IIWA_COMMAND";
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
    int cur_plan_number = plan_number_;
    int64_t cur_time_us = -1;
    int64_t start_time_us = -1;

    // Initialize the timestamp to an invalid number so we can detect the first message.
    iiwa_status_.utime = cur_time_us;

    lcmt_iiwa_command iiwa_command;
    iiwa_command.num_joints = kNumJoints;
    iiwa_command.joint_position.resize(kNumJoints, 0.);
    iiwa_command.num_torques = kNumJoints;
    iiwa_command.joint_torque.resize(kNumJoints, 0.);

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

        const int kNumDof = 7;
        const double cur_traj_time_s = static_cast<double>(cur_time_us - start_time_us) / 1e6;

        const auto q_ref = qtraj_->value(cur_traj_time_s);
        const auto qd_ref = qdtraj_->value(cur_traj_time_s);
        const auto qdd_ref = qddtraj_->value(cur_traj_time_s);

        iiwa_command.utime = iiwa_status_.utime;

        // Inverse Dynamics
        // Set desired joint position and velocity
        Eigen::VectorXd joint_position_desired(kNumDof); 
        joint_position_desired << q_ref;
        Eigen::VectorXd joint_velocity_desired(kNumDof); 
        joint_velocity_desired << qd_ref;
        Eigen::VectorXd joint_accel_desired(kNumDof); 
        joint_accel_desired << qdd_ref;
        
        double *qptr = &iiwa_status_.joint_position_measured[0];
        Eigen::Map<Eigen::VectorXd> q(qptr, kNumDof);
        double *qdptr = &iiwa_status_.joint_velocity_estimated[0];
        Eigen::Map<Eigen::VectorXd> qd(qdptr, kNumDof);

        // Computing inverse dynamics torque command
        InverseDynamicsController* torqueCtrl_;
        Eigen::VectorXd torque_command = torqueCtrl_->computeTorqueCmd(joint_position_desired, joint_velocity_desired, joint_accel_desired, q, qd, kNumDof, tree_);

        // -------->(For Safety) Set up iiwa position command<----------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_position[joint] = q_ref(joint);
        }

        // -------->Set up iiwa torque command<-------------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_torque[joint] = torque_command(joint);
          iiwa_command.joint_torque[joint] = std::max(-150.0, std::min(150.0, iiwa_command.joint_torque[joint]));
        }

        lcm_.publish(kLcmCommandChannel, &iiwa_command);
      }
    }
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
  RigidBodyTree<double> tree(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      drake::multibody::joints::kFixed);

  RobotPlanRunner runner(tree);
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