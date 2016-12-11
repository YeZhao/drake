/// @file
///
/// kuka_control is designed to compute the torque command based on 
/// desired joint position, velocity, and acceleration, and measured joint position and velocity.
/// Currently, the final torque command is composed of inverse dynamics torque command and joint position PD
/// controller command.
/// (TODO: Generalize this ID controller to more general types of feedback controllers)
/// Messages are sent via LCM channels.

#include <lcm/lcm-cpp.hpp>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_path.h"
#include "drake/examples/kuka_iiwa_arm/iiwa_common.h"
#include "drake/multibody/rigid_body_tree.h"

#include "drake/lcmt_robot_controller_reference.hpp"
#include "drake/lcmt_iiwa_status.hpp"
#include "drake/lcmt_iiwa_command.hpp"

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
const char* const kLcmControlRefChannel = "CONTROLLER_REFERENCE";
const char* const kLcmCommandChannel = "IIWA_COMMAND";
const int kNumJoints = 7;

class RobotController {
 public:
  /// tree is aliased
  explicit RobotController(const RigidBodyTree<double>& tree, const RigidBodyTree<double>& gravity_tree,
                           const RigidBodyTree<double>& new_gravity_tree)
      : tree_(tree), gravity_tree_(gravity_tree), new_gravity_tree_(new_gravity_tree), controller_trigger_(false) {
    VerifyIiwaTree(tree);
    lcm_.subscribe(kLcmStatusChannel,
                    &RobotController::HandleStatus, this);
    lcm_.subscribe(kLcmControlRefChannel,
                    &RobotController::HandleControl, this);
  }

  void Run() {
    int64_t cur_time_us = -1;

    // Initialize the timestamp to an invalid number so we can detect the first message.
    iiwa_status_.utime = cur_time_us;
    robot_controller_reference_.utime = cur_time_us;

    lcmt_iiwa_command iiwa_command;
    iiwa_command.num_joints = kNumJoints;
    iiwa_command.joint_position.resize(kNumJoints, 0.);
    iiwa_command.num_torques = kNumJoints;
    iiwa_command.joint_torque.resize(kNumJoints, 0.);

    Eigen::VectorXd joint_position_desired(kNumJoints); 
    Eigen::VectorXd joint_velocity_desired(kNumJoints); 
    Eigen::VectorXd joint_accel_desired(kNumJoints); 

    bool half_servo_rate_flag_ = true; // make the iiwa command get published every two servo loops

    while (true) {      
      // Call lcm handle until at least one message is processed
      while (0 == lcm_.handleTimeout(10)) { }
      DRAKE_ASSERT(iiwa_status_.utime != -1);
      DRAKE_ASSERT(robot_controller_reference_.utime != -1);
      cur_time_us = iiwa_status_.utime;

      if (controller_trigger_) {  
        const int kNumDof = 7;
        iiwa_command.utime = iiwa_status_.utime;

        // Kuka-Controller (Currently implement an inverse dynamics controller)
        // Set desired joint position, velocity and acceleration
        for (int joint = 0; joint < kNumDof; joint++){
          joint_position_desired(joint) = robot_controller_reference_.joint_position_desired[joint];
          joint_velocity_desired(joint) = robot_controller_reference_.joint_velocity_desired[joint];
          joint_accel_desired(joint) = robot_controller_reference_.joint_accel_desired[joint];
        }

        double *qptr = &iiwa_status_.joint_position_measured[0];
        Eigen::Map<Eigen::VectorXd> q(qptr, kNumDof);
        double *qdptr = &iiwa_status_.joint_velocity_estimated[0];
        Eigen::Map<Eigen::VectorXd> qd(qdptr, kNumDof);

        // Computing inverse dynamics torque command
        KinematicsCache<double> cache = tree_.doKinematics(q, qd);
        const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;

        // note that, the gravity term in the inverse dynamics is set to zero.
        Eigen::VectorXd torque_command = tree_.inverseDynamics(cache, no_external_wrenches, joint_accel_desired, false);
        
        Eigen::VectorXd z = Eigen::VectorXd::Zero(kNumDof); 
        KinematicsCache<double> cache0 = gravity_tree_.doKinematics(q, qd);
        Eigen::VectorXd gravity_torque = gravity_tree_.inverseDynamics(cache0, no_external_wrenches, z, false);
        torque_command -= gravity_torque;




        Eigen::VectorXd new_z = Eigen::VectorXd::Zero(kNumDof); 
        KinematicsCache<double> new_cache0 = new_gravity_tree_.doKinematics(q, qd);
        Eigen::VectorXd new_gravity_torque = new_gravity_tree_.inverseDynamics(new_cache0, no_external_wrenches, new_z, false);

        Eigen::VectorXd new_gripper_gravity_torque = new_gravity_torque - gravity_torque;
        std::cout << "new_gripper_gravity_torque: " << new_gripper_gravity_torque << std::endl;

        //std::cout << "id: " << tree_.get_model_instance_id() << std::endl;
        //const std::set<int> end_effector_id_set={3};
        double gripper_mass = 0.425;//tree_.getMass(end_effector_id_set);
        //std::cout << "gripper_mass" << gripper_mass << std::endl;
        //Eigen::VectorXd J = tree_.centerOfMassJacobian(cache, {7}, true);

        int base_body_or_frame_ind = 2;
        int end_effector_body_or_frame_ind = 8;
        int expressed_in_body_or_frame_ind = 0;
        bool in_terms_of_qdot = true;
        std::vector<int> v_or_qdot_indices;

        Eigen::MatrixXd J = tree_.geometricJacobian(cache, base_body_or_frame_ind, end_effector_body_or_frame_ind, 
                           expressed_in_body_or_frame_ind, in_terms_of_qdot, &v_or_qdot_indices);

        //std::cout << "Jacobian" << J << std::endl;
        Eigen::VectorXd gravity_force(6);
        gravity_force << 0, 0, -gripper_mass * 9.8;
        Eigen::VectorXd center_of_mass(3);
        center_of_mass << 0, 0, 0.0584;
        Eigen::Vector3d gripper_rotational_torque;
        double ele1 = -center_of_mass(2) * gravity_force(1) + center_of_mass(1) * gravity_force(2);
        double ele2 = center_of_mass(2) * gravity_force(0) - center_of_mass(0) * gravity_force(2);
        double ele3 = -center_of_mass(1) * gravity_force(0) + center_of_mass(0) * gravity_force(1);
        gripper_rotational_torque << ele1, ele2, ele3;
        std::cout << "gripper_rotational_torque" << gripper_rotational_torque << std::endl; 
        //std::cout << center_of_mass.cross(gravity_force) << std::endl; 
        //gripper_rotational_torque << center_of_mass.cross(gravity_force);
        //std::cout << "gravity_vector: " << gravity_vector << std::endl;
        //std::cout << "J column" << J.cols() << std::endl;
        Eigen::VectorXd torque_gripper = J.bottomRows(3).transpose() * gravity_force + J.topRows(3).transpose() * gripper_rotational_torque;
        Eigen::VectorXd torque_gripper_diff = torque_gripper - new_gripper_gravity_torque;
        std::cout << "torque_gripper" << torque_gripper << std::endl;

        //std::cout << "torque_gripper" << torque_gripper << std::endl;
        std::cout << "torque_gripper_diff" << torque_gripper_diff << std::endl;


        // PD position control
        Eigen::VectorXd position_ctrl_torque_command(kNumDof);
        Eigen::VectorXd Kp_pos_ctrl(kNumDof); // 7 joints
        Kp_pos_ctrl << 225, 361, 144, 81, 324, 36, 49;// best gains for December 9th kuka demo
        //Kp_pos_ctrl << 100, 100, 100, 100, 100, 81, 50;// original gains
        Eigen::VectorXd Kd_pos_ctrl(kNumDof); // 7 joints
        Kd_pos_ctrl << 25, 33, 20, 15, 36, 5, 14;// best gains for December 9th kuka demo, tune down damping gains from dummy critically damped gains
        //Kd_pos_ctrl << 19, 19, 19, 19, 19, 18, 14;// original gains
        // (TODOs) Add integral control (anti-windup)
        for (int joint = 0; joint < kNumJoints; joint++) {
          position_ctrl_torque_command(joint) = Kp_pos_ctrl(joint)*(joint_position_desired(joint) - iiwa_status_.joint_position_measured[joint])
                                              + Kd_pos_ctrl(joint)*(joint_velocity_desired(joint) - iiwa_status_.joint_velocity_estimated[joint]);
        }
        //Combination of ID torque control and IK position control
        torque_command += position_ctrl_torque_command;

        // -------->(For Safety) Set up iiwa position command<----------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_position[joint] = joint_position_desired(joint);
        }

        // -------->Set up iiwa torque command<-------------
        for (int joint = 0; joint < kNumJoints; joint++) {
          iiwa_command.joint_torque[joint] = torque_command(joint);
          iiwa_command.joint_torque[joint] = std::max(-150.0, std::min(150.0, iiwa_command.joint_torque[joint]));
        }

        if (half_servo_rate_flag_){
          half_servo_rate_flag_ = false;
          lcm_.publish(kLcmCommandChannel, &iiwa_command);
        }else{
          half_servo_rate_flag_ = true;
        }
      }
    
    }
  }

 private:
  void HandleStatus(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_iiwa_status* status) {
    iiwa_status_ = *status;
  }

  void HandleControl(const lcm::ReceiveBuffer* rbuf, const std::string& chan,
                    const lcmt_robot_controller_reference* input) {
    robot_controller_reference_ = *input;
    controller_trigger_ = true;
  }

  lcm::LCM lcm_;
  const RigidBodyTree<double>& tree_;
  const RigidBodyTree<double>& gravity_tree_;
  const RigidBodyTree<double>& new_gravity_tree_;
  bool controller_trigger_;// control runner wait for the first message from plan runner 
  lcmt_iiwa_status iiwa_status_;
  lcmt_robot_controller_reference robot_controller_reference_;
};

int DoMain(int argc, const char* argv[]) {

  RigidBodyTree<double> gravity_tree(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      drake::multibody::joints::kFixed);

  RigidBodyTree<double> new_gravity_tree(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_fixed_gripper.urdf",
      drake::multibody::joints::kFixed);

  RigidBodyTree<double> tree(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14_fixed_gripper.urdf",
      drake::multibody::joints::kFixed);

  RobotController runner(tree, gravity_tree, new_gravity_tree);
  runner.Run();
  return 0;
}

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::DoMain(argc, argv);
}