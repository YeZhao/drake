/// @file
///
/// kuka_inverse_dynamics_controller is designed to compute the torque command based on.
/// desired joint position, velocity, and acceleration, and measured joint position and velocity.
/// The final torque command is composed of inverse dynamics torque command and joint position PD 
/// controller commands.

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

class InverseDynamicsController{
 public:
  //Construct a InverseDynamicsController
  explicit InverseDynamicsController(){}

  //Compute torque command given desired and measured joint states.
  Eigen::VectorXd computeTorqueCmd(Eigen::VectorXd joint_position_desired, Eigen::VectorXd joint_velocity_desired, Eigen::VectorXd joint_accel_desired, 
                                Eigen::Map<Eigen::VectorXd> q, Eigen::Map<Eigen::VectorXd> qd, const int kNumDof, const RigidBodyTree<double>& tree_){
	  // Computing inverse dynamics torque command
	  KinematicsCache<double> cache = tree_.doKinematics(q, qd);
	  const RigidBodyTree<double>::BodyToWrenchMap no_external_wrenches;

	  // note that, the gravity term in the inverse dynamics is set to zero.
	  Eigen::VectorXd torque_command = tree_.inverseDynamics(cache, no_external_wrenches, joint_accel_desired, false);
	  
	  Eigen::VectorXd z = Eigen::VectorXd::Zero(kNumDof); 
	  KinematicsCache<double> cache0 = tree_.doKinematics(q, qd);
	  Eigen::VectorXd gravity_torque = tree_.inverseDynamics(cache0, no_external_wrenches, z, false);
	  torque_command -= gravity_torque;

	  // PD position control
	  Eigen::VectorXd position_ctrl_torque_command(kNumDof);
	  Eigen::VectorXd Kp_pos_ctrl(kNumDof); // 7 joints
	  // Kp_pos_ctrl << 225, 289, 144, 49, 324, 49, 49;
	  Kp_pos_ctrl << 100, 100, 100, 100, 100, 50, 50;
	  Eigen::VectorXd Kd_pos_ctrl(kNumDof); // 7 joints
	  // Kd_pos_ctrl << 30, 34, 24, 14, 36, 14, 14;
	  Kd_pos_ctrl << 19, 19, 19, 19, 19, 14, 14;
	  // (TODOs) Add integral control (anti-windup)
	  for (int joint = 0; joint < kNumDof; joint++) {
	    position_ctrl_torque_command(joint) = Kp_pos_ctrl(joint)*(joint_position_desired(joint) - q[joint])
	                                        + Kd_pos_ctrl(joint)*(joint_velocity_desired(joint) - qd[joint]);
	  }
	  //Combination of ID torque control and IK position control
	  torque_command += position_ctrl_torque_command;

  return torque_command;
}

 private:
  Eigen::VectorXd torque_command_;
};

}  // namespace
}  // namespace kuka_iiwa_arm
}  // namespace examples
}  // namespace drake