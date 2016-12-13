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

#include "drake/util/drakeGeometryUtil.h"

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
        //std::cout << "torque_command: " << torque_command << std::endl;


        Eigen::VectorXd z = Eigen::VectorXd::Zero(kNumDof); 
        //KinematicsCache<double> cache0 = gravity_tree_.doKinematics(q, qd);
        //Eigen::VectorXd gravity_torque = gravity_tree_.inverseDynamics(cache0, no_external_wrenches, z, false);
        //std::cout << "gravity_torque: " << gravity_torque << std::endl;

/*        int kTwistSize = 6;
        Matrix6X<Scalar> net_wrench(kTwistSize, 1);
        const auto& cache_element = cache.getElement(*tree_.bodies[12]);
        const auto& body_inertia = cache_element.inertia_in_world;
        net_wrench.noalias() = body_inertia * body_acceleration;*/
        Eigen::VectorXd torque_command_comp = inverseDynamics_gravity(cache, z, false, tree_);
        std::cout << "torque_command_comp: " << torque_command_comp << std::endl;
        torque_command -= torque_command_comp;
        


        Eigen::VectorXd new_z = Eigen::VectorXd::Zero(kNumDof); 
        KinematicsCache<double> new_cache0 = new_gravity_tree_.doKinematics(q, qd);
        Eigen::VectorXd new_gravity_torque = new_gravity_tree_.inverseDynamics(new_cache0, no_external_wrenches, new_z, false);

        Eigen::VectorXd new_gripper_gravity_torque = new_gravity_torque - gravity_torque;
        //std::cout << "new_gripper_gravity_torque: " << new_gripper_gravity_torque << std::endl;
        //std::cout << "end effector index =====================" << tree_.FindBodyIndex("iiwa_link_ee") << std::endl;
        //std::cout << "base index =====================" << tree_.FindBodyIndex("base") << std::endl;
        //std::cout << "right_outer_finger =====================" << tree_.FindBodyIndex("right_outer_finger") << std::endl;
        //std::cout << "robotiq_85_base_link =====================" << tree_.FindBodyIndex("robotiq_85_base_link") << std::endl;


        //std::cout << "id: " << tree_.get_model_instance_id() << std::endl;
        //const std::set<int> end_effector_id_set={3};
        double gripper_mass = 0.425;//tree_.getMass(end_effector_id_set);
        //std::cout << "gripper_mass" << gripper_mass << std::endl;
        //Eigen::VectorXd J = tree_.centerOfMassJacobian(cache, {7}, true);

        int base_body_or_frame_ind = 1;//2;
        int end_effector_body_or_frame_ind = 11;//10;
        int expressed_in_body_or_frame_ind = 11;//10;
        bool in_terms_of_qdot = true;
        std::vector<int> v_or_qdot_indices;

        Eigen::MatrixXd J = tree_.geometricJacobian(cache, base_body_or_frame_ind, end_effector_body_or_frame_ind, 
                           expressed_in_body_or_frame_ind, in_terms_of_qdot, &v_or_qdot_indices);

        //std::cout << "Jacobian" << J << std::endl;
        Eigen::VectorXd gravity_force(3);
        gravity_force << 0, 0, -gripper_mass * 9.8;
        Eigen::VectorXd center_of_mass(3);
        center_of_mass << 0, 0, 0.0584;

        Eigen::VectorXd com_spatial_vector(6);
        com_spatial_vector << 0, 0, 0.0584, 0, 0, 0;
        
        //Eigen::Vector3d gripper_rotational_torque;
        //gripper_rotational_torque << ele1, ele2, ele3;
        //std::cout << "gripper_rotational_torque" << gripper_rotational_torque << std::endl; 

        auto T_world_to_ee = tree_.relativeTransform(cache, base_body_or_frame_ind,
                        expressed_in_body_or_frame_ind);
        auto transformed_gripper_pos = transformSpatialMotion(T_world_to_ee, com_spatial_vector);
        //std::cout << "transformed_gripper_pos: " << transformed_gripper_pos << std::endl; 
        double ele1 = -transformed_gripper_pos(2) * gravity_force(1) + transformed_gripper_pos(1) * gravity_force(2);
        double ele2 = transformed_gripper_pos(2) * gravity_force(0) - transformed_gripper_pos(0) * gravity_force(2);
        double ele3 = -transformed_gripper_pos(1) * gravity_force(0) + transformed_gripper_pos(0) * gravity_force(1);
        Eigen::Vector3d gripper_cross_torque;
        gripper_cross_torque << ele1, ele2, ele3;



        //std::cout << center_of_mass.cross(gravity_force) << std::endl; 
        //gripper_rotational_torque << center_of_mass.cross(gravity_force);
        //std::cout << "gravity_vector: " << gravity_vector << std::endl;
        //std::cout << "J column" << J.cols() << std::endl;
        //std::cout << "J.bottomRows(3).transpose():" << J.bottomRows(3).transpose() << std::endl;
        Eigen::VectorXd torque_gripper = J.bottomRows(3).transpose() * gravity_force + J.topRows(3).transpose() * gripper_cross_torque;
        Eigen::VectorXd torque_gripper_diff = torque_gripper - new_gripper_gravity_torque;
        Eigen::VectorXd torque_gripper_diff2 = torque_command_comp - gravity_torque;
        
        //std::cout << "torque_gripper" << torque_gripper << std::endl;
        std::cout << "torque_gripper_diff2" << torque_gripper_diff2 << std::endl;


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



  Eigen::VectorXd inverseDynamics_gravity(
      // TODO(#2274) Fix NOLINTNEXTLINE(runtime/references).
      KinematicsCache<double>& cache,
      const Eigen::VectorXd& vd,
      bool include_velocity_terms, const RigidBodyTree<double>& tree) const {
    cache.checkCachedKinematicsSettings(
        include_velocity_terms, include_velocity_terms, "inverseDynamics");

/*    //updateCompositeRigidBodyInertias(cache);

    // TODO(#3114) pass this in as an argument
    const bool include_acceleration_terms = true;

    // Compute spatial accelerations and net wrenches that should be exerted to
    // achieve those accelerations.
    // TODO(tkoolen) should preallocate:
    int kTwistSize = 6;
    Matrix6X<double> body_accelerations(kTwistSize, tree.bodies.size());
    Matrix6X<double> net_wrenches(kTwistSize, tree.bodies.size());
    for (size_t i = 0; i < tree.bodies.size(); ++i) {
      const RigidBody<double>& body = *tree.bodies[i];

      if (body.has_parent_body()) {
        const RigidBody<double>& parent_body = *(body.get_parent());
        const auto& cache_element = cache.getElement(body);

        auto body_acceleration = body_accelerations.col(i);

        // Initialize body acceleration to acceleration of parent body.
        auto parent_acceleration =
            body_accelerations.col(parent_body.get_body_index());
        body_acceleration = parent_acceleration;

        // Add component due to joint acceleration.
        if (include_acceleration_terms) {
          const DrakeJoint& joint = body.getJoint();
          auto vd_joint = vd.middleRows(body.get_velocity_start_index(),
                                        joint.get_num_velocities());
          body_acceleration.noalias() +=
              cache_element.motion_subspace_in_world * vd_joint;
        }

        auto net_wrench = net_wrenches.col(i);
        const auto& body_inertia = cache_element.inertia_in_world;
        net_wrench.noalias() = body_inertia * body_acceleration;

      } else {
        drake::TwistVector<double> a_grav;
        a_grav << 0, 0, 0, 0, 0, -9.81;
        body_accelerations.col(i) = -a_grav.cast<double>();
        net_wrenches.col(i).setZero();
      }
    }

    // Do a backwards pass to compute joint wrenches from net wrenches (update in
    // place), and project the joint wrenches onto the joint's motion subspace to
    // find the joint torque.
    auto& joint_wrenches = net_wrenches;
    // the following will eliminate the need for another explicit instantiation:
    const auto& joint_wrenches_const = net_wrenches;

    int num_joints = 7;
    VectorX<double> torques(num_joints, 1);
    //Eigen::VectorXd torques(num_joints);
    for (ptrdiff_t i = tree.bodies.size() - 1; i >= 0; --i) {
      RigidBody<double>& body = *tree.bodies[i];
      //std::cout << "tree.body.name: " << body.get_name() << std::endl;
      if (body.has_parent_body()) {
        const auto& cache_element = cache.getElement(body);
        const auto& joint = body.getJoint();
        
        auto joint_wrench = joint_wrenches_const.col(i);

        // Compute joint torques associated with the joint wrench. Joint torque
        // tau_i can be computed from the wrench across the joint between i and
        // lambda(i) as
        //
        //    J^lambda(i)_i ' * W_i
        //
        // This can be derived from a power balance.
        const auto& motion_subspace = cache_element.motion_subspace_in_world;
        auto joint_torques = torques.middleRows(body.get_velocity_start_index(),
                                                joint.get_num_velocities());
        joint_torques.noalias() = motion_subspace.transpose() * joint_wrench;

        // The joint wrench W exerted upon body i through the joint between i
        // and lambda(i) acts in the opposite direction on lambda(i) (Newton's
        // third law). This wrench, -W, should be subtracted from the net
        // wrench exerted upon lambda(i) (similar to external wrenches), so W
        // should be *added* to the net wrench.
        const RigidBody<double>& parent_body = *(body.get_parent());
        auto parent_joint_wrench =
            joint_wrenches.col(parent_body.get_body_index());
        parent_joint_wrench += joint_wrench;
      }
    }*/

    const bool include_acceleration_terms = true;
    int num_joints = 7;

    // Compute spatial accelerations and net wrenches that should be exerted to
    // achieve those accelerations.
    // TODO(tkoolen) should preallocate:
    int kTwistSize = 6;

    unsigned int body_size_no_gripper = 12;
    Matrix6X<double> body_accelerations_no_gripper(kTwistSize, body_size_no_gripper);
    Matrix6X<double> net_wrenches_no_gripper(kTwistSize, body_size_no_gripper);
    for (size_t i = 0; i < body_size_no_gripper; ++i) {
      const RigidBody<double>& body = *tree.bodies[i];
      //std::cout << "tree.body.name (no gripper): " << body.get_name() << std::endl;
      if (body.has_parent_body()) {
        const RigidBody<double>& parent_body = *(body.get_parent());
        const auto& cache_element = cache.getElement(body);

        auto body_acceleration = body_accelerations_no_gripper.col(i);

        // Initialize body acceleration to acceleration of parent body.
        auto parent_acceleration =
            body_accelerations_no_gripper.col(parent_body.get_body_index());
        body_acceleration = parent_acceleration;
        // Add component due to joint acceleration.
        if (include_acceleration_terms) {
          const DrakeJoint& joint = body.getJoint();
          auto vd_joint = vd.middleRows(body.get_velocity_start_index(),
                                        joint.get_num_velocities());
          body_acceleration.noalias() +=
              cache_element.motion_subspace_in_world * vd_joint;
        }

        auto net_wrench_no_gripper = net_wrenches_no_gripper.col(i);
        const auto& body_inertia = cache_element.inertia_in_world;
        net_wrench_no_gripper.noalias() = body_inertia * body_acceleration;

      } else {
        drake::TwistVector<double> a_grav;
        a_grav << 0, 0, 0, 0, 0, -9.81;
        body_accelerations_no_gripper.col(i) = -a_grav.cast<double>();
        net_wrenches_no_gripper.col(i).setZero();
      }
    }
    // Do a backwards pass to compute joint wrenches from net wrenches (update in
    // place), and project the joint wrenches onto the joint's motion subspace to
    // find the joint torque.
    auto& joint_wrenches_no_gripper = net_wrenches_no_gripper;
    // the following will eliminate the need for another explicit instantiation:
    const auto& joint_wrenches_const_no_gripper = net_wrenches_no_gripper;
    VectorX<double> torques_no_gripper(num_joints, 1);
    //Eigen::VectorXd torques_no_gripper(num_joints);

    for (ptrdiff_t i = body_size_no_gripper - 1; i >= 0; --i) {
      RigidBody<double>& body = *tree.bodies[i];

      if (body.has_parent_body()) {
        const auto& cache_element = cache.getElement(body);
        const auto& joint = body.getJoint();
        
        auto joint_wrench_no_gripper = joint_wrenches_const_no_gripper.col(i);
        // Compute joint torques associated with the joint wrench. Joint torque
        // tau_i can be computed from the wrench across the joint between i and
        // lambda(i) as
        //
        //    J^lambda(i)_i ' * W_i
        //
        // This can be derived from a power balance.

        const auto& motion_subspace = cache_element.motion_subspace_in_world;
        auto joint_torques_no_gripper = torques_no_gripper.middleRows(body.get_velocity_start_index(),
                                                joint.get_num_velocities());

        joint_torques_no_gripper.noalias() = motion_subspace.transpose() * joint_wrench_no_gripper;

        // The joint wrench W exerted upon body i through the joint between i
        // and lambda(i) acts in the opposite direction on lambda(i) (Newton's
        // third law). This wrench, -W, should be subtracted from the net
        // wrench exerted upon lambda(i) (similar to external wrenches), so W
        // should be *added* to the net wrench.
        const RigidBody<double>& parent_body = *(body.get_parent());
        auto parent_joint_wrench_no_gripper =
            joint_wrenches_no_gripper.col(parent_body.get_body_index());
        parent_joint_wrench_no_gripper += joint_wrench_no_gripper;
      }
    }

    //VectorX<double> torques_gripper(num_joints, 1);
    //Eigen::VectorXd torques_gripper(num_joints);

    //torques_gripper << torques - torques_no_gripper;
    /*for(int i = 0; i < num_joints; i++){
      torques_gripper(i) = torques(i) - torques_no_gripper(i);
    }
*/
    //std::cout << "torques_gripper: " << torques_gripper << std::endl;

    return torques_no_gripper;
  }

/*
  Eigen::VectorXd frictionTorques(
      Eigen::MatrixBase<DerivedV> const& v, const RigidBodyTree<double>& tree) const {
    typedef typename DerivedV::Scalar Scalar;
    Eigen::VectorXd ret(7, 1);

    for (auto it = tree.bodies.begin(); it != tree.bodies.end(); ++it) {
      RigidBody<T>& body = **it;
      if (tree.body.has_parent_body()) {
        const DrakeJoint& joint = body.getJoint();
        int nv_joint = joint.get_num_velocities();
        int v_start_joint = body.get_velocity_start_index();
        auto v_body = v.middleRows(v_start_joint, nv_joint);
        ret.middleRows(v_start_joint, nv_joint) = joint.frictionTorque(v_body);
      }
    }

    return ret;
  }*/


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