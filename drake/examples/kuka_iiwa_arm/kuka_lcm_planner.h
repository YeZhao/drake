#include <lcm/lcm-cpp.hpp>
#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>


#include "robotlocomotion/robot_plan_t.hpp"
#include "drake/lcmt_generic_planner_request.hpp"
#include "drake/lcmt_matlab_plan_request.hpp"
#include "drake/lcmt_matlab_plan_response.hpp"


#include "drake/multibody/ik_options.h"
#include "drake/multibody/rigid_body_ik.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/common/drake_path.h"

#include "rigid_body_plan_parser.h"
// #include "kuka_dircol_params.h"

namespace drake{
namespace examples{
namespace kuka_iiwa_arm{
namespace{
bot_core::robot_state_t lcmRobotState(double, Eigen::VectorXd, RigidBodyTree<double>*);

class KukaPlanner{
  public:
    static const char* IK_REQUEST_CHANNEL;
    static const char* PLAN_REQUEST_CHANNEL;
    static const char* IK_RESPONSE_CHANNEL;
    static const char* PLAN_RESPONSE_CHANNEL;

    KukaPlanner(RigidBodyTree<double>* kuka,
                std::shared_ptr<lcm::LCM> lcm){
      lcm_ = lcm;
      kuka_ = kuka;

      lcm_->subscribe(IK_REQUEST_CHANNEL,
                      &KukaPlanner::handleRequest,
                      this);
      lcm_->subscribe(PLAN_REQUEST_CHANNEL,
                      &KukaPlanner::handleRequest,
                      this);
    }

    void run(){
      std::cout << "Starting the main loop\n";
      while(true){
        lcm_->handleTimeout(1000);
      }
    }

    void handleRequest(const lcm::ReceiveBuffer* rbuf,
                       const std::string& chan,
                       const lcmt_generic_planner_request* status){
      std::cout << "Messge from " << chan << std::endl;
      std::cout << "Handling request in KukaPlanner, passing control to appropriate function\n";
      
      if (chan == IK_REQUEST_CHANNEL){
        this->handleIkRequest(status);
      }else if (chan == PLAN_REQUEST_CHANNEL){
        this->handlePlanRequest(status);
      }
    }

  protected:
    RigidBodyTree<double>* kuka_;
    std::shared_ptr<lcm::LCM> lcm_;

    virtual void handleIkRequest(const lcmt_generic_planner_request*) = 0;

    virtual void handlePlanRequest(const lcmt_generic_planner_request*) = 0;

    void publishIkResponse(Eigen::VectorXd* pose, int info){
      // build the state
      robotlocomotion::robot_plan_t plan;
      plan.utime = 0;
      plan.robot_name = "iiwa"; 
      plan.num_states = 1;
      std::cout << "tag1" << std::endl;
      Eigen::VectorXd state(kuka_->get_num_positions() + kuka_->get_num_velocities());
      state.head(kuka_->get_num_positions()) = *pose;
      state.tail(kuka_->get_num_velocities()) = Eigen::VectorXd::Zero(kuka_->get_num_velocities());
      std::cout << "tag2" << std::endl;
      
      plan.plan.push_back(lcmRobotState(0.0,state,kuka_));
      std::cout << "tag3" << std::endl;

      plan.plan_info.push_back(info);
      // publish
      std::cout << "tag4" << std::endl;
      lcm_->publish(IK_RESPONSE_CHANNEL, &plan);
      std::cout << "tag5" << std::endl;
    }

    void publishPlanResponse(Eigen::VectorXd* time, Eigen::MatrixXd* trajectory, std::vector<int> info){
      robotlocomotion::robot_plan_t plan;
      plan.utime=0;
      plan.robot_name="iiwa";
      plan.num_states = time->size();
      for (int i=0; i<plan.num_states; i++){
        std::cout << "test " << i << std::endl;
        plan.plan.push_back(lcmRobotState((*time)[i], trajectory->col(i), kuka_));
        plan.plan_info.push_back(info[i]);
      }
      plan.num_grasp_transitions = 0;
      plan.left_arm_control_type = plan.NONE;
      plan.left_leg_control_type = plan.NONE;
      plan.right_arm_control_type = plan.NONE;
      plan.right_leg_control_type = plan.NONE;

      // TODO: build the message from the time and trajectory
      lcm_->publish(PLAN_RESPONSE_CHANNEL, &plan);
    }

};

const char* KukaPlanner::IK_REQUEST_CHANNEL = "IK_REQUEST";
const char* KukaPlanner::PLAN_REQUEST_CHANNEL = "PLANNER_REQUEST";
const char* KukaPlanner::IK_RESPONSE_CHANNEL = "CANDIDATE_MANIP_IKPLAN";
const char* KukaPlanner::PLAN_RESPONSE_CHANNEL = "CANDIDATE_MANIP_PLAN";

class KukaIkPlanner : public KukaPlanner{
  public:
    KukaIkPlanner(RigidBodyTree<double>* kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaPlanner(kuka, lcm){}

  protected:
    virtual void handleIkRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling ik request from KukaIkPlanner\n";

      // unpack the LCM message
      auto constraints = parse_constraints(status->constraints, kuka_);
      std::vector<RigidBodyConstraint*> constraint_ptrs(constraints.size());
      for (unsigned int i=0; i<constraints.size(); i++){
        constraint_ptrs[i] = constraints[i].get();
      }
      auto poses = parse_json_object(status->poses);
      auto joint_names = parse_json_list(status->joint_names);
      auto seed_pose_full = parse_json_list(poses[status->seed_pose]);
      auto nominal_pose_full = parse_json_list(poses[status->nominal_pose]);
      Eigen::VectorXd seed_pose(kuka_->get_num_positions());
      Eigen::VectorXd nominal_pose(kuka_->get_num_positions());
      auto pos_idx_map = kuka_->computePositionNameToIndexMap();
      for (unsigned int i=0; i<joint_names.size(); i++){
        auto joint = joint_names[i]; 
        // ignore all of the positions associated with the floating base
        if (contains(joint,"base"))
          continue;
        int idx = pos_idx_map[joint];
        seed_pose[idx] = parse_json_double(seed_pose_full[i]);
        nominal_pose[idx] = parse_json_double(nominal_pose_full[i]);
      }
      // ignore the specified ikoptions for now
      // TODO: read the ikoptions from the message
      IKoptions ikoptions(kuka_);

      //// print some of the message components for debugging
      std::cout << "---------- Message Values ------------" << std::endl;
      std::cout << "Poses: " << status->poses << std::endl;
      std::cout << "Seed Pose: " << status->seed_pose << std::endl;
      std::cout << "Nominal Pose: " << status->nominal_pose << std::endl;
      std::cout << "End Pose: " << status->end_pose << std::endl;
      std::cout << "Joint Names: " << status->joint_names << std::endl;
      std::cout << "Options: " << status->options << std::endl;
      // print some of the computed values
      std::cout << "---------- Computed Values ------------" << std::endl;
      std::cout << "Seed Pose: \n" << seed_pose << std::endl;
      std::cout << "Nominal Pose: \n" << nominal_pose << std::endl;

      auto results = inverseKinSimple(kuka_, seed_pose, nominal_pose, constraint_ptrs, ikoptions);
      std::cout << "Computed IK" << std::endl;
      publishIkResponse(&(results.q_sol[0]),results.info[0]);
    }

    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling plan request in KukaIKPlanner\n";  
      const int num_timesteps = 20;
      Eigen::VectorXd time_vec = Eigen::VectorXd::Zero(num_timesteps);
      Eigen::MatrixXd traj = Eigen::MatrixXd::Zero(kuka_->get_num_positions(), num_timesteps);
      std::vector<int> info(num_timesteps);
      for (int i=0; i<num_timesteps; i++){
        info[i]=1;
      }
      // TODO compute IK and publish response
      std::cout << "Publishing an empty trajectory:\n" << traj << std::endl;
      publishPlanResponse(&time_vec, &traj, info);
      std::cout << "Finished publishing\n";
    }
};


class KukaDircolPlanner : public KukaIkPlanner {
  public:
    KukaDircolPlanner(RigidBodyTree<double>* kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaIkPlanner(kuka, lcm){}

  protected:
    // override the plan request, but not the IK request
    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling plan request from KukaDircolPlanner\n";
      // unpack the LCM message
      // no constraints in the plan request
      // auto constraints = parse_constraints(status->constraints, kuka_);
      // std::cout << "Num constraints " << constraints.size() << std::endl;

      auto poses = parse_json_object(status->poses);
      auto joint_names = parse_json_list(status->joint_names);
      auto start_pose_full = parse_json_list(poses[status->seed_pose]);
      auto nominal_pose_full = parse_json_list(poses[status->nominal_pose]);
      auto goal_pose_full = parse_json_list(poses[status->end_pose]);
      Eigen::VectorXd start_pose(kuka_->get_num_positions());
      Eigen::VectorXd nominal_pose(kuka_->get_num_positions());
      Eigen::VectorXd goal_pose(kuka_->get_num_positions());
      auto pos_idx_map = kuka_->computePositionNameToIndexMap();
      for (unsigned int i=0; i<joint_names.size(); i++){
        auto joint = joint_names[i]; 
        // ignore all of the positions associated with the floating base
        if (contains(joint,"base"))
          continue;
        int idx = pos_idx_map[joint];
        start_pose[idx] = parse_json_double(start_pose_full[i]);
        nominal_pose[idx] = parse_json_double(nominal_pose_full[i]);
        goal_pose[idx] = parse_json_double(goal_pose_full[i]);
      }

      
      // TODO: compute the dynamic trajectory
      
      const int num_timesteps = 20;
      // const int time_lower_bound = 2;
      // const int time_upper_bound = 6;
      
      // solvers::DircolTrajectoryOptimization dircol_traj(
      //     7 /* num inputs */, 14 /* num_states */,
      //     num_timesteps, time_lower_bound,time_upper_bound);
      // AddTrajectoryParams(kNumTimeSamples, x0, xG, &dircol_traj, std::make_unique(&kuka_));

      // print some of the message components for debugging
      std::cout << "---------- Message Values ------------" << std::endl;
      std::cout << "Constraints" << status->constraints << std::endl;
      std::cout << "Poses: " << status->poses << std::endl;
      std::cout << "Seed Pose: " << status->seed_pose << std::endl;
      std::cout << "Nominal Pose: " << status->nominal_pose << std::endl;
      std::cout << "End Pose: " << status->end_pose << std::endl;
      std::cout << "Joint Names: " << status->joint_names << std::endl;
      std::cout << "Options: " << status->options << std::endl;
      // print some of the computed values
      std::cout << "---------- Computed Values ------------" << std::endl;
      std::cout << "Start Pose: \n" << start_pose << std::endl;
      std::cout << "Nominal Pose: \n" << nominal_pose << std::endl;
      std::cout << "Goal Pose: \n" << goal_pose << std::endl;

      Eigen::VectorXd time_vec = Eigen::VectorXd::Zero(num_timesteps);
      Eigen::MatrixXd traj = Eigen::MatrixXd::Zero(kuka_->get_num_positions(), num_timesteps);
      std::vector<int> info(num_timesteps);
      for (int i=0; i<num_timesteps; i++){
        info[i]=1;
      }

      // TODO compute IK and publish response
      std::cout << "Publishing an empty trajectory:\n" << traj << std::endl;
      publishPlanResponse(&time_vec, &traj, info);
      std::cout << "Finished publishing\n";
    }
};

bot_core::robot_state_t lcmRobotState(double t, Eigen::VectorXd q, RigidBodyTree<double>* tree){
  auto pos = q.head(tree->get_num_positions());
  auto vel = q.tail(tree->get_num_velocities());

  bot_core::robot_state_t msg;
  
  msg.utime = (int) t*1e6;
  // TODO: figure out what the pose means in this context
  // msg.pose = ???
  msg.num_joints = tree->get_num_positions();
  for (int i=0; i<msg.num_joints; i++){
    msg.joint_name.push_back(tree->get_position_name(i));
    msg.joint_position.push_back(pos[i]);
    msg.joint_velocity.push_back(vel[i]);
    msg.joint_effort.push_back(0.0);
  }
  
  // msg.twist = bot_core::twist_t();
  // msg.twist.linear_velocity = bot_core::vector_3d_t()
  // msg.twist.angular_velocity = bot_core::vector_3d_t()
  // msg.force_torque = bot_core::vector_3d_t();
  for (int i=0; i<3; i++){
    msg.force_torque.l_hand_force[i] = 0.0;  
    msg.force_torque.l_hand_torque[i] = 0.0;
    msg.force_torque.r_hand_force[i] = 0.0;  
    msg.force_torque.r_hand_torque[i] = 0.0; 
  }
  return msg;

}

class KukaMatlabDircolPlanner : public KukaIkPlanner {
  public:
    static const char* MATLAB_PLAN_REQUEST_CHANNEL;
    static const char* MATLAB_PLAN_RESPONSE_CHANNEL;

    KukaMatlabDircolPlanner(RigidBodyTree<double>* kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaIkPlanner(kuka, lcm){
      
      lcm_->subscribe(MATLAB_PLAN_RESPONSE_CHANNEL,
        &KukaMatlabDircolPlanner::handleMatlabResponse, this);
    }

    void handleMatlabResponse(const lcm::ReceiveBuffer* rbuf,
                              const std::string& chan,
                              const lcmt_matlab_plan_response* plan){
      std::cout << "handling matlab response" << std::endl;
      Eigen::VectorXd t(plan->num_timesteps);
      Eigen::MatrixXd traj(kuka_->get_num_positions(), plan->num_timesteps);
      std::cout << traj << std::endl;
      std::vector<int> info;

      for (int i=0; i<plan->num_timesteps; i++){
        std::cout << "parsing time " << i << std::endl;
        t[i] = plan->time[i];
        std::cout << "parsing position " << i << std::endl;
        info.push_back(1);//plan->status;
        for (int j=0; j < kuka_->get_num_positions(); j++){
          traj(j,i) = plan->state[j][i];
        }
      }

      publishPlanResponse(&t, &traj, info);
    }
  protected:
    // override the plan request, but not the IK request
    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      auto poses = parse_json_object(status->poses);
      auto joint_names = parse_json_list(status->joint_names);
      auto start_pose_full = parse_json_list(poses[status->seed_pose]);
      auto nominal_pose_full = parse_json_list(poses[status->nominal_pose]);
      auto goal_pose_full = parse_json_list(poses[status->end_pose]);
      Eigen::VectorXd start_pose(kuka_->get_num_positions());
      Eigen::VectorXd nominal_pose(kuka_->get_num_positions());
      Eigen::VectorXd goal_pose(kuka_->get_num_positions());
      auto pos_idx_map = kuka_->computePositionNameToIndexMap();
      for (unsigned int i=0; i<joint_names.size(); i++){
        auto joint = joint_names[i]; 
        // ignore all of the positions associated with the floating base
        if (contains(joint,"base"))
          continue;
        int idx = pos_idx_map[joint];
        start_pose[idx] = parse_json_double(start_pose_full[i]);
        nominal_pose[idx] = parse_json_double(nominal_pose_full[i]);
        goal_pose[idx] = parse_json_double(goal_pose_full[i]);
      }

      lcmt_matlab_plan_request msg;
      msg.timestamp = time(NULL);
      msg.state_size = kuka_->get_num_positions()*2;
      // set position
      for (int i=0; i<kuka_->get_num_positions(); i++){
        msg.start_state.push_back(goal_pose[i]);
        msg.goal_state.push_back(goal_pose[i]);
      }
      // set velocity
      for (int i=0; i<kuka_->get_num_positions(); i++){
        msg.start_state.push_back(0);
        msg.goal_state.push_back(0);
      }

      // send request
      lcm_->publish(MATLAB_PLAN_REQUEST_CHANNEL, &msg);

    }
};

const char* KukaMatlabDircolPlanner::MATLAB_PLAN_REQUEST_CHANNEL = "MATLAB_KUKA_DIRCOL_REQUEST";
const char* KukaMatlabDircolPlanner::MATLAB_PLAN_RESPONSE_CHANNEL = "MATLAB_KUKA_DIRCOL_RESPONSE";

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake