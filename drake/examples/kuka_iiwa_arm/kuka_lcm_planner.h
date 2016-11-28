#include <lcm/lcm-cpp.hpp>
#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>


#include "drake/lcmt_plan.hpp"
#include "drake/lcmt_robot_state.hpp"
#include "drake/lcmt_generic_planner_request.hpp"


#include "drake/multibody/ik_options.h"
#include "drake/multibody/rigid_body_ik.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/common/drake_path.h"

#include "rigid_body_plan_parser.h"

namespace drake{
namespace examples{
namespace kuka_iiwa_arm{

class KukaPlanner{
  public:
    static const char* IK_REQUEST_CHANNEL;
    static const char* PLAN_REQUEST_CHANNEL;
    static const char* IK_RESPONSE_CHANNEL;
    static const char* PLAN_RESPONSE_CHANNEL;

    KukaPlanner(std::shared_ptr<RigidBodyTree<double>> kuka,
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
      std::cout << "Handling request in KukaPlanner, passing control to appropriate function\n";

      auto constraints = parse_constraints(status->constraints, kuka_.get());
      
      if (chan == IK_REQUEST_CHANNEL){
        this->handleIkRequest(status);
      }else if (chan == IK_REQUEST_CHANNEL){
        this->handlePlanRequest(status);
      }
    }

  protected:
    std::shared_ptr<RigidBodyTree<double>> kuka_;

    virtual void handleIkRequest(const lcmt_generic_planner_request*) = 0;

    virtual void handlePlanRequest(const lcmt_generic_planner_request*) = 0;

    void publishIkResponse(Eigen::VectorXd* pose, int info){
      lcmt_plan plan;

      plan.utime = 0; // we're not using the utime on this message
      plan.num_states = 1;

      // build the state
      lcmt_robot_state state;
      state.timestamp = 0;
      state.num_robots = 1;
      state.robot_name.push_back(std::string("iiwa")); 
      state.num_joints = kuka_->get_num_positions();

      for (int i=0; i<kuka_->get_num_positions(); i++){
        state.joint_name.push_back(std::string("iiwa_joint_") + std::to_string(i));
        state.joint_position.push_back((*pose)[i]);
        state.joint_velocity.push_back(0.0);
        state.joint_robot.push_back(0);
      }
      plan.plan.push_back(state);
      plan.info.push_back(info);

      lcm_->publish(IK_RESPONSE_CHANNEL, &plan);
    }

    void publishPlanResponse(Eigen::VectorXd* time, Eigen::MatrixXd* trajectory){
      lcmt_plan plan;
      // TODO: build the message from the time and trajectory
      lcm_->publish(PLAN_RESPONSE_CHANNEL, &plan);
    }

  private:
    std::shared_ptr<lcm::LCM> lcm_;

};

const char* KukaPlanner::IK_REQUEST_CHANNEL = "IK_REQUEST";
const char* KukaPlanner::PLAN_REQUEST_CHANNEL = "PLANNER_REQUEST";
const char* KukaPlanner::IK_RESPONSE_CHANNEL = "CANDIDATE_MANIP_IKPLAN";
const char* KukaPlanner::PLAN_RESPONSE_CHANNEL = "CANDIDATE_MANIP_PLAN";

class KukaIkPlanner : public KukaPlanner{
  public:
    KukaIkPlanner(std::shared_ptr<RigidBodyTree<double>> kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaPlanner(kuka, lcm){}

  protected:
    virtual void handleIkRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling ik request from KukaIkPlanner\n";

      // status->poses
      // TODO compute IK and publish response
      Eigen::VectorXd pose = Eigen::VectorXd::Zero(kuka_->get_num_positions());
      std::cout << "Returning pose:\n" << pose << std::endl;
      int info = 0;
      publishIkResponse(&pose,info);
      std::cout << "Finished publishing\n";
    }

    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {
      std::cout << "Handling plan request from KukaDircolPlanner\n";  
      const int num_timesteps = 20;
      Eigen::VectorXd time_vec(num_timesteps);
      Eigen::MatrixXd traj(kuka_->get_num_positions(), num_timesteps);

      // TODO compute IK and publish response
      std::cout << "Publishing an empty trajectory:\n" << traj << std::endl;
      publishPlanResponse(&time_vec, &traj);
      std::cout << "Finished publishing\n";
    }
};


class KukaDircolPlanner : public KukaIkPlanner {
  public:
    KukaDircolPlanner(std::shared_ptr<RigidBodyTree<double>> kuka,
                 std::shared_ptr<lcm::LCM> lcm) : KukaIkPlanner(kuka, lcm){}

  protected:
    // override the plan request, but not the IK request
    virtual void handlePlanRequest(const lcmt_generic_planner_request* status) {

      std::cout << "Handling plan request from KukaDircolPlanner\n";
      const int num_timesteps = 20;
      Eigen::VectorXd time_vec(num_timesteps);
      Eigen::MatrixXd traj(kuka_->get_num_positions(), num_timesteps);
      // TODO compute IK and publish response
      std::cout << "Publishing empty response\n";
      publishPlanResponse(&time_vec, &traj);
    }
};

} // kuka_iiwa_arm
} // examples
} // drake