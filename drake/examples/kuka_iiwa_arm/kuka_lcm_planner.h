#include <lcm/lcm-cpp.hpp>
#include <memory>
#include <iostream>
#include <string>

#include "drake/lcmt_plan.hpp"
#include "drake/lcmt_robot_state.hpp"
#include "drake/lcmt_generic_planner_request.hpp"


#include "drake/multibody/ik_options.h"
#include "drake/multibody/rigid_body_ik.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/common/drake_path.h"

namespace drake{
namespace examples{
namespace kuka_iiwa_arm{

class KukaPlanner{
  public:
    static const char* IK_REQUEST_CHANNEL;
    static const char* PLAN_REQUEST_CHANNEL;
    static const char* IK_RESPONSE_CHANNEL;
    static const char* PLAN_RESPONSE_CHANNEL;

    KukaPlanner(std::shared_ptr<lcm::LCM> lcm = std::make_shared<lcm::LCM>(),
                std::shared_ptr<RigidBodyTree<double>> kuka = std::make_shared<RigidBodyTree<double>>()){
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
      if (chan == IK_REQUEST_CHANNEL){
        this->handleIkRequest(rbuf, status);
      }else if (chan == IK_REQUEST_CHANNEL){
        this->handlePlanRequest(rbuf, status);
      }
    }

  protected:
    std::shared_ptr<RigidBodyTree<double>> kuka_;

    virtual void handleIkRequest(const lcm::ReceiveBuffer*,
                                 const lcmt_generic_planner_request*) = 0;

    virtual void handlePlanRequest(const lcm::ReceiveBuffer*,
                                   const lcmt_generic_planner_request*) = 0;

    void publishIkResponse(Eigen::VectorXd* pose, int info){
      lcmt_plan plan;

      // TODO: build the message from the pose vector
      plan.utime = 0; // we're not using the utime on this message
      plan.num_states = 1;
      lcmt_robot_state state;
      state.timestamp = 0;
      state.num_robots = 1;
      state.robot_name.push_back(std::string("iiwa"));
      state.num_joints = kuka_->get_num_positions();
      for (int i=0; i<kuka_->get_num_positions(); i++){
        state.joint_name.push_back(std::string("joint %i") + std::to_string(i));
        state.joint_position.push_back((*pose)[i]);
        state.joint_position.push_back(0.0);
      }
      plan.plan.push_back(state);
      plan.info.push_back(info);


//         int64_t    timestamp;
//         int32_t    num_robots;
//         std::vector< std::string > robot_name;
//         int32_t    num_joints;
//         std::vector< int32_t > joint_robot;
//         std::vector< std::string > joint_name;
//         std::vector< float > joint_position;
//         std::vector< float > joint_velocity;
      std::cout << "Finished building pose\n";
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
  protected:
    virtual void handleIkRequest(const lcm::ReceiveBuffer* rbuf,
                                 const lcmt_generic_planner_request* status) {
      std::cout << "Handling ik request from KukaIkPlanner\n";
      // TODO compute IK and publish response
      Eigen::VectorXd pose = Eigen::VectorXd::Zero(kuka_->get_num_positions());
      int info = 0;
      std::cout << "Publishing an empty pose\n";
      publishIkResponse(&pose,info);
      std::cout << "Finished publishing\n";
    }

    virtual void handlePlanRequest(const lcm::ReceiveBuffer* rbuf,
                           const lcmt_generic_planner_request* status) {
      std::cout << "Handling plan request from KukaDircolPlanner\n";  
      const int num_timesteps = 20;
      Eigen::VectorXd time_vec(num_timesteps);
      Eigen::MatrixXd traj(kuka_->get_num_positions(), num_timesteps);

      // TODO compute IK and publish response
      std::cout << "Publishing an empty trajectory\n";
      publishPlanResponse(&time_vec, &traj);
      std::cout << "Finished publishing\n";
    }
};


class KukaDircolPlanner : public KukaIkPlanner {
  protected:
    // override the plan request, but not the IK request
    virtual void handlePlanRequest(const lcm::ReceiveBuffer* rbuf,
                                   const lcmt_generic_planner_request* status) {

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