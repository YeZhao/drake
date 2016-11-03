#include <iostream>
#include <time.h>

#include <lcm/lcm-cpp.hpp>

#include "drake/lcmt_iiwa_command.hpp"
#include "drake/lcmt_iiwa_status.hpp"

#include "drake/examples/kuka_iiwa_arm/iiwa_status.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/common/drake_path.h"


namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace{

const char* kLcmDarkCommandChannel = "IIWA_DARK_COMMAND";
const int kNumJoints = 7;

/** Handles sending and receiving of LCM messages
 *
 */
class MessageHandler {
  public:
    /** Constructor
     */
    MessageHandler(std::shared_ptr<lcm::LCM> lcm): lcm_(lcm){
      // subscribe the handler to IIWA_STATUS
      lcm_->subscribe(IiwaStatus<double>::channel(),
                      &MessageHandler::handleMessage,
                      this);
      printf("Listening to channel: %s\n",IiwaStatus<double>::channel().c_str());
      // record the time of initialization
      time(&last_received_time_);
    }

    /** Checks if it has been more than 3 seconds since the last message
     */
    bool hasTimedOut(){
      // check if it has been more than 0.5 seconds since the last message
      return difftime(time(NULL),last_received_time_) > 3;
    }

    /** Checks for new commands from LCM
     *  This must be called in the event loop in order for new messages to come in
     */
    bool handle(){
      // wait 10ms before moving releasing control of the process
      return lcm_->handleTimeout(10);
    }

    /** Waits for the first message to come in up to a timeout
     *  Returns true if a message came in successfully and false if it failed
     */
    bool waitForFirstMessage(){
      bool success = true;
      while(!this->received_first_message_ && success){
        this->handle();
        success = !(this->hasTimedOut());
      }
      return success;
    }

    void printPositions(){
      int numJoints = iiwa_status_.num_joints;
      for (int i=0; i<numJoints; i++){
        printf("Joint %i measured at: %f\n", i, iiwa_status_.joint_position_measured[i]); 
      }
    }

    lcmt_iiwa_status getStatus(){
      return iiwa_status_;
    }

    void publish(lcmt_iiwa_command iiwa_command){
      lcm_->publish(kLcmDarkCommandChannel, &iiwa_command);
    }

  private:
    lcmt_iiwa_status iiwa_status_;
    time_t last_received_time_;
    bool received_first_message_ = false;
    std::shared_ptr<lcm::LCM> lcm_;

    void handleMessage(const lcm::ReceiveBuffer* rbuf,
                       const std::string& chan,
                       const lcmt_iiwa_status* status){
      iiwa_status_ = *status;
      time(&last_received_time_);
      received_first_message_ = true;
      // TODO: do other things on status update
    }
};

int main(int argc, const char* argv[]){
  // allocate a pointer to an LCM object and create the message handler
  std::shared_ptr<lcm::LCM> lcm = std::make_shared<lcm::LCM>();
  MessageHandler handler(lcm);

  // load the model for the kuka arm
  RigidBodyTree tree(
    drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
    drake::systems::plants::joints::kFixed);

  // wait for the first status to come in
  bool success = handler.waitForFirstMessage();

  if (success){
    printf("Status has come in\n");
    printf("Publishing to channel %s\n", kLcmDarkCommandChannel);    
  }else{
    printf("Timed out waiting for status\n");
    return 1;
  }

  // initialize the iiwa command 
  lcmt_iiwa_command command;
  command.num_joints = kNumJoints;
  command.joint_position.resize(kNumJoints, 0.);
  command.num_torques = kNumJoints;
  command.joint_torque.resize(kNumJoints, 0.);
  while(!handler.hasTimedOut()){
    if (handler.handle()){
      handler.printPositions();
      command.timestamp = handler.getStatus().timestamp;
      //TODO: get the state from the MessageHandler and compute the inverse dynamics required to resist gravity 
      handler.publish(command);
    }
  }
  printf("Timed out waiting for status\n");

  return 0;
}

} // namespace
} // kuka_iiwa_arm
} // examples
} // drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::main(argc, argv);
}