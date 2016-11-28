#import <iostream>
#import <memory>

#import "kuka_lcm_planner.h"
#include "drake/common/drake_path.h"
#include "drake/multibody/rigid_body_tree.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

int main(int argc, const char* argv[]) {
  std::string nested_braces ="{{{}{{}{}}{}}{{}{}}{{}}}";
  std::string nested_brackets ="[[[[][][]][][]]]";
  std::string test_string = "\"This is a JSON string\", other junk";
  std::string test_map= "{\"1\":\"hello\",\"2\":\"world\",\"test\":4,\"test\":true,\"test2\":false,\"asdg\":[1,2,3,4,5]}";
  std::string test_list= "[1,2,3,true,false,\"test\",{\"1\":1,\"2\":2},[1,2,3,4]]";
  // std::string test_numeric_list= "[1,2,3,.1,0.2,-1.3,-Infinity,Infinity,NaN,0.0,-0.0]";
  // should be 12
  std::cout << "Closing brace position: " << drake::examples::kuka_iiwa_arm::find_closing_delim(nested_braces,1) << std::endl;
  // should be 22
  std::cout << "Closing brace position: " << drake::examples::kuka_iiwa_arm::find_closing_delim(nested_brackets,2,'[',']') << std::endl;
  // should be 9
  std::cout << "Closing string position: " << drake::examples::kuka_iiwa_arm::find_closing_delim(test_string,0,'\"','\"') << std::endl;
  // 
  drake::examples::kuka_iiwa_arm::parse_json_object(test_map);
  std::cout << "completed map test" << std::endl;
  drake::examples::kuka_iiwa_arm::parse_json_list(test_list);

  auto kuka = std::make_shared<RigidBodyTree<double>>(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      drake::multibody::joints::kFixed);
  auto lcm = std::make_shared<lcm::LCM>();

  KukaDircolPlanner planner(kuka, lcm);
  planner.run();
  return 0;
}

} // namespace
} // namespace kuka_iiwa_arm
} // namespace examples
} // namespace drake

int main(int argc, const char* argv[]) {
  return drake::examples::kuka_iiwa_arm::main(argc, argv);
}