#include <iostream>
#include <memory>

#include "kuka_lcm_planner.h"
#include "drake/common/drake_path.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/multibody/parsers/urdf_parser.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

int main(int argc, const char* argv[]) {

  auto kuka = std::make_shared<RigidBodyTree<double>>();
  parsers::urdf::AddModelInstanceFromUrdfFileToWorld(
      GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      multibody::joints::kFixed, kuka.get());
  auto lcm = std::make_shared<lcm::LCM>();

  //(TODO: fix the kuka RigidBodyTree reference issue)  
  KukaMatlabDircolPlanner planner(&*kuka, lcm);
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