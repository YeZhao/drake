#include <iostream>
#include <memory>

#include "kuka_lcm_planner.h"
#include "drake/common/drake_path.h"
#include "drake/multibody/rigid_body_tree.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

int main(int argc, const char* argv[]) {

  RigidBodyTree<double> kuka(
      drake::GetDrakePath() + "/examples/kuka_iiwa_arm/urdf/iiwa14.urdf",
      drake::multibody::joints::kFixed);
  auto lcm = std::make_shared<lcm::LCM>();

  KukaMatlabDircolPlanner planner(&kuka, lcm);
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