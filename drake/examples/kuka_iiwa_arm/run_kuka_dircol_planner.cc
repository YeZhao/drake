#import <iostream>

#import "kuka_lcm_planner.h"

namespace drake {
namespace examples {
namespace kuka_iiwa_arm {
namespace {

int main(int argc, const char* argv[]) {
  std::cout << "Test1 \n";
  KukaDircolPlanner planner;
  std::cout << "Test\n";
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