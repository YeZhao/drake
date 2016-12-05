// standard libraries
#include <iostream>
#include <string>
#include <memory>
// externals
#include <Eigen/Core>
// drake
#include "drake/common/drake_path.h"
#include "drake/multibody/rigid_body_tree.h"
#include "drake/multibody/rigid_body_plant/rigid_body_plant.h"
#include "drake/multibody/rigid_body_plant/contact_results.h"
#include "drake/systems/plants/constraint/direct_collocation_constraint.h"
#include "drake/common/eigen_autodiff_types.h"


namespace drake{
namespace examples{
namespace kuka_iiwa_arm{
namespace{
using drake::systems::RigidBodyPlant;
using Eigen::AutoDiffScalar;
using Eigen::VectorXd;

int main(){
  auto urdf_path = drake::GetDrakePath()+"/examples/kuka_iiwa_arm/urdf/iiwa14.urdf";
  AutoDiffScalar<VectorXd> test(10);
  // RigidBodyPlant<AutoDiffScalar<VectorXd>> kuka_plant(std::make_unique<RigidBodyTree<AutoDiffScalar<VectorXd>>>(urdf_path, drake::multibody::joints::kFixed));


  const Eigen::VectorXd x0 = Eigen::VectorXd::Zero(14);
  const Eigen::VectorXd xG = Eigen::VectorXd::Constant(14,0.5);

  return 0;
}

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake

int main(){
  return drake::examples::kuka_iiwa_arm::main();
}