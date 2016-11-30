#include <Eigen/Core>

#include "kuka_plant.h"
#include "drake/common/eigen_autodiff_types.h"
#include "drake/systems/plants/constraint/direct_collocation_constraint.h"

namespace drake{
namespace examples{
namespace kuka_iiwa_arm{
namespace{

class KukaRunningCost{
  public:
    // TODO: take this from the arm get_num_positions
    static size_t numInputs(){return 14;}
    static size_t numOutputs(){return 1;} 

    template <typename ScalarType>
    void eval(const VecIn<ScalarType>& x, VecOut<ScalarType>& y) const {
      y = x.transpose()*x;
    }
}

class KukaFinalCost{
  public:
    KukaFinalCost(Eigen::VectorXd goal): goal(goal) {}
  // TODO: take this from the arm get_num_positions
  static size_t numInputs(){return 14;}
  static size_t numOutputs(){return 1;} 

  template <typename ScalarType>
  void eval(const VecIn<ScalarType>& x, VecOut<ScalarType>& y) const {
    y = (x-goal).transpose()*(x-goal);
  }
  private:
    Eigen::VectorXd goal;
}

void AddTrajectoryParams(int num_time_samples, const Eigen::VectorXd& x0 const Eigen::VectorXd& xG,
                         solvers::DircolTrajectoryOptimization* dircol_traj, std::unique_ptr<RigidBodyTree<AutoDiffXd>> tree){
  const int kTorqueLimit = 3;
  // TODO: change 7 to get_num_positions
  const Eigen::VectorXd u_min = Eigen::VectorXd::Constant(7,-kTorqueLimit);
  const Eigen::VectorXd u_max = Eigen::VectorXd::Constant(7,kTorqueLimit);

  dircol_traj->AddInputBounds(u_min, u_max);

  dircol_traj->AddStateConstraint(
      std::make_shared<LinearEqualityConstraint>(
          Eigen::MatrixXd::Identity(14), x0), {0});
  dircol_traj->AddStateConstraint(
      std::make_shared<LinearEqualityConstraint>(
          Eigen::MatrixXd::Identity(14), xG), {num_time_samples - 1});

  dircol_traj->AddRunningCostFunc(KukaRunningCost());
  dircol_traj->AddFinalCostFunc(KukaFinalCost(xG));
  dircol_traj->AddDynamicConstraint(
      std::make_shared<
      drake::systems::System2DirectCollocationConstraint>(
          std::make_unique<RigidBodyPlant>(tree)));
}
}

} // anonymous
} // kuka_iiwa_arm
} // examples
} // drake