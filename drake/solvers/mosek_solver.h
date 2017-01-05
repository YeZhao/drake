// A wrapper file for MosekLP and mosekQP that handles constraint and
// objective marshalling
#pragma once

#include <Eigen/Core>

#include "drake/solvers/mathematical_program.h"

namespace drake {
namespace solvers {

class MosekSolver : public MathematicalProgramSolverInterface {
 public:
  /** available()
  * Defined true if Mosek was included during compilation, false otherwise.
  */
  bool available() const override;

  SolutionResult Solve(MathematicalProgram& prog) const override;
};

}  // namespace solvers
}  // namespace drake
