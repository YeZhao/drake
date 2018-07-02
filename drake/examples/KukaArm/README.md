# run robust traj opt demo

This folder contains a demo guidance to run robust traj opt and several debugging comments.

## Pick and place demo

This demo shows a pick and place motion.
 * run runTrajOpt.m and tune several key parameters
    * number of knot points N.
    * specify a uncertainty type r.uncertainty_source: 'friction_coeff', 'object_initial_position', 'friction_coeff+object_initial_position'. (Currently set to be empty, no uncertainty)
    * robust cost parameter contact_robust_cost_coeff (for the arm example, usually 1 is the best, 10 causes 41 error, 0.1 solves successfully but causes nonphysical motions)
    * total time: T0 (usually 2-4seconds)
    * Traj Opt tolerances: 'MajorFeasibilityTolerance', 'MinorFeasibilityTolerance', 'MinorOptimalityTolerance', 'MajorOptimalityTolerance' (usually [1e-6,5e-3])
    * Cost coeff in running cost running_cost_fun (currently use R = 1e-6, Q = blkdiag(10*eye(8),0*eye(6),10*eye(14)), state term penalizes deviation from final state).

 * robust cost function robustVariancecost_ML() is defined in RobustContactImplicitTrajectoryOptimization_Kuka.m
    * Select different sampling mechanism: 1: unscented transform, option 2: random sampling with a smaller number of sample point (currently, still test single sample point on the nominal trajectory) [Line578-586](https://github.com/YeZhao/drake/blob/3ead1b39928b2c8270c6fe1f88193642dfff8e21/drake/matlab/solvers/trajectoryOptimization/RobustContactImplicitTrajectoryOptimization_Kuka.m#L578-L586).
   * noise_sample_type: 1: each sample point has one specific environment instance; 2: each sample point generates M next-time-step sample points via propagating through all environment instances. (currently, use 1 for computation efficiency) [Line609](https://github.com/YeZhao/drake/blob/3ead1b39928b2c8270c6fe1f88193642dfff8e21/drake/matlab/solvers/trajectoryOptimization/RobustContactImplicitTrajectoryOptimization_Kuka.m#L609).
   * Debugging for QP solveer numerical gradient: [Line691-712](https://github.com/YeZhao/drake/blob/3ead1b39928b2c8270c6fe1f88193642dfff8e21/drake/matlab/solvers/trajectoryOptimization/RobustContactImplicitTrajectoryOptimization_Kuka.m#L691-L712).

 * QP solver in TimeSteppingRigidBodyManipulator_Kuka.m
    * Select the QP regularizer R matrix on [Line842](https://github.com/YeZhao/drake/blob/3ead1b39928b2c8270c6fe1f88193642dfff8e21/drake/matlab/systems/plants/TimeSteppingRigidBodyManipulator_Kuka.m#L842),: 1: nonlinear scaling function, 2: linear scaling function, 3: constant, e.g., 1e-3 (currently use constant one) 
    * Set a constant term in Q matrix, currently use 1e-6*eye(num_params). See [Line864](https://github.com/YeZhao/drake/blob/3ead1b39928b2c8270c6fe1f88193642dfff8e21/drake/matlab/systems/plants/TimeSteppingRigidBodyManipulator_Kuka.m#L864).
    * Set v_min = 0.05/h (reduce gradient sensitivity. used to be phi(i)/h, i is the index of contact point). See
    [Line797](https://github.com/YeZhao/drake/blob/3ead1b39928b2c8270c6fe1f88193642dfff8e21/drake/matlab/systems/plants/TimeSteppingRigidBodyManipulator_Kuka.m#L797),


## Comments

After fixing out several tricks and bugs in QP solver, I realize that the numerical and analytical gradients match with each other to an satisfactorily accurate level. Now all the traj opts are run merely with analytical gradients, which is fast. There is no need to use numerical graidents, which are quite slow. 

Number of knot points should be at least 20. If it is fewer than that, it is likely to generate non-physical motion even though it gives a success sovle.

There is a second test only involving first phase motion: pick/grasping (no A-to-B motion). It is in runTrajOpt2.m 

Run original Michael’s implementation. In RobustContactImplicitTrajectoryOptimization_Kuka.m file, comment out Line 243 and uncomment the following lines: Lines 130 (dynamics constraint), Lines 141-144 (two nonlinear complementarity constraints), Lines 177-178 (linear complementarity constraint), Lines 245-247 (complementarity constraint slack variable cost). Then it should give a non-robust run.