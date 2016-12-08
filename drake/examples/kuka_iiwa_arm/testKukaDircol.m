x0 = zeros(14,1);
xG = [0;-0.68;0;-1.688;0;0.5635;0;zeros(7,1)];

planner = KukaPlanner();
[t, x, u, info] = planner.plan_dircol(x0, xG);