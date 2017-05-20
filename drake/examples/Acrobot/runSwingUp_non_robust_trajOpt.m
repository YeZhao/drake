function runSwingUp()
%% runs trajectory optimization and animates open-loop playback

p_nominal = AcrobotPlant;
v = AcrobotVisualizer(p_nominal);
[utraj,xtraj,z,prog] = swingUpTrajectory(p_nominal);
%      sys = cascade(utraj,p);
%      xtraj=simulate(sys,utraj.tspan,zeros(4,1));
v.playback(xtraj);

h_nominal = z(prog.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj.eval(t_nominal);% this is exactly same as z_nominal components
u_nominal = utraj.eval(t_nominal)';

color_line_type_set = {'r-','g-','k-','m-','c-','b-.','r-.','g-.','k-.','m-.'};
m = 1;
nominal_linewidth = 1;

color_line_type = 'b-';
figure(10)
hold on;
plot(t_nominal', u_nominal, color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('u');
hold on;

figure(21)
subplot(2,2,1)
hold on;
plot(t_nominal', x_nominal(1,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('x_1')
hold on;

subplot(2,2,2)
hold on;
plot(t_nominal', x_nominal(2,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('x_2')
hold on;

subplot(2,2,3)
hold on;
plot(t_nominal', x_nominal(3,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('xdot_1')
hold on;

subplot(2,2,4)
hold on;
plot(t_nominal', x_nominal(4,:), color_line_type, 'LineWidth',nominal_linewidth);
xlabel('t');
ylabel('xdot_2')
hold on;

% LQR Cost Matrices
Q = diag([10 10 1 1]);
R = .1;
Qf = 100*eye(4);

% perturb model parameters
% run 100 times simulations to test the success rate
TestNum = 100;
NumofSuccessTest = 0;
paramerr_test = [];
paramstd = 1/5;
count = 0;
paramerr_test = load('PARAMERR2.dat');
fail_index = [];
residual = [];
for i = 1:TestNum
    p_perturb = p_nominal;
    ltvsys = tvlqr(p_perturb,xtraj,utraj,Q,R,Qf);

    % Perturb original parameter estimates with random percentage error
    % normally distributed with standard dev = paramstd, and greater than -1
%     paramerr_test(i,:) = randn(1,10)*paramstd;
%     while sum(paramerr_test(i,:)<=-1)~=0
%         paramerr_test(paramerr_test(i,:)<-1) = randn(1,sum(paramerr_test(i,:)<-1))*paramstd;
%     end     
    p_perturb.m1 = p_perturb.m1 + p_perturb.m1*paramerr_test(i,3);
    p_perturb.m2 = p_perturb.m2 + p_perturb.m2*paramerr_test(i,4);
    
    p_perturb.l1 = p_perturb.l1 + p_perturb.l1*paramerr_test(i,1);
    p_perturb.l2 = p_perturb.l2 + p_perturb.l2*paramerr_test(i,2);
        
    p_perturb.lc1 = p_perturb.lc1 + p_perturb.lc1*paramerr_test(i,7);
    p_perturb.lc2 = p_perturb.lc2 + p_perturb.lc2*paramerr_test(i,8);
    p_perturb.Ic1 = p_perturb.Ic1 + p_perturb.Ic1*paramerr_test(i,9);
    p_perturb.Ic2 = p_perturb.Ic2 + p_perturb.Ic2*paramerr_test(i,10);
    
    sys=feedback(p_perturb,ltvsys);
    
    xtraj_new = simulate(sys,xtraj.tspan, [0;0;0;0]);%+0.05*randn(4,1)
%     v = AcrobotVisualizer(p_perturb);
%     v.playback(xtraj_new,struct('slider',true));
    
    xtraj_new_eval = [];
    t_new_eval = [];
    
    fieldname_str = strcmp(fieldnames(xtraj_new),'traj');
    if fieldname_str(1) == 1
        for j = 1:size(xtraj_new.traj,2)
            t_new_eval = [t_new_eval, linspace(xtraj_new.traj{j}.tspan(1),xtraj_new.traj{j}.tspan(2),size(xtraj_new.traj{j}.pp.breaks,2))];
            xtraj_new_eval = [xtraj_new_eval, ppval(xtraj_new.traj{j}.pp.breaks,xtraj_new.traj{j}.pp)];
        end
    else
        t_new_eval = linspace(xtraj_new.tspan(1),xtraj_new.tspan(2),size(xtraj_new.pp.breaks,2));
        xtraj_new_eval = ppval(xtraj_new.pp.breaks,xtraj_new.pp);
    end

    figure(21)
    subplot(2,2,1)
    hold on;
    plot(t_new_eval, xtraj_new_eval(1,:),color_line_type_set{m}, 'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('x_1')
    hold on;
    
    subplot(2,2,2)
    hold on;
    plot(t_new_eval, xtraj_new_eval(2,:),color_line_type_set{m},'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('x_2')
    hold on;
    
    subplot(2,2,3)
    hold on;
    plot(t_new_eval, xtraj_new_eval(3,:),color_line_type_set{m},'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('xdot_1')
    hold on;
    
    subplot(2,2,4)
    hold on;
    plot(t_new_eval, xtraj_new_eval(4,:),color_line_type_set{m},'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('xdot_2')
    hold on;
    
    for j = 1:size(p_nominal.xG,1)
        x_final_desired(j) = p_nominal.xG(j);
    end
    
    x_final_diff = xtraj_new_eval(:,end) - x_final_desired';
    
    if (sum(abs(x_final_diff)) < 1)
        NumofSuccessTest = NumofSuccessTest + 1;
    else
        disp('fail this time');
        fail_index = [fail_index, count];
    end
    residual = [residual sum(abs(x_final_diff))];
    
    count
    i
    NumofSuccessTest
    count = count + 1;
end

Q = diag([10 10 1 1]);
R = .1;
Qf = 100*eye(4);

ltvsys = tvlqr(p,xtraj,utraj,Q,R,Qf);
sys=feedback(p,ltvsys);
xtraj_new = simulate(sys,xtraj.tspan, [0;0;0;0]);%+0.05*randn(4,1)


v = AcrobotVisualizer(p);
v.playback(xtraj_new,struct('slider',true));
end
