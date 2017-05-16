function runRobustSwingUp()
%% runs robust trajectory optimization and animates open-loop playback
clear all;clc;

p_perturb = AcrobotPlant;
p_perturb_averg = AcrobotPlant;
p_nominal = AcrobotPlant;
v = AcrobotVisualizer(p_perturb);
v_averg = AcrobotVisualizer(p_perturb_averg);
N = 41;

% --- step 1: generate optimal trajs and LQR gains of nominal model ----
[utraj_nominal,xtraj_nominal,z_nominal,prog_nominal,K_nominal] = swingUpTrajectory(p_nominal,N);
v_nominal = AcrobotVisualizer(p_nominal);
v_nominal.playback(xtraj_nominal,struct('slider',true));

h_nominal = z_nominal(prog_nominal.h_inds);
t_nominal = [0; cumsum(h_nominal)];
x_nominal = xtraj_nominal.eval(t_nominal);% this is exactly same as z_nominal components
u_nominal = utraj_nominal.eval(t_nominal)';

% plot nominal model trajs
nominal_linewidth = 2.5;
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

% simulate with LQR gains
% LQR Cost Matrices
Q = diag([10 10 1 1]);
R = .1;
Qf = 100*eye(4);

ltvsys = tvlqr(p_nominal,xtraj_nominal,utraj_nominal,Q,R,Qf);

sys=feedback(p_nominal,ltvsys);

xtraj_new = simulate(sys,xtraj_nominal.tspan, [0;0;0;0]);%+0.05*randn(4,1)
v = AcrobotVisualizer(p_nominal);
v.playback(xtraj_new,struct('slider',true));

% --- step 2: generate optimal trajs of perturbed model ----

% Synthetic parameter estimation mode
% 'base':       Base case parameter estimation - intial params = true params
% 'paramerr':   w/ parameter error but no measurement noise
% 'measnoise':  w/ parameter error and w/ measurement noise
% 'delay':      w/ parameter error and w/ measurement noise and delay (Not complete)
mode = 'paramerr';

paramstd = 1/10; % Standard deviation of the parameter value percent error
SampleNum = 20; % number of sampled trajectories
Qr = diag([10 10 1 1]);
Rr = .1;
Qrf = 100*eye(4);
color_line_type_set = {'r-','g-','k-','m-','c-','b-.','r-.','g-.','k-.','m-.'};

% perturb model parameters
paramerr = [];
for i = 1:SampleNum
    if ~strcmp(mode,'base')
        % Perturb original parameter estimates with random percentage error
        % normally distributed with standard dev = paramstd, and greater than -1
        paramerr(i,:) = randn(1,10)*paramstd;
        while sum(paramerr(i,:)<=-1)~=0
            paramerr(paramerr(i,:)<-1) = randn(1,sum(paramerr(i,:)<-1))*paramstd;
        end        
    end    
end

Max_iter = 10;
m = 1;
x_traj_max_diff_sum_percent = 100;
x_nominal_new = [];
u_nominal_new = [];
K_nominal_new = [];

while(m <= Max_iter && x_traj_max_diff_sum_percent > 10)
    
    utraj_eval = [];
    xtraj_eval = [];
    utrajArray = [];
    xtrajArray = [];
    if m > 1
        u_nominal = u_nominal_new;
        x_nominal = x_nominal_new;
        K_nominal = K_nominal_new;
    end
    
    for i = 1:SampleNum
        p_perturb = p_nominal;
        if ~strcmp(mode,'base')
            % Perturb original parameter estimates with random percentage error
            
            % p_perturb.l1 = p_perturb.l1 + p_perturb.l1*paramerr(i,1);
            % p_perturb.l2 = p_perturb.l2 + p_perturb.l2*paramerr(i,2);
            p_perturb.m1 = p_perturb.m1 + p_perturb.m1*paramerr(i,3);
            p_perturb.m2 = p_perturb.m2 + p_perturb.m2*paramerr(i,4);
            % p_perturb.b1  = p_perturb.b1 + p_perturb.b1*paramerr(i,5);
            % p_perturb.b2  = p_perturb.b2 + p_perturb.b2*paramerr(i,6);
            % p_perturb.lc1 = p_perturb.lc1 + p_perturb.lc1*paramerr(i,7);
            % p_perturb.lc2 = p_perturb.lc2 + p_perturb.lc2*paramerr(i,8);
            % p_perturb.Ic1 = p_perturb.Ic1 + p_perturb.Ic1*paramerr(i,9);
            % p_perturb.Ic2 = p_perturb.Ic2 + p_perturb.Ic2*paramerr(i,10);
        end
        
        [utraj,xtraj,z_perturb,prog_perturb] = robustSwingUpTrajectory(p_perturb,u_nominal,x_nominal,K_nominal,Qr,Qrf,Rr,N);
        v_perturb = AcrobotVisualizer(p_perturb);
        v_perturb.playback(xtraj,struct('slider',true));
        
        utraj_eval = [];
        xtraj_eval = [];
        h_perturb = z_perturb(prog_perturb.h_inds);
        t_perturb = [0; cumsum(h_perturb)];
        xtraj_eval = xtraj.eval(t_perturb);% this is exactly same as z_nominal components
        utraj_eval = utraj.eval(t_perturb)';

        utrajArray(:,:,i) = utraj_eval;
        xtrajArray(:,:,i) = xtraj_eval;
        
%         figure(10)
%         hold on;
%         plot(t_perturb', utraj_eval,color_line_type_set{1});
%         xlabel('t');
%         ylabel('u')
%         hold on;
%         
%         figure(20)
%         subplot(2,2,1)
%         hold on;
%         plot(t_perturb', xtraj_eval(1,:),color_line_type_set{1});
%         xlabel('t');
%         ylabel('x_1')
%         hold on;
%         
%         subplot(2,2,2)
%         hold on;
%         plot(t_perturb', xtraj_eval(2,:),color_line_type_set{1});
%         xlabel('t');
%         ylabel('x_2')
%         hold on;
%         
%         subplot(2,2,3)
%         hold on;
%         plot(t_perturb', xtraj_eval(3,:),color_line_type_set{1});
%         xlabel('t');
%         ylabel('xdot_1')
%         hold on;
%         
%         subplot(2,2,4)
%         hold on;
%         plot(t_perturb', xtraj_eval(4,:),color_line_type_set{1});
%         xlabel('t');
%         ylabel('xdot_2')
%         hold on;
    end
    
    % --- step 3: regenerate optimal trajs of nominal model ----
    
    [utraj_nominal_new,xtraj_nominal_new,z_nominal_new,prog_nominal_new,K_nominal_new] = swingUpTrajectory(p_nominal,N,utrajArray,xtrajArray);
    v_nominal_new = AcrobotVisualizer(p_nominal);
    v_nominal_new.playback(xtraj_nominal_new,struct('slider',true));
    
    h_nominal_new = z_nominal_new(prog_nominal_new.h_inds);
    t_nominal_new = [0; cumsum(h_nominal_new)];
    x_nominal_new = xtraj_nominal_new.eval(t_nominal_new);% this is exactly same as z_nominal components
    u_nominal_new = utraj_nominal_new.eval(t_nominal_new)';
    
    % plot the nominal model trajs
    figure(11)
    hold on;
    plot(t_nominal_new', u_nominal_new,color_line_type_set{m}, 'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('u')
    hold on;
    
    figure(21)
    subplot(2,2,1)
    hold on;
    plot(t_nominal_new', x_nominal_new(1,:),color_line_type_set{m}, 'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('x_1')
    hold on;
    
    subplot(2,2,2)
    hold on;
    plot(t_nominal_new', x_nominal_new(2,:),color_line_type_set{m},'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('x_2')
    hold on;
    
    subplot(2,2,3)
    hold on;
    plot(t_nominal_new', x_nominal_new(3,:),color_line_type_set{m},'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('xdot_1')
    hold on;
    
    subplot(2,2,4)
    hold on;
    plot(t_nominal_new', x_nominal_new(4,:),color_line_type_set{m},'LineWidth',nominal_linewidth);
    xlabel('t');
    ylabel('xdot_2')
    hold on;
    %legend('initial nominal','1st iterative nominal','2nd iterative nominal','3rd iterative nominal','4th iterative nominal','5th iterative nominal','6th iterative nominal',);

    % termination condition
    % accumulate state traj changes
    x_traj_diff = abs((x_nominal - x_nominal_new)./x_nominal);
    [row, col] = find(isnan(x_traj_diff));% set nan elements to zero
    x_traj_diff(row,col) = 0;
    [row, col] = find(isinf(x_traj_diff));% set inf elements to zero
    x_traj_diff(row,col) = 0;
    
    x_traj_diff_sum = sum(x_traj_diff');
    x_traj_diff_sum_percent = 100*x_traj_diff_sum/N;
    x_traj_max_diff_sum_percent = max(x_traj_diff_sum_percent);
    
    % increment interation index
    m = m + 1;
    
    % save data of new nominal model
    save U_NOMINAL_BEST.dat u_nominal_new -ASCII
    save X_NOMINAL_BEST.dat x_nominal_new -ASCII
    save K_NOMINAL_BEST.dat K_nominal_new -ASCII

    save XTRAJ_PPTRAJ.mat xtraj_nominal_new
    save UTRAJ_PPTRAJ.mat utraj_nominal_new
    
%     u_nominal_new = load('U_NOMINAL_BEST.dat'); 
%     x_nominal_new = load('X_NOMINAL_BEST.dat'); 
%     K_nominal_new = load('K_NOMINAL_BEST.dat'); 
%     
%     xtraj_pptraj = load('XTRAJ_PPTRAJ.mat');
%     utraj_pptraj = load('UTRAJ_PPTRAJ.mat');
end

% LQR Cost Matrices
Q = diag([10 10 1 1]);
R = .1;
Qf = 100*eye(4);

% perturb model parameters
% run 100 times simulations to test the success rate
TestNum = 100;
NumofSuccessTest = 0;
paramerr_test = [];
for i = 1:TestNum
    p_perturb = p_nominal;
    if ~strcmp(mode,'base')
        % Perturb original parameter estimates with random percentage error
        % normally distributed with standard dev = paramstd, and greater than -1
        paramerr_test(i,:) = randn(1,10)*paramstd;
        while sum(paramerr_test(i,:)<=-1)~=0
            paramerr_test(paramerr_test(i,:)<-1) = randn(1,sum(paramerr_test(i,:)<-1))*paramstd;
        end        
        p_perturb.m1 = p_perturb.m1 + p_perturb.m1*paramerr_test(i,3);
        p_perturb.m2 = p_perturb.m2 + p_perturb.m2*paramerr_test(i,4);
    end
    
    ltvsys = tvlqr(p_perturb,xtraj_nominal_new,utraj_nominal_new,Q,R,Qf);
    
    sys=feedback(p_perturb,ltvsys);
    
    xtraj_new = simulate(sys,xtraj_nominal_new.tspan, [0;0;0;0]);%+0.05*randn(4,1)
    v = AcrobotVisualizer(p_perturb);
    v.playback(xtraj_new,struct('slider',true));
    
    xtraj_new_eval = [];
    t_new_eval = [];
    for j = 1:size(xtraj_new.traj,2)
        t_new_eval = [t_new_eval, linspace(xtraj_new.traj{j}.tspan(1),xtraj_new.traj{j}.tspan(2),size(xtraj_new.traj{j}.pp.breaks,2))];
        xtraj_new_eval = [xtraj_new_eval, ppval(xtraj_new.traj{j}.pp.breaks,xtraj_new.traj{j}.pp)];
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
    
    for i = 1:size(p_nominal.xG,1)
        x_final_desired(i) = p_nominal.xG(i);
    end
    
    x_final_diff = xtraj_new_eval(:,end) - x_final_desired';
    
    if (sum(abs(x_final_diff)) < 0.01)
        NumofSuccessTest = NumofSuccessTest + 1;
    end
    
end

end