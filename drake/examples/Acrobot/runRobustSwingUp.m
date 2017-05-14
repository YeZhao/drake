function runRobustSwingUp()
%% runs robust trajectory optimization and animates open-loop playback
clear all;clc;

p_perturb = AcrobotPlant;
p_perturb_averg = AcrobotPlant;
p_nominal = AcrobotPlant;
v = AcrobotVisualizer(p_perturb);
v_averg = AcrobotVisualizer(p_perturb_averg);

% --- step 1: generate optimal trajs and LQR gains of nominal model ----
[utraj,xtraj,K] = swingUpTrajectory(p_nominal);

% Synthetic parameter estimation mode
% 'base':       Base case parameter estimation - intial params = true params
% 'paramerr':   w/ parameter error but no measurement noise
% 'measnoise':  w/ parameter error and w/ measurement noise
% 'delay':      w/ parameter error and w/ measurement noise and delay (Not complete)
mode = 'paramerr';

% Standard deviation of the parameter value percent error
paramstd = 1/5;
N = 1; % number of sampled trajectories

% perturb model parameters
for i = 1:N
    p_perturb = p_nominal;
    if ~strcmp(mode,'base')
        % Perturb original parameter estimates with random percentage error
        % normally distributed with standard dev = paramstd, and greater than -1
        paramerr = 0*randn(1,10)*paramstd;
        while sum(paramerr<=-1)~=0
            paramerr(paramerr<-1) = randn(1,sum(paramerr<-1))*paramstd;
        end
        
        % p_perturb.l1 = p_perturb.l1 + p_perturb.l1*paramerr(1);
        % p_perturb.l2 = p_perturb.l2 + p_perturb.l2*paramerr(2);
        % p_perturb.m1 = p_perturb.m1 + p_perturb.m1*paramerr(3);
        % p_perturb.m2 = p_perturb.m2 + p_perturb.m2*paramerr(4);
        
%         p_perturb.b1  = p_perturb.b1 + p_perturb.b1*paramerr(5);
%         p_perturb.b2  = p_perturb.b2 + p_perturb.b2*paramerr(6);
%         p_perturb.lc1 = p_perturb.lc1 + p_perturb.lc1*paramerr(7);
%         p_perturb.lc2 = p_perturb.lc2 + p_perturb.lc2*paramerr(8);
%         p_perturb.Ic1 = p_perturb.Ic1 + p_perturb.Ic1*paramerr(9);
%         p_perturb.Ic2 = p_perturb.Ic2 + p_perturb.Ic2*paramerr(10);
        
%         if i == 1
%             p_perturb_averg.b1  = p_perturb.b1;
%             p_perturb_averg.b2  = p_perturb.b2;
%             p_perturb_averg.lc1 = p_perturb.lc1;
%             p_perturb_averg.lc2 = p_perturb.lc2;
%             p_perturb_averg.Ic1 = p_perturb.Ic1;
%             p_perturb_averg.Ic2 = p_perturb.Ic2;
%         else
%             p_perturb_averg.b1  = p_perturb_averg.b1 + p_perturb.b1;
%             p_perturb_averg.b2  = p_perturb_averg.b2 + p_perturb.b2;
%             p_perturb_averg.lc1 = p_perturb_averg.lc1 + p_perturb.lc1;
%             p_perturb_averg.lc2 = p_perturb_averg.lc2 + p_perturb.lc2;
%             p_perturb_averg.Ic1 = p_perturb_averg.Ic1 + p_perturb.Ic1;
%             p_perturb_averg.Ic2 = p_perturb_averg.Ic2 + p_perturb.Ic2;
%             
%             if i == N
%                 p_perturb_averg.b1  = p_perturb_averg.b1/N;
%                 p_perturb_averg.b2  = p_perturb_averg.b2/N;
%                 p_perturb_averg.lc1 = p_perturb_averg.lc1/N;
%                 p_perturb_averg.lc2 = p_perturb_averg.lc2/N;
%                 p_perturb_averg.Ic1 = p_perturb_averg.Ic1/N;
%                 p_perturb_averg.Ic2 = p_perturb_averg.Ic2/N;
%             end
%         end
    end
    
    [utraj,xtraj] = swingUpTrajectory(p_perturb);
    v = AcrobotVisualizer(p_perturb);
    v.playback(xtraj,struct('slider',true));

    utraj_eval = [];
    xtraj_eval = [];
    utraj_eval = ppval(utraj.pp.breaks,utraj.pp)';
    xtraj_eval = ppval(xtraj.pp.breaks,xtraj.pp)';
        
    figure(1)
    hold on;
    plot(utraj_eval);
    ylabel('u')
    hold on;
    
    figure(2)
    hold on;
    plot(xtraj_eval(:,1));
    ylabel('x_1')
    hold on;
    
    figure(3)
    hold on;
    plot(xtraj_eval(:,2));
    ylabel('x_2')
    hold on;
    
    figure(4)
    hold on;
    plot(xtraj_eval(:,3));
    ylabel('xdot_1')
    hold on;
    
    figure(5)
    hold on;
    plot(xtraj_eval(:,4));
    ylabel('xdot_2')
    hold on;
end

Qf=diag([1000*(1/0.05)^2 1000*(1/0.05)^2 10 10]);
Q = diag([10 10 10 10]);  R=0.1; % LQR Cost Matrices

%[utraj,xtraj] = swingUpTrajectory(p_perturb);

% hacky way to re-implement the average plant
p_perturb.b1  = p_perturb_averg.b1;
p_perturb.b2  = p_perturb_averg.b2;
p_perturb.lc1 = p_perturb_averg.lc1;
p_perturb.lc2 = p_perturb_averg.lc2;
p_perturb.Ic1 = p_perturb_averg.Ic1;
p_perturb.Ic2 = p_perturb_averg.Ic2;

ltvsys = tvlqr(p_perturb,xtraj,utraj,Q,R,Qf);

sys=feedback(p_perturb,ltvsys);

xtraj_new = simulate(sys,xtraj.tspan, [0;0;0;0]);%+0.05*randn(4,1)
v = AcrobotVisualizer(p_perturb);

v.playback(xtraj_new,struct('slider',true));

% xtraj_new_eval = ppval(xtraj_new.pp.breaks,xtraj_new.pp)';
% 
% figure(2)
% hold on;
% plot(xtraj_new_eval(:,1));
% ylabel('x_1')
% hold on;
% 
% figure(3)
% hold on;
% plot(xtraj_new_eval(:,2));
% ylabel('x_2')
% hold on;
% 
% figure(4)
% hold on;
% plot(xtraj_new_eval(:,3));
% ylabel('xdot_1')
% hold on;
% 
% figure(5)
% hold on;
% plot(xtraj_new_eval(:,4));
% ylabel('xdot_2')
% hold on;

%v.playback(xtraj);

end