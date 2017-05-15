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

% plot the nominal model trajs
figure(1)
hold on;
plot(u_nominal);
ylabel('u')
hold on;

figure(2)
subplot(2,2,1)
hold on;
plot(x_nominal(1,:));
ylabel('x_1')
hold on;

subplot(2,2,2)
hold on;
plot(x_nominal(2,:));
ylabel('x_2')
hold on;

subplot(2,2,3)
hold on;
plot(x_nominal(3,:));
ylabel('xdot_1')
hold on;

subplot(2,2,4)
hold on;
plot(x_nominal(4,:));
ylabel('xdot_2')
hold on;

% --- step 2: generate optimal trajs of perturbed model ----

% Synthetic parameter estimation mode
% 'base':       Base case parameter estimation - intial params = true params
% 'paramerr':   w/ parameter error but no measurement noise
% 'measnoise':  w/ parameter error and w/ measurement noise
% 'delay':      w/ parameter error and w/ measurement noise and delay (Not complete)
mode = 'paramerr';

% Standard deviation of the parameter value percent error
paramstd = 1/5;
SampleNum = 1; % number of sampled trajectories
Qr = diag([10 10 1 1]);
Rr = .1;
Qrf = 100*eye(4);
color_line_type1 = 'r--';
color_line_type2 = 'g--';
color_line_type3 = 'k--';
color_line_type4 = 'b-.';
color_line_type5 = 'r-.';

% perturb model parameters
for i = 1:SampleNum
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
        p_perturb.m1 = p_perturb.m1 + p_perturb.m1*paramerr(3);
        p_perturb.m2 = p_perturb.m2 + p_perturb.m2*paramerr(4);
        % p_perturb.b1  = p_perturb.b1 + p_perturb.b1*paramerr(5);
        % p_perturb.b2  = p_perturb.b2 + p_perturb.b2*paramerr(6);
        % p_perturb.lc1 = p_perturb.lc1 + p_perturb.lc1*paramerr(7);
        % p_perturb.lc2 = p_perturb.lc2 + p_perturb.lc2*paramerr(8);
        % p_perturb.Ic1 = p_perturb.Ic1 + p_perturb.Ic1*paramerr(9);
        % p_perturb.Ic2 = p_perturb.Ic2 + p_perturb.Ic2*paramerr(10);
        
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
    
    [utraj,xtraj,z,prog] = robustswinguptrajectory(p_perturb,u_nominal,x_nominal,K_nominal,Qr,Qrf,Rr,N);
    v_perturb = AcrobotVisualizer(p_perturb);
    v_perturb.playback(xtraj,struct('slider',true));

    utraj_eval = [];
    xtraj_eval = [];
    utraj_eval = ppval(utraj.pp.breaks,utraj.pp)';
    xtraj_eval = ppval(xtraj.pp.breaks,xtraj.pp)';
        
    figure(1)
    hold on;
    plot(utraj_eval,color_line_type3);
    ylabel('u')
    hold on;
    
    figure(2)
    subplot(2,2,1)
    hold on;
    plot(xtraj_eval(:,1),color_line_type3);
    ylabel('x_1')
    hold on;
    
    subplot(2,2,2)
    hold on;
    plot(xtraj_eval(:,2),color_line_type3);
    ylabel('x_2')
    hold on;
    
    subplot(2,2,3)
    hold on;
    plot(xtraj_eval(:,3),color_line_type3);
    ylabel('xdot_1')
    hold on;
    
    subplot(2,2,4)
    hold on;
    plot(xtraj_eval(:,4),color_line_type3);
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