function [c,dc] = robustVariancecost_check_xonly(obj, x, X0)
x_full = x;%X0(1:obj.nx*obj.N);
Fext_full = X0(obj.nx*obj.N+1:end);

x = reshape(x_full, obj.nx, obj.N);
u = reshape(Fext_full, obj.nFext, obj.N);% note that, in this bricking example, we treat external force as control input
nq = obj.plant.getNumPositions;
nv = obj.plant.getNumVelocities;
nu = obj.nFext;%obj.plant.getNumInputs;

% sigma points
Px = zeros(obj.nx,obj.nx,obj.N);
Px(:,:,1) = obj.cached_Px(:,:,1);

% disturbance variance
% currently only consider terrain height and friction coefficient
Pw = diag([0.01, 0.04]); %[to be tuned]
w_phi = normrnd(zeros(1,obj.N),sqrt(Pw(1,1)),1,obj.N);%height noise
w_mu = normrnd(zeros(1,obj.N),sqrt(Pw(2,2)),1,obj.N);%friction coefficient noise
w_noise = [w_phi;w_mu];

scale = .01;% [to be tuned]
nw = size(Pw,1);
K = [10*ones(nu,nq),ones(nu,nv)];

%initialize c and dc
c = 0;
dc = zeros(1, 1+obj.N*(obj.nx+1));% hand coding number of inputs
kappa = 1;
x_mean = zeros(obj.nx, obj.N);

% mean residual cost at first time step is 0, variance matrix is c(k=1) = Px(1);
%c = kappa*trace(Px(:,:,1));

% initialize gradient of Tr(V) w.r.t state vector x
dTrVdx(:,:,1) = zeros(obj.N-1,obj.nx);
dTrVdu(:,:,1) = zeros(obj.N-1,nu);

tStart = tic;

for k = 1:obj.N-1%[Ye: double check the index]
    %Generate sigma points from Px(i+1)
    %[the sequential way to be modified]
    % currently, only use the initial variance matrix for the
    % propogation
    if k == 1
        [S,d] = chol(blkdiag(Px(:,:,k), Pw), 'lower');
        if d
            diverge  = k;
            return;
        end
        S = scale*S;
        Sig(:,:,k) = [S -S];
        for j = 1:(2*(obj.nx+nw))
            Sig(:,j,k) = Sig(:,j,k) + [x(:,k); w_noise(:,k)];
        end
        w_averg = 1/(2*(obj.nx+nw));
        x_mean(:,k) = zeros(obj.nx,1);
        for j = 1:(2*(obj.nx+nw))
            x_mean(:,k) = x_mean(:,k) + w_averg*Sig(1:obj.nx,j,k);
        end
        c = c + norm(x(:,k)-x_mean(:,k))^2;
    end
    
    %Propagate sigma points through nonlinear dynamics
    for j = 1:(2*(obj.nx+nw))
        % a hacky way to implement the control input
        [H,C,B,dH,dC,dB] = obj.plant.manipulatorDynamics(Sig(1:obj.nx/2,j,k),Sig(obj.nx/2+1:obj.nx,j,k));
        Hinv(:,:,j,k) = inv(H);
        Bmatrix(:,:,j,k) = [1;zeros(5,1)];%B;hand coding
        
        %                     [phi,normal,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(Sig(1:obj.nx/2,j),false,obj.options.active_collision_options);
        %                     J = zeros(nl,nq);
        %                     J(1:1+obj.nD:end,:) = n;
        %                     dJ = zeros(nl*nq,nq);
        %                     dJ(1:1+obj.nD:end,:) = dn;%[double check how dn is factorized]
        %
        %                     for k=1:length(D)
        %                         J(1+k:1+obj.nD:end,:) = D{k};
        %                         dJ(1+k:1+obj.nD:end,:) = dD{k};
        %                     end
        
        % add feedback control
        t = obj.plant.timestep*(k-1);%[double make sure obj.h is updated correctly]
        u_fdb_k = u(:,k) - K*(Sig(1:obj.nx,j,k) - x(:,k));
        [xdn,df] = obj.plant.update(t,Sig(1:obj.nx,j,k),u_fdb_k);
        
        Sig(1:obj.nx/2,j,k+1) = xdn(1:obj.nx/2);
        Sig(obj.nx/2+1:obj.nx,j,k+1) = xdn(obj.nx/2+1:obj.nx);
        dfdu(:,:,j,k+1) = [obj.plant.timestep^2*Hinv(:,:,j,k)*Bmatrix(:,:,j,k);obj.plant.timestep*Hinv(:,:,j,k)*Bmatrix(:,:,j,k)];
        dfdSig(:,:,j,k+1) = df(:,2:obj.nx+1) - dfdu(:,:,j,k+1)*K;
        dfdx(:,:,j,k+1) = dfdu(:,:,j,k+1)*K;
        
        % useless
        %                         if k == 2 % define gradient component when k=1 (initial gradient)
        %                             dfdu(:,:,j,1) = [obj.plant.timestep^2*Hinv(:,:,j,k)*Bmatrix(:,:,j,k);obj.plant.timestep*Hinv(:,:,j,k)*Bmatrix(:,:,j,k)];
        %                             dfdx(:,:,j,1) = dfdu(:,:,j,1)*K;
        %                         end
    end
    
    %Calculate mean and variance w.r.t. [x_k] from sigma points
    w_averg = 1/(2*(obj.nx+nw));
    x_mean(:,k+1) = zeros(obj.nx,1);
    for j = 1:(2*(obj.nx+nw))
        x_mean(:,k+1) = x_mean(:,k+1) + w_averg*Sig(1:obj.nx,j,k+1);
    end
    Px(:,:,k+1) = zeros(obj.nx);
    %alpha = 1e-3;
    %w_coeff = (1/(2*alpha^2*(obj.nx+nw)));
    w = 0.5/scale^2;
    for j = 1:(2*(obj.nx+nw))
        Px(:,:,k+1) = Px(:,:,k+1) + w*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))';
    end
    
    % accumulate returned cost
    c = c + norm(x(:,k+1)-x_mean(:,k+1))^2;
    %                     for j = 1:(2*(obj.nx+nw))
    %                         V_comp = (Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))*(Sig(1:obj.nx,j,k+1)-x_mean(:,k+1))';
    %                         c = c + kappa*trace(w*V_comp);
    %                     end
    
    % derivative of variance matrix
    % gradient of Tr(V) w.r.t state vector x
    dTrVdx(:,:,k+1) = zeros(obj.N-1,obj.nx);
    dTrVdu(:,:,k+1) = zeros(obj.N-1,nu);
    dmeanRdx(:,:,k+1) = zeros(obj.N,obj.nx);
    dmeanRdu(:,:,k+1) = zeros(obj.N-1,nu);
    
    %dmeanRdx(k,:,k) = 2*(x(:,k)-x_mean(:,k))';
    
    for j=k:-1:1
        dTrVdx(j,:,k+1) = zeros(1,obj.nx);
        dTrVu(j,:,k+1) = zeros(1,nu);
        dmeanRdx(j,:,k+1) = zeros(1,obj.nx);
        dmeanRdu(j,:,k+1) = zeros(1,nu);
        
        % gradient w.r.t state x
        dSig_m_kplus1_dx_sum = zeros(obj.nx);
        % gradient w.r.t control u
        dSig_m_kplus1_du_sum = zeros(obj.nx,1);
        %                         dSig_m_kplus1_dSigma1_sum = zeros(obj.nx);
        
        for i=1:2*(obj.nx+nw)
            if i == 1
                for m = 1:(2*(obj.nx+nw))% this for-loop is for \bar{x}_{k+1}, and only needs to go through once since the mean remains the same for different sigma points
                    % gradient of Tr(V_{k+1}) w.r.t control x and u
                    dSig_m_kplus1_dx = zeros(obj.nx);
                    dSig_m_kplus1_du = zeros(obj.nx,1);
                    %                                     dSig_m_kplus1_dSigma1 = zeros(obj.nx);
                    
                    chain_rule_indx = k-j;
                    if j ~= 1
                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1);
                    else% if j == 1, there is an extra gradient to take w.r.t dSigma1_m_dx1 due to sampling mechanism
                        dSig_m_kplus1_dx = dfdx(:,:,m,j+1) + dfdSig(:,:,m,2)*eye(obj.nx);
                    end
                    dSig_m_kplus1_du = dfdu(:,:,m,j+1);% [double check that du is not affected]
                    
                    while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                        dSig_m_kplus1_dx = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_dx;
                        dSig_m_kplus1_du = dfdSig(:,:,m,k+2-chain_rule_indx)*dSig_m_kplus1_du;
                        chain_rule_indx = chain_rule_indx - 1;
                    end
                    dSig_m_kplus1_dx_sum = dSig_m_kplus1_dx_sum+dSig_m_kplus1_dx;
                    dSig_m_kplus1_du_sum = dSig_m_kplus1_du_sum+dSig_m_kplus1_du;
                    
                    %                                     if j == 1% consider a special chain rule w.r.t initial sigma point, i.e., k = 1
                    %                                         dSig_m_kplus1_dSigma1_sum = dSig_m_kplus1_dSigma1_sum + dfdSig(:,:,m,k+1);
                    %                                     end
                end
            end
            
            % run 2*(obj.nx+nw) times in total to obtain
            % gradient w.r.t sigma points
            dSig_i_kplus1_dx = zeros(obj.nx);
            dSig_i_kplus1_du = zeros(obj.nx,1);
            chain_rule_indx = k-j;
            dSig_i_kplus1_dx = dfdx(:,:,i,j+1);
            dSig_i_kplus1_du = dfdu(:,:,i,j+1);
            while(chain_rule_indx>0)% apply the chain rule w.r.t. sigma points
                dSig_i_kplus1_dx = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_dx;
                dSig_i_kplus1_du = dfdSig(:,:,i,k+2-chain_rule_indx)*dSig_i_kplus1_du;
                chain_rule_indx = chain_rule_indx - 1;
            end
            
            dTrVdx(j,:,k+1) = dTrVdx(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_dx - w_averg*dSig_m_kplus1_dx_sum);
            dTrVdu(j,:,k+1) = dTrVdu(j,:,k+1) + 2*w*(Sig(1:obj.nx,i,k+1)-x_mean(:,k+1))'*(dSig_i_kplus1_du - w_averg*dSig_m_kplus1_du_sum);
        end
        
        % gradient of mean residual w.r.t state x and control u, assume norm 2
        dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_averg*dSig_m_kplus1_dx_sum);
        dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 2*(x(:,k+1)-x_mean(:,k+1))'*(-w_averg*dSig_m_kplus1_du_sum);
        
        %dmeanRdx(j,:,k+1) = dmeanRdx(j,:,k+1) + 2*(x(:,k+1))'*(eye(obj.nx));
        %dmeanRdu(j,:,k+1) = dmeanRdu(j,:,k+1) + 0;
    end
end
tElapsed = toc(tStart);

dc = [];
% cost gradient w.r.t x at first time step is zero
for k=1:obj.N % index for x_k
    dTrV_sum_dx_k = zeros(1, obj.nx);
    if (k == 1)
        dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))'*(eye(obj.nx)-eye(obj.nx));% equal to zero vector
    else
        dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))';
    end
    %                     dmeanR_sum_dx_k = zeros(1, obj.nx);
    %                     if (k == 1)
    %                         dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))'*(eye(obj.nx)-eye(obj.nx));
    %                         dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(1,:,2) + dmeanRdx(1,:,3);
    %                     elseif k == 2
    %                         dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))';
    %                         dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(2,:,3);
    %                     else
    %                         dmeanR_sum_dx_k = 2*(x(:,k)-x_mean(:,k))';
    %                     end
    for kk = k+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
        dTrV_sum_dx_k = dTrV_sum_dx_k + dTrVdx(k,:,kk);
        dmeanR_sum_dx_k = dmeanR_sum_dx_k + dmeanRdx(k,:,kk);
    end
    %                     dc = [dc, dmeanR_sum_dx_k+kappa*dTrV_sum_dx_k];
    dc = [dc, dmeanR_sum_dx_k];
end

% cost gradient w.r.t u at first time step is zero, since
% c(k=1) = Px(:,:,1)
% for k=1:obj.N % index for u_k
%     dTrV_sum_du_k = zeros(1, nu);
%     dmeanR_sum_du_k = zeros(1, nu);
%     %                     if k == 1
%     %                     dmeanR_sum_du_k = dmeanRdu(1,:,2)+dmeanRdu(1,:,3);
%     %                     elseif k == 2
%     %                         dmeanR_sum_du_k = dmeanRdu(2,:,3);
%     %                     end
%     for kk = k+1:obj.N % index for TrV_kk and residual of ||x_k - \bar{x}_k||
%         dTrV_sum_du_k = dTrV_sum_du_k + dTrVdu(k,:,kk);
%         dmeanR_sum_du_k = dmeanR_sum_du_k + dmeanRdu(k,:,kk);
%     end
%     %                     dc = [dc, dmeanR_sum_du_k+kappa*dTrV_sum_du_k];
%     dc = [dc, dmeanR_sum_du_k];
% end

dc = dc(:);
obj.cached_Px = Px;
end
