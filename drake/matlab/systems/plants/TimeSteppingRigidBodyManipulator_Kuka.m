classdef TimeSteppingRigidBodyManipulator_Kuka < DrakeSystem
    % A discrete time system which simulates (an Euler approximation of) the
    % manipulator equations, with contact / limits resolved using the linear
    % complementarity problem formulation of contact in Stewart96.
    
    properties (Access=protected)
        manip  % the CT manipulator
        sensor % additional TimeSteppingRigidBodySensors (beyond the sensors attached to manip)
        dirty=true;
    end
    
    properties (SetAccess=protected)
        timestep 
        twoD=false
        position_control=false;
        LCP_cache;
        enable_fastqp; % whether we use the active set LCP
        lcmgl_contact_forces_scale = 0;  % <=0 implies no lcmgl
        z_inactive_guess_tol = .01;
        multiple_contacts = false;
        gurobi_present = false;
        update_convex = true;
        phi_max = 0.1; % m, max contact force distance %this threhold should be small
        phiL_max = 0.05; % m, max contact force distance %this threhold should be different than phi_max
        active_threshold = 0.1; % height below which contact forces are calculated
        contact_threshold = 1e-3; % threshold where force penalties are eliminated (modulo regularization)
        active_collision_options; % used in contactConstraint
        body
    end
     
    methods
        function obj=TimeSteppingRigidBodyManipulator_Kuka(manipulator_or_urdf_filename,timestep,options)
            if (nargin<3) options=struct(); end
            if ~isfield(options,'twoD') options.twoD = false; end
            
            if ispc
                error('Drake:MissingDependency:PathLCP','PathLCP is known to fail on windows.  See https://github.com/RobotLocomotion/drake/issues/569');
            end
            
            typecheck(timestep,'double');
            sizecheck(timestep,1);
            
            if isempty(manipulator_or_urdf_filename) || ischar(manipulator_or_urdf_filename)
                % then make the corresponding manipulator
                w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
                if options.twoD
                    manip = PlanarRigidBodyManipulator(manipulator_or_urdf_filename,options);
                else
                    manip = RigidBodyManipulator(manipulator_or_urdf_filename,options);
                end
                warning(w);
            else
                manip = manipulator_or_urdf_filename;
            end
            
            typecheck(manip,'RigidBodyManipulator');
            obj = obj@DrakeSystem(0,manip.getNumStates(),manip.getNumInputs(),manip.getNumOutputs(),manip.isDirectFeedthrough(),manip.isTI());
            obj.manip = manip;
            if isa(manip,'PlanarRigidBodyManipulator')
                obj.twoD = true;
            end
            
            if isfield(options, 'multiple_contacts')
                typecheck(options.multiple_contacts, 'logical');
                obj.multiple_contacts = options.multiple_contacts;
            end
            
            if ~isfield(options,'enable_fastqp')
                obj.enable_fastqp = checkDependency('fastqp');
            else
                typecheck(options.enable_fastqp,'logical');
                obj.enable_fastqp = options.enable_fastqp;
                if obj.enable_fastqp && ~checkDependency('fastqp')
                    warning('Drake:TimeSteppingRigidBodyManipulator_Kuka:MissingDependency','You seem to be missing fastQP. Disabling active-set LCP update.')
                    obj.enable_fastqp = false;
                end
            end
            
            if isfield(options,'z_inactive_guess_tol')
                % you might consider setting this if the system consistently
                % gives the "ResolvingLCP" warning
                obj.z_inactive_guess_tol = options.z_inactive_guess_tol;
            end
            
            if isfield(options,'lcmgl_contact_forces_scale')
                obj.lcmgl_contact_forces_scale = options.lcmgl_contact_forces_scale;
                if obj.lcmgl_contact_forces_scale>0,
                    checkDependency('lcmgl');
                end
            end
            
            obj.timestep = timestep;
            obj.LCP_cache = SharedDataHandle(struct('t',[],'x',[],'u',[],'nargout',[], ...
                'z',[],'Mqdn',[],'wqdn',[], 'possible_contact_indices',[],'possible_limit_indices',[], ...
                'dz',[],'dMqdn',[],'dwqdn',[],'contact_data',[],'fastqp_active_set',[]));
            
            obj = setSampleTime(obj,[timestep;0]);
            
            obj.body = manip.body;
            obj = compile(obj);
        end
        
        function checkDirty(obj)
            if (obj.dirty)
                error('You''ve changed something about this model and need to manually compile it.  Use obj=compile(obj).');
            end
        end
        
        function manip = getManipulator(obj)
            manip = obj.manip;
        end
        
        function y = output(obj,t,x,u)
            checkDirty(obj);
            
            if ~isDirectFeedthrough(obj)
                u=[];
            end
            if isa(obj.getStateFrame(),'MultiCoordinateFrame')
                x_manip = double(Point(obj.getStateFrame(),x).inFrame(obj.manip.getStateFrame()));
            else
                x_manip = x;
            end
            y = obj.manip.output(t,x_manip,u);
            for i=1:length(obj.sensor)
                y = [y; obj.sensor{i}.output(obj,i+1,t,x,u)];
            end
        end
        
        function model = enableIdealizedPositionControl(model, flag)
            index = getActuatedJoints(model.manip);
            if length(unique(index))~=length(index)
                error('idealized position control currently assumes one actuator per joint');
            end
            model.position_control = logical(flag);
            model.dirty = true;
        end
        
        function model = compile(model)
            w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            model.gurobi_present = checkDependency('gurobi');
            model.manip = model.manip.compile();
            warning(w);
            
            model = setNumDiscStates(model,model.manip.getNumContStates());
            model = setNumInputs(model,model.manip.getNumInputs());
            
            if (model.position_control)
                index = getActuatedJoints(model.manip);
                pdff = pdcontrol(model,eye(model.num_u),eye(model.num_u));
                model = setInputLimits(model,model.manip.joint_limit_min(index),model.manip.joint_limit_max(index));
                model = setInputFrame(model,getInputFrame(pdff));
            else
                model = setInputLimits(model,model.manip.umin,model.manip.umax);
                model = setInputFrame(model,getInputFrame(model.manip));
            end
            model = setStateFrame(model,getStateFrame(model.manip));
            
            if ~isempty(model.sensor)
                feedthrough = model.manip.isDirectFeedthrough;
                outframe{1} = getOutputFrame(model.manip);
                stateframe{1} = getStateFrame(model.manip);
                for i=1:length(model.sensor)
                    model.sensor{i} = model.sensor{i}.compile(model,model.manip);
                    outframe{i+1} = model.sensor{i}.constructFrame(model);
                    if isa(model.sensor{i},'TimeSteppingRigidBodySensorWithState')
                        stateframe{i+1} = model.sensor{i}.constructStateFrame(model);
                    end
                    feedthrough = feedthrough || model.sensor{i}.isDirectFeedthrough;
                end
                fr = MultiCoordinateFrame.constructFrame(outframe);
                state_fr = MultiCoordinateFrame.constructFrame(stateframe);
                if ~isequal_modulo_transforms(fr,getOutputFrame(model))
                    model = setNumOutputs(model,fr.dim);
                    model = setOutputFrame(model,fr);
                end
                if ~isequal_modulo_transforms(state_fr,getStateFrame(model))
                    model = setNumDiscStates(model,state_fr.dim);
                    model = setStateFrame(model,state_fr);
                end
                model = setDirectFeedthrough(model,feedthrough);
            else
                model = setNumOutputs(model,getNumOutputs(model.manip));
                model = setOutputFrame(model,getOutputFrame(model.manip));
                model = setDirectFeedthrough(model,model.manip.isDirectFeedthrough);
            end
            model.LCP_cache.data.t = NaN;
            model.LCP_cache.data.x = NaN(model.getNumStates(),1);
            model.LCP_cache.data.u = NaN(model.getNumInputs(),1);
            model.LCP_cache.data.nargout = NaN;
            model.dirty = false;
        end
        
        function x0 = getInitialState(obj)
            if ~isempty(obj.initial_state)
                x0 = obj.initial_state;
                return;
            end
            
            x0 = obj.manip.getInitialState();
            for i=1:length(obj.sensor)
                if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
                    x0 = [x0; obj.sensor{i}.getInitialState(obj)];
                end
            end
        end
        
        function [phiC,normal,V,n,D,xA,xB,idxA,idxB] = getContactTerms(obj,q,kinsol)
            if nargin<3
                kinematics_options.compute_gradients = 1;
                kinsol = doKinematics(obj, q, [], kinematics_options);
            end
            
            [phiC,normal,d,xA,xB,idxA,idxB,mu,n,D] = obj.contactConstraints_manual(kinsol,obj.multiple_contacts);
            
            % reconstruct perturbed mu
            if ~isempty(obj.uncertain_mu)
                mu = obj.uncertain_mu*ones(length(mu),1);
            end
            % [double make sure that mu is not interweaving contactConstraints]
            
            % TODO: clean up
            nk = length(d);
            V = cell(1,2*nk);
            muI = sparse(diag(mu));
            norm_mat = sparse(diag(1./sqrt(1 + mu.^2)));
            for k=1:nk,
                V{k} = (normal + d{k}*muI)*norm_mat;
                V{nk+k} = (normal - d{k}*muI)*norm_mat;
            end
        end
        
        function [V] = getV(obj,q,kinsol)
            if nargin<3
                kinematics_options.compute_gradients = 1;
                kinsol = doKinematics(obj, q, [], kinematics_options);
            end
            
            [normal,d,mu] = obj.contactConstraints_manual_v2(kinsol);
            
            % reconstruct perturbed mu
            if ~isempty(obj.uncertain_mu)
                mu = obj.uncertain_mu*ones(length(mu),1);
            end
            % [double make sure that mu is not interweaving contactConstraints]
            
            % TODO: clean up
            nk = length(d);
            V = cell(1,2*nk);
            muI = sparse(diag(mu));
            norm_mat = sparse(diag(1./sqrt(1 + mu.^2)));
            for k=1:nk,
                V{k} = (normal + d{k}*muI)*norm_mat;
                V{nk+k} = (normal - d{k}*muI)*norm_mat;
            end
        end
        
        function [xdn,df] = solveQPOpt(obj,h,x,u)
            [xdn,df] = geval(@obj.solveQP,h,x,u,struct('grad_method','numerical'));
        end
        
        function [xdn,df] = solveQP(obj,X0)
            global x_previous
            global df_previous
            %global timestep_updated
              
            num_q = obj.manip.getNumPositions;
            num_v = obj.manip.getNumVelocities;
            num_x = num_q + num_v;
            num_u = obj.manip.getNumInputs;
            arm_dim = 8;
            obj_dim = num_q - arm_dim;
            passive_joint_dim = 1;
            
            h = X0(1);
            x = X0(2:1+num_x);
            u = X0(2+num_x:1+num_u+num_x);
            
            % try
            %     if (nargout>1)
            %         [obj,z,Mvn,wvn,dz,dMvn,dwvn] = solveLCP(obj,t,x,u);
            %     else
            %         [obj,z,Mvn,wvn] = solveLCP(obj,t,x,u);
            %     end
            %
            %     num_q = obj.manip.num_positions;
            %     q=x(1:num_q); v=x((num_q+1):end);
            %     h = obj.timestep;
            %
            %     if isempty(z)
            %         vn = wvn;
            %     else
            %         vn = Mvn*z + wvn;
            %     end
            %
            %     kinsol = obj.manip.doKinematics(q);
            %     vToqdot = obj.manip.vToqdot(kinsol);
            %     qdn = vToqdot*vn;
            %     qn = q+ h*qdn;
            %     % Find quaternion indices
            %     quat_bodies = obj.manip.body([obj.manip.body.floating] == 2);
            %     quat_positions = [quat_bodies.position_num];
            %     for i=1:size(quat_positions,2)
            %         quat_dot = qdn(quat_positions(4:7,i));
            %         if norm(quat_dot) > 0
            %             % Update quaternion by following geodesic
            %             qn(quat_positions(4:7,i)) = q(quat_positions(4:7,i)) + quat_dot/norm(quat_dot)*tan(norm(h*quat_dot));
            %             qn(quat_positions(4:7,i)) = qn(quat_positions(4:7,i))/norm(qn(quat_positions(4:7,i)));
            %         end
            %     end
            %     xdn = [qn;vn];
            %
            %     if (nargout>1)  % compute gradients
            %         if isempty(z)
            %             dqdn = dwvn;
            %         else
            %             dqdn = matGradMult(dMvn,z) + Mvn*dz + dwvn;
            %         end
            %         df = [ [zeros(num_q,1), eye(num_q), zeros(num_q,num_q+obj.num_u)]+h*dqdn; dqdn ];
            %     end
            %
            %     for i=1:length(obj.sensor)
            %         if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
            %             if (nargout>1)
            %                 [obj,xdn_sensor,df_sensor] = update(obj.sensor{i},obj,t,x,u);
            %             else
            %                 [obj,xdn_sensor] = update(obj.sensor{i},obj,t,x,u);
            %             end
            %             xdn = [xdn;xdn_sensor];
            %             if (nargout>1)
            %                 df = [df; df_sensor];
            %             end
            %         end
            %     end
            % catch
            % this function implement an update based on Todorov 2011, where
            % instead of solving the full SOCP, we make use of polyhedral
            % friction cone approximations and solve a QP.
            
            % w is a num_c * num_d disturbance vector
            % assume for now that it has the same basis vectors as contact forces
            
            % q_{k+1} = q_{k} + qd_{k+1}*h;
            % qd_{k+1} = qd_{k} + H^{-1}*(B*u-C)*h + J'*f;
            
            if obj.twoD
                num_d = 2;
            else
                num_d = 4;
            end
            dim = 3;%[3D]
            %h = timestep_updated;

            q=x(1:num_q);
            v=x(num_q+(1:num_v));
            
            kinematics_options.compute_gradients = 1;
            kinsol = doKinematics(obj, q, [], kinematics_options);
            vToqdot = obj.manip.vToqdot(kinsol);
            
            [H,C,B,dH,dC,dB] = manipulatorDynamics(obj.manip,q,v);
            
            dq = diag(max(sqrt(eps(q)), 1e-7));
            num_diff_index_set = [1,8,9,10,11,12,13,14];
            for m = 1:length(num_diff_index_set)
                index = num_diff_index_set(m);                
                [H_p,C_p,~,~,~,~] = manipulatorDynamics(obj.manip,q+dq(:,index),v);
                [H_m,C_m,~,~,~,~] = manipulatorDynamics(obj.manip,q-dq(:,index),v);
                dC(:,index) = (C_p-C_m)/(2*dq(index,index));
                dH(:,index) = reshape((H_p-H_m)/(2*dq(index,index)),[],1);
            end
            
            if obj.num_u > 0
                [phiC,normal,V,n,D,xA,xB,idxA,idxB] = getContactTerms(obj,q,kinsol);
            else
                [phiC,normal,V,n,D,xA,xB,idxA,idxB] = getContactTerms(obj,q,kinsol);
            end
            
            num_c = length(phiC);
            
            if nargin<5
                w = zeros(num_c*num_d,1);
            end
            
            %[3D]
            active = [1:length(phiC)]';%find(phiC + h*n*vToqdot*v < obj.active_threshold);
            
            phiC = phiC(active);
            normal = normal(:,active);
            
            Apts = xA(:,active);
            Bpts = xB(:,active);
            Aidx = idxA(active);
            Bidx = idxB(active);
            
            index = [1:3];% defined for 3D kuka arm
            
            [phiL,JL] = obj.manip.jointLimitConstraints(q);
            if length(phiL) ~= 2*(passive_joint_dim + arm_dim) % kuka arm
                keyboard
            end
            
            %compute for full dimension joint limit
            possible_limit_indices = true(length(phiL),1);%(phiL + h*JL*vToqdot*v) < obj.active_threshold;
            nL = sum(possible_limit_indices);
            JL = JL(possible_limit_indices,:);
            
            num_active = length(active);
            num_beta = num_active*num_d; % coefficients for friction poly
            num_full_dim = num_active*dim+nL;
            
            JA = []; JAx = []; JAy = []; JAz = [];
            dJA = []; dJAx = []; dJAy = []; dJAz = [];
            world_pts = [];
            for i=1:length(Aidx)
                [pp,J_,dJ_] = forwardKin(obj.manip,kinsol,Aidx(i),Apts(:,i));%[Ye: reshaping of dJ is commented out in forwardKin() ]
                JA = [JA; J_];
                JAx = [JAx;J_(1,:)]; JAy = [JAy;J_(2,:)]; JAz = [JAz;J_(3,:)];
                dJA = [dJA;dJ_];
                world_pts = [world_pts, pp];
            end
            
            JB = []; JBx = []; JBy = []; JBz = [];
            dJB = []; dJBx = []; dJBy = []; dJBz = [];
            for i=1:length(Bidx)
                [~,J_,dJ_] = forwardKin(obj.manip,kinsol,Bidx(i),Bpts(:,i));
                JB = [JB; J_];
                JBx = [JBx;J_(1,:)]; JBy = [JBy;J_(2,:)]; JBz = [JBz;J_(3,:)];
                dJB = [dJB;dJ_];
            end
            
            J = JA-JB;
            dJ_original = dJA-dJB;
            Jx = JAx - JBx; Jy = JAy - JBy; Jz = JAz - JBz;
            
            %Jz is
            %[p1q1z, p1q2z, ..., p1q14z
            % p2q1z, p2q2z, ..., p2q14z
            %  ...     ...  ...,  ...
            % p12q1z, p12q2z, ..., p12q14z]
            dJ = dJ_original;%dJ is reformulated later.
            %reshape(dJ(:,1),3,14,12) is
            % [p1q1x p1q2x ..., p1q14x             [p2q1x p2q2x ..., p2q14x        [p12q1x p12q2x ..., p12q14x
            %  p1q1y p1q2y ..., p1q14y       ....   p2q1y p2q2y ..., p2q14y         p12q1y p12q2y ..., p12q14y
            %  p1q1z p1q2z ..., p1q14z]             p2q1z p2q2z ..., p2q14z]        p12q1z p12q2z ..., p12q14z]
            %% compute numerical Jacobian gradient for dJ and dV
            function [J_num, V_num] = object_gradient_numerical(X0)
                h_num = X0(1);
                x_num = X0(2:1+num_q+num_v);
                u_num = X0(2+num_q+num_v:num_u+num_q+num_v);
                
                q_num=x_num(1:num_q);
                v_num=x_num(num_q+(1:obj.manip.getNumVelocities));
                
                kinematics_options.compute_gradients = 1;
                kinsol_num = doKinematics(obj, q_num, [], kinematics_options);
                
                if obj.num_u > 0
                    [~,~,V_num,~,~,xA_num,xB_num,idxA_num,idxB_num] = getContactTerms(obj,q_num,kinsol_num);
                else
                    [~,~,V_num,~,~,xA_num,xB_num,idxA_num,idxB_num] = getContactTerms(obj,q_num,kinsol_num);
                end
                
                V_num = horzcat(V_num{:});
                I = eye(num_c*num_d);
                V_cell_num = cell(1,num_active);
                for ii=1:num_c+nL
                    if ii<=num_active
                        % is a contact point
                        idx_beta_num = active(ii):num_c:num_c*num_d;
                        V_cell_num{ii} = V_num*I(idx_beta_num,:)'; % basis vectors for ith contact
                    end
                end
                V_num = blkdiag(V_cell_num{:},eye(nL));
                
                Apts_num = xA_num(:,active);
                Bpts_num = xB_num(:,active);
                Aidx_num = idxA_num(active);
                Bidx_num = idxB_num(active);
                
                JA_num = [];
                world_pts = [];
                for k=1:length(Aidx_num)
                    [~,J_num_,~] = forwardKin(obj.manip,kinsol_num,Aidx_num(k),Apts_num(:,k));%[Ye: reshaping of dJ is commented out in forwardKin() ]
                    JA_num = [JA_num; reshape(J_num_(index,:),[],1)];
                end
                
                JB_num = [];
                for k=1:length(Bidx_num)
                    [~,J_num_,~] = forwardKin(obj.manip,kinsol_num,Bidx_num(k),Bpts_num(:,k));
                    JB_num = [JB_num; reshape(J_num_(index,:),[],1)];
                end
                
                J_num = JA_num-JB_num;
            end
            
            function [V_num] = V_gradient_numerical(X0)
                h_num = X0(1);
                x_num = X0(2:1+num_q+num_v);
                u_num = X0(2+num_q+num_v:num_u+num_q+num_v);
                
                q_num=x_num(1:num_q);
                v_num=x_num(num_q+(1:obj.manip.getNumVelocities));
                
                kinematics_options.compute_gradients = 1;
                kinsol_num = doKinematics(obj, q_num, [], kinematics_options);
                     
                if obj.num_u > 0
                    [V_num] = getV(obj,q_num,kinsol_num);
                else
                    [V_num] = getV(obj,q_num,kinsol_num);
                end
                
                V_num = horzcat(V_num{:});
                I = eye(num_c*num_d);
                V_cell_num = cell(1,num_active);
                for ii=1:num_c+nL
                    if ii<=num_active
                        % is a contact point
                        idx_beta_num = active(ii):num_c:num_c*num_d;
                        V_cell_num{ii} = V_num*I(idx_beta_num,:)'; % basis vectors for ith contact
                    end
                end
                V_num = blkdiag(V_cell_num{:},eye(nL));
            end
            
            function [JA_ground_num, JB_num] = AB_gradient_numerical(X0)
                x_num = X0(2:1+num_q+num_v);                
                q_num=x_num(1:num_q);
                
                kinematics_options.compute_gradients = 1;
                kinsol_num = doKinematics(obj, q_num, [], kinematics_options);
                
                xB_num = obj.computexB(kinsol_num);
                
                for k=1:groundcontact_num
                    JB_ground_obj_ori = dJB_finger_obj_ori_analytical(q_num,xB_num(1,k),xB_num(2,k),xB_num(3,k));
                    J_B_num_new = [zeros(3,8),eye(3),JB_ground_obj_ori];
                    J_A_num_new = zeros(3,length(J_B_num_new(3,:)));%J_B_num_new;%
                    J_A_num_new(3,:) = zeros(1,length(J_A_num_new(3,:)));% set it to zero, since point A is always on the ground
                    JB_num(:,k) = reshape(J_B_num_new,[],1);
                    %as a by-product, we can derive JA_ground_num
                    JA_ground_num(:,k) = reshape(J_A_num_new,[],1);
                end
                
                for k=1:fingercontact_num
                    JB_finger_obj_ori = dJB_finger_obj_ori_analytical(q_num,xB_num(1,groundcontact_num+k),xB_num(2,groundcontact_num+k),xB_num(3,groundcontact_num+k));
                    J_B_num_new = [zeros(3,8),eye(3),JB_finger_obj_ori];
                    JB_num(:,k+groundcontact_num) = reshape(J_B_num_new,[],1);
                end
            end
            
            function [JA_ground_num, JB_num] = AB_gradient_numerical_new(X0)
                x_num = X0(2:1+num_q+num_v);                
                q_num=x_num(1:num_q);
                
                kinematics_options.compute_gradients = 1;
                kinsol_num = doKinematics(obj, q_num, [], kinematics_options);
                
                xB_num = obj.computexB(kinsol_num);
                
                for k=1:groundcontact_num
                    %JB_ground_obj_ori = dJB_finger_obj_ori_analytical(q_num,xB_num(1,k),xB_num(2,k),xB_num(3,k));
                    %J_B_num_new = [zeros(3,8),eye(3),JB_ground_obj_ori];
                    [~,J_B_num_new,~] = forwardKin(obj.manip,kinsol_num,Bidx(k),xB_num(:,k));
                    J_A_num_new = zeros(3,length(J_B_num_new(3,:)));%J_B_num_new;%
                    J_A_num_new(3,:) = zeros(1,length(J_A_num_new(3,:)));% set it to zero, since point A is always on the ground
                    JB_num(:,k) = reshape(J_B_num_new,[],1);
                    %as a by-product, we can derive JA_ground_num
                    JA_ground_num(:,k) = reshape(J_A_num_new,[],1);
                end
                
                for k=1:fingercontact_num
                    %JB_finger_obj_ori = dJB_finger_obj_ori_analytical(q_num,xB_num(1,groundcontact_num+k),xB_num(2,groundcontact_num+k),xB_num(3,groundcontact_num+k));
                    %J_B_num_new = [zeros(3,8),eye(3),JB_finger_obj_ori];
                    [~,J_B_num_new,~] = forwardKin(obj.manip,kinsol_num,Bidx(groundcontact_num+k),xB_num(:,groundcontact_num+k));
                    JB_num(:,k+groundcontact_num) = reshape(J_B_num_new,[],1);
                end
            end
            
            fcn = @object_gradient_numerical;
            %fcnAB = @AB_gradient_numerical;
            fcnAB = @AB_gradient_numerical_new;
            fcnV = @V_gradient_numerical;
            
            %tic
            dV = zeros(num_full_dim,(num_d*num_c+nL),num_q);
            % for i=1:num_q
            %     rr = max(sqrt(eps(X0)),1e-7);
            %     rr(1:i) = 0;%corresponding to time step h
            %     m = i+2;% only q component affects Jacobian J
            %     rr(m:end) = rr(m:end) - rr(m:end);% set all other elements (unrelated to q state) to be zero
            %     X0_p = X0 + rr;
            %     X0_m = X0 - rr;
            %
            %     % % old approach for comparison
            %     % [J_object_p,V_num_p] = feval(fcn,X0_p);
            %     % [J_object_m,V_num_m] = feval(fcn,X0_m);
            %     %
            %     % dJ_ori(:,i) = (J_object_p - J_object_m)/rr(i+1);
            %     % dV_ori(:,:,i) = (V_num_p - V_num_m)/rr(i+1);
            %
            %     %% new approach, 4 times faster
            %
            %     [V_num_p] = feval(fcnV,X0_p);
            %     [V_num_m] = feval(fcnV,X0_m);
            %     dV(:,:,i) = (V_num_p - V_num_m)/(2*rr(i+1));
            %
            %     q = X0(2:15);
            %     arm_dim = 8;
            %     obj_dim = 6;
            %     groundcontact_num = 4;
            %     fingercontact_num = 8;
            %
            %     %for this part, we only compute 4 ground contact points for
            %     %JA, and full contact points (4 ground contact points + 8
            %     %finger contact points) for JB
            %     [JA_ground_num_p, JB_ground_finger_num_p] = feval(fcnAB,X0_p);
            %     [JA_ground_num_m, JB_ground_finger_num_m] = feval(fcnAB,X0_m);
            %
            %     kinematics_options.compute_gradients = 1;
            %     kinsol_p = doKinematics(obj, X0_p(2:15), [], kinematics_options);
            %     kinsol_m = doKinematics(obj, X0_m(2:15), [], kinematics_options);
            %
            %     [~,~,~,~,~,xA_p,xB_p,~,~] = getContactTerms(obj,X0_p(2:15),kinsol_p);
            %     [~,~,~,~,~,xA_m,xB_m,~,~] = getContactTerms(obj,X0_m(2:15),kinsol_m);
            %
            %     if i <= arm_dim
            %         dJA_ground = zeros(groundcontact_num*3*length(q),1);
            %         dJA_finger = [];
            %         for fingercontactIndx = 1:fingercontact_num
            %             dJA_single_finger = [dJA_analytical(q,fingercontactIndx,i);zeros(obj_dim*3,1)];
            %
            %             kinematics_options.compute_gradients = 1;
            %             kinsol_p = doKinematics(obj, X0_p(2:15), [], kinematics_options);
            %             [~,JA_single_finger_p,~] = forwardKin(obj.manip,kinsol_p,Aidx(groundcontact_num+fingercontactIndx),xA_p(:,groundcontact_num+fingercontactIndx));
            %
            %             kinematics_options.compute_gradients = 1;
            %             kinsol_m = doKinematics(obj, X0_m(2:15), [], kinematics_options);
            %             [~,JA_single_finger_m,~] = forwardKin(obj.manip,kinsol_m,Aidx(groundcontact_num+fingercontactIndx),xA_m(:,groundcontact_num+fingercontactIndx));
            %             dJA_single_finger_new = (JA_single_finger_p - JA_single_finger_m)/(2*rr(i+1));
            %             dJA_single_finger_new = [reshape(dJA_single_finger_new,[],1)];
            %
            %             %[~,JA_single_finger,dJA_single_finger] = forwardKin(obj.manip,kinsol,Aidx(groundcontact_num+fingercontactIndx),xA(:,groundcontact_num+fingercontactIndx));
            %             %dJA_single_finger_new = [reshape(JA_single_finger,[],1)];
            %
            %             dJA_finger = [dJA_finger;dJA_single_finger_new];
            %         end
            %         dJA(:,i) = [dJA_ground;dJA_finger];
            %     elseif i <= arm_dim+obj_dim/2
            %         dJA(:,i) = zeros((groundcontact_num+fingercontact_num)*3*length(q),1);
            %     else
            %         %dJA_ground only becomes non-zero for 3-DOF object orientation pertubation
            %         dJA_ground = (JA_ground_num_p - JA_ground_num_m)/(2*rr(i+1));
            %         dJA_ground = reshape(dJA_ground,[],1);
            %         dJA(:,i) = [dJA_ground;zeros(fingercontact_num*3*length(q),1)];
            %     end
            %
            %     dJB_ground_finger = (JB_ground_finger_num_p - JB_ground_finger_num_m)/(2*rr(i+1));
            %     dJB(:,i) = reshape(dJB_ground_finger,[],1);
            %
            %     dJ_old(:,i) = dJA(:,i) - dJB(:,i);
            %     %run permute(reshape(dJ(:,1),3,14,12),[2,1,3]) to see
            %     % [p1q1x p1q1y p1q1z             [p2q1x p2q1y p2q1z        [p12q1x p12q1y p12q1z
            %     %  p1q2x p1q2y p1q2z       ....   p2q2x p2q2y p2q2z         p12q2x p12q2y p12q2z
            %     %  ...   ...   ...                ...   ...   ...           ...   ...   ...
            %     %  p1q14x p1q14y p1q14z]          p2q14x p2q14y p2q14z]    p12q14x p12q14y p12q14z]
            % end

            groundcontact_num = 4;
            fingercontact_num = 8;
            for i=1:num_q
                rr = max(sqrt(eps(X0)),1e-7);
                rr(1:i) = 0;%corresponding to time step h
                m = i+2;% only q component affects Jacobian J
                rr(m:end) = rr(m:end) - rr(m:end);% set all other elements (unrelated to q state) to be zero
                X0_p = X0 + rr;
                X0_m = X0 - rr;
                
                % % old approach for comparison
                % [J_object_p,V_num_p] = feval(fcn,X0_p);
                % [J_object_m,V_num_m] = feval(fcn,X0_m);
                %
                % dJ_ori(:,i) = (J_object_p - J_object_m)/rr(i+1);
                % dV_ori(:,:,i) = (V_num_p - V_num_m)/rr(i+1);
                
                %% new approach, 4 times faster
                
                [V_num_p] = feval(fcnV,X0_p);
                [V_num_m] = feval(fcnV,X0_m);
                dV(:,:,i) = (V_num_p - V_num_m)/(2*rr(i+1));
                
                q = X0(2:1+num_q);
                
                %for this part, we only compute 4 ground contact points for
                %JA, and full contact points (4 ground contact points + 8
                %finger contact points) for JB
                %[JA_ground_num_p, JB_ground_finger_num_p] = feval(fcnAB,X0_p);
                %[JA_ground_num_m, JB_ground_finger_num_m] = feval(fcnAB,X0_m);

                kinematics_options.compute_gradients = 1;
                kinsol_p = doKinematics(obj, X0_p(2:1+num_q), [], kinematics_options);
                kinsol_m = doKinematics(obj, X0_m(2:1+num_q), [], kinematics_options);
                
                [~,~,~,~,~,xA_p,xB_p,~,~] = getContactTerms(obj,X0_p(2:1+num_q),kinsol_p);
                [~,~,~,~,~,xA_m,xB_m,~,~] = getContactTerms(obj,X0_m(2:1+num_q),kinsol_m);
                
                J_B_p_num = []; J_B_m_num = []; J_A_p_num = []; J_A_m_num = [];
                
                kinsol_p = doKinematics(obj, X0_p(2:1+num_q), [], kinematics_options);
                for k=1:groundcontact_num+fingercontact_num
                    [~,J_B_p_num_new,~] = forwardKin(obj.manip,kinsol_p,Bidx(k),xB_p(:,k));
                    J_B_p_num = [J_B_p_num;J_B_p_num_new];
                end
                
                kinsol_m = doKinematics(obj, X0_m(2:1+num_q), [], kinematics_options);
                for k=1:groundcontact_num+fingercontact_num
                    [~,J_B_m_num_new,~] = forwardKin(obj.manip,kinsol_m,Bidx(k),xB_m(:,k));
                    J_B_m_num = [J_B_m_num;J_B_m_num_new];
                end
                
                kinsol_p = doKinematics(obj, X0_p(2:1+num_q), [], kinematics_options);
                for k=1:groundcontact_num+fingercontact_num
                    [~,J_A_p_num_new,~] = forwardKin(obj.manip,kinsol_p,Aidx(k),xA_p(:,k));
                    J_A_p_num = [J_A_p_num;J_A_p_num_new];
                end
                
                kinsol_m = doKinematics(obj, X0_m(2:1+num_q), [], kinematics_options);
                for k=1:groundcontact_num+fingercontact_num
                    [~,J_A_m_num_new,~] = forwardKin(obj.manip,kinsol_m,Aidx(k),xA_m(:,k));
                    J_A_m_num = [J_A_m_num;J_A_m_num_new];
                end
                
                dJA_num = (J_A_p_num - J_A_m_num)/(2*rr(i+1));
                dJB_num = (J_B_p_num - J_B_m_num)/(2*rr(i+1));
                
                dJ_original2 = dJA_num - dJB_num; 
                
                %dJ(:,i) = dJA(:,i) - dJB(:,i); 
                dJ(:,i) = reshape(permute(reshape(dJ_original2,3,12,num_q),[1,3,2]),[],1);
                %run permute(reshape(dJ(:,1),3,14,12),[2,1,3]) to see
                % [p1q1x p1q1y p1q1z             [p2q1x p2q1y p2q1z        [p12q1x p12q1y p12q1z
                %  p1q2x p1q2y p1q2z       ....   p2q2x p2q2y p2q2z         p12q2x p12q2y p12q2z
                %  ...   ...   ...                ...   ...   ...           ...   ...   ...
                %  p1q14x p1q14y p1q14z]          p2q14x p2q14y p2q14z]    p12q14x p12q14y p12q14z]
                
            end
            
            %% end of numerical Jacobian gradient for dJ and dV
            
            % dJ new sequence, stack Aidx first and then num_q
            dJx = []; dJy = []; dJz = [];
            for j=1:length(Aidx)
                for i=1:num_q
                    dJx = [dJx;dJ(dim*(i-1)+dim*num_q*(j-1)+1,:)];
                    dJy = [dJy;dJ(dim*(i-1)+dim*num_q*(j-1)+2,:)];%[3D]
                    dJz = [dJz;dJ(dim*(i-1)+dim*num_q*(j-1)+3,:)];
                end
            end
            dJD{1} = dJx;
            dJD{2} = dJy;%[3D]
            
            % dJ new sequence, stack Aidx first and then num_q, only used later for dMvn derivation
            dJx_new = []; dJy_new = []; dJz_new = [];
            for i=1:num_q
                for j=1:length(Aidx)
                    dJx_new = [dJx_new;dJ(dim*(i-1)+dim*num_q*(j-1)+1,:)];
                    dJy_new = [dJy_new;dJ(dim*(i-1)+dim*num_q*(j-1)+2,:)];%[3D]
                    dJz_new = [dJz_new;dJ(dim*(i-1)+dim*num_q*(j-1)+3,:)];
                end
            end
            dJD_new{1} = dJx_new;
            dJD_new{2} = dJy_new;%[3D]
            
            %joint limit
            J_contact = J;
            J = [J;JL];
            phiL = phiL(possible_limit_indices);
            phi = [phiC;phiL];
            
            global external_force_index;
             
            addExternalForce = 0;
            if addExternalForce == 1
                w_force = load('external_force.dat');
                tau_external = [zeros(8+1,1);w_force(external_force_index);zeros(1+3,1)];
            else
                tau_external = zeros(num_q,1);
            end
            
            if (obj.num_u>0)
                tau = B*u - C + tau_external;
                dtau = [zeros(num_v,1), matGradMult(dB,u) - dC, B];
            else
                tau = -C;
                dtau = [zeros(num_v,1), -dC, zeros(size(B))];
            end
            Hinv = inv(H);
            
            V = horzcat(V{:});
            I = eye(num_c*num_d);
            V_cell = cell(1,num_active);
            v_min = zeros(length(phi),1);
            for i=1:length(phi)
                if i<=num_active
                    % is a contact point
                    idx_beta = active(i):num_c:num_c*num_d;
                    V_cell{i} = V*I(idx_beta,:)'; % basis vectors for ith contact
                end
                v_min(i) = -phi(i)/h;
            end
            V = blkdiag(V_cell{:},eye(nL));
            
            A = J*vToqdot*Hinv*vToqdot'*J';
            c = J*vToqdot*v + J*vToqdot*Hinv*tau*h;
            
            % contact smoothing matrix
            R_min = 1e-3;
            R_max = 1e-1;
            r = zeros(num_active,1);
            r(phiC>=obj.phi_max) = R_max;
            r(phiC<=obj.contact_threshold) = R_min;
            ind = (phiC > obj.contact_threshold) & (phiC < obj.phi_max);
            
            %nonlinear function
            %y = (phiC(ind)-obj.contact_threshold)./(obj.phi_max - obj.contact_threshold)*2 - 1; % scale between -1,1
            %r(ind) = R_min + R_max./(1+exp(-10*y));
            
            % linear function
            r(ind) = R_min + (R_max - R_min)*(phiC(ind)-obj.contact_threshold)./(obj.phi_max - obj.contact_threshold);
            
            % new continuous nonlinear function
            %ind = ones(length(phiC),1);
            %y = (phiC(ind)-obj.contact_threshold)./(obj.phi_max - obj.contact_threshold)*2 - 1;
            %r(ind) = R_min + R_max./(1+exp(-20*y));
            
            r = repmat(r,1,dim)';
            %R = diag([r(:)',r(:)']);
            R = diag(r(:));
            
            for i=1:num_q
                for j=1:length(ind)
                    %linear function
                    if ind(j) == 0
                      drdq(j,i) = 0;
                    else
                      drdq(j,i) = (R_max - R_min)*n(j,i)./(obj.phi_max - obj.contact_threshold);
                    end
                    %nonlinear scaling
                    %drdq(j,i) = (R_max*10*exp(-20*y(j))./(1+exp(-20*y(j)))^2)*(2*n(j,i)/(obj.phi_max - obj.contact_threshold));
                end
                dRdq_tmp = repmat(drdq(:,i),1,dim)';
                dRsubdq(:,:,i) = diag(reshape(dRdq_tmp,[],1));
            end
            
            if any(ind > 0)
                S_weighting_unit = diag(ones(3,1));%[3D]
            else
                S_weighting_unit = diag(ones(3,1));%[3D]
                % right now, use a unified weighting coeff for x
                % direction regardless whether the brick enters the
                % safety region. But it could be tuned to different
                % values such that different weighting can be used
            end
            
            % joint limit smoothing matrix
            W_min = 1e-3;
            W_max = 1e-1;
            w = zeros(nL,1);
            w(phiL>=obj.phiL_max) = W_max; 
            w(phiL<=obj.contact_threshold) = W_min;
            ind = (phiL > obj.contact_threshold) & (phiL < obj.phiL_max);
            
            % nonlinear function
            %y = (phiL(ind)-obj.contact_threshold)./(obj.phiL_max - obj.contact_threshold)*2 - 1; % scale between -1,1
            %w(ind) = W_min + W_max./(1+exp(-10*y));
            
            % linear function
            w(ind) = W_min + (W_max - W_min)*(phiL(ind)-obj.contact_threshold)./(obj.phiL_max - obj.contact_threshold);
            
            %new continuous function
            %ind = ones(length(phiL),1);
            %y = (phiL(ind)-obj.contact_threshold)./(obj.phiL_max - obj.contact_threshold)*2 - 1;
            %w(ind) = W_min + W_max./(1+exp(-20*y));
            
            W = diag(w(:));
            
            for i=1:num_q
                for j=1:length(ind)
                    %linear function
                    if ind(j) == 0
                       dwdq(j,i) = 0;
                    else
                       dwdq(j,i) = (W_max - W_min)*JL(j,i)./(obj.phiL_max - obj.contact_threshold);
                    end
                    %nonlinear function
                    %dwdq(j,i) = (W_max*10*exp(-20*y(j))./(1+exp(-20*y(j)))^2)*(2*JL(j,i)/(obj.phiL_max - obj.contact_threshold));
                end
                dWdq(:,:,i) = diag(dwdq(:,i));
            end
            
            R = blkdiag(R,W);%1e-1*eye(52);%
            
            for i=1:num_q
                dRdq(:,:,i) = blkdiag(dRsubdq(:,:,i),dWdq(:,:,i));%zeros(52);%
            end
            
            num_params = num_beta+nL;
            % lambda_ub = zeros(num_params,1);
            % scale_fact = 1e3;
            % phiC_pos = phiC;
            % phiC_pos(phiC<0)=0;
            % lambda_ub(1:num_beta) = repmat(max(0.01, scale_fact*(obj.phi_max./phiC_pos - 1.0)),1,num_d)';
            % phiL_pos = phiL;
            % phiL_pos(phiL<0)=0;
            % lambda_ub(num_beta+(1:nL)) = max(0.01, scale_fact*(obj.phi_max./phiL_pos - 1.0));
            lambda_ub = inf*ones(num_params,1);% disable this unnecessary bounding constraint
            
            % define weighting parameters
            for i=1:num_active
                S_weighting_array{i} = S_weighting_unit;
            end
            for i=1:nL
                S_weighting_array{num_active+i} = 1;
            end
            S_weighting = blkdiag(S_weighting_array{:});
            
            Q = 0.5*V'*S_weighting*(A+R)*S_weighting*V + 1e-6*eye(num_params);
             
            % N*(A*z + c) - v_min \ge 0
            Ain = zeros(num_active+nL,num_params);
            bin = zeros(num_active+nL,1);
            for i=1:num_active
                idx = (i-1)*dim + (1:dim);
                Ain(i,:) = normal(:,i)'*A(idx,:)*V;
                bin(i) = v_min(i) - normal(:,i)'*c(idx);
            end
            for i=1:nL
                idx = num_active*dim + i;
                Ain(i+num_active,:) = A(idx,:)*V;
                bin(i+num_active) = v_min(i+num_active) - c(idx);
            end
            
            Ain_fqp = full([-Ain; -eye(num_params); eye(num_params)]);
            bin_fqp = [-bin; zeros(num_params,1); lambda_ub];
            
            %[result_qp,info_fqp] = fastQPmex({Q},V'*c,Ain_fqp,bin_fqp,[],[],obj.LCP_cache.data.fastqp_active_set);
            
            %if 1 % info_fqp<0
            %disp('calling gurobi');
            model.LCP_cache.data.fastqp_active_set = [];
            gurobi_options.outputflag = 0; % verbose flag
            gurobi_options.method = 1; % -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier
            
            try
                model.Q = sparse(Q);
                model.obj = V'*S_weighting*c;
                model.A = sparse(Ain);
                model.rhs = bin;
                model.sense = repmat('>',length(bin),1);
                model.lb = zeros(num_params,1);
                model.ub = lambda_ub;
                result = gurobi(model,gurobi_options);
                result_qp = result.x;
            catch
                display('gurobi solve failure');
                %keyboard
                xdn = x_previous;
                df = df_previous;
                return
            end
            f = S_weighting*V*(result_qp);% each 3x1 block is for one contact point, x, y, and z direction are all negative values, since it points from B to A.
            %f is in the sequence of p1q1x, p1q2x,...,p1q14x,p1q1y,p1q2y,...,p1q14y,p1q1z,p1q2z,...,p1q14z,
            %p2q1x, p2q2x,...,p2q14x,p2q1y,p2q2y,...,p2q14y,p2q1z,p2q2z,...,p2q14z,...  
            active_set = find(abs(Ain_fqp*result_qp - bin_fqp)<1e-6);
            obj.LCP_cache.data.fastqp_active_set = active_set;
            
            % mu_tolerance = 1e-5;
            % for i=1:num_active
            %     if abs(sqrt(f(1+(i-1)*3)^2+f(2+(i-1)*3)^2)/f(i*3)) > mu(1) + mu_tolerance% need to be modified for non-terrain based surface
            %         keyboard
            %     end
            % end
            %f_vec = [f_vec,f];
            
            %% checker that the analytical solution from KKT condition
            % gives correct contact force solution
            Ain_fqp_active = Ain_fqp(active_set,:);
            bin_fqp_active = bin_fqp(active_set);
            
            Qinv = pinv(2*Q);
            % this analytical solution is slightly different than the
            % one in SJ Wright's numerical optimizaiton
            G = Qinv - Qinv*Ain_fqp_active' ...
                *pinv(Ain_fqp_active*Qinv*Ain_fqp_active')*Ain_fqp_active*Qinv;
            E = Qinv*Ain_fqp_active'*pinv(Ain_fqp_active*Qinv*Ain_fqp_active');
            E = -E;
            % note that, Ain in Ain_fqp_active is negative, thus the sign of E is negative since there are three
            % Ain_fqp_active multiplied in E.
            F = pinv(Ain_fqp_active*Qinv*Ain_fqp_active');
            
            result_analytical = - G*V'*S_weighting*c - E*bin_fqp_active;
            %f = S_weighting*V*(result_analytical);
            
            % check contact force solved by analytical solution and
            % gurobi solver, they should be the same
            lambda_diff = result_analytical - result_qp;
            lambda_diff_sum = sum(abs(lambda_diff));
            if lambda_diff_sum > 1e-2
                disp('Error: contact force solved by analytical solution and gurobi solver are different')
            end
            %debugging
            %% end of checker
            
            % update state at next time step
            vn = v + Hinv*(tau*h + vToqdot'*J'*f);
            qdn = vToqdot*vn;
            qn = q + qdn*h;
            
            wvn = v + h*Hinv*tau;
            Mvn = Hinv*vToqdot'*J';% note that the sign of J is negative
            
            %% compute gradient component
            total_possible_contact_point = num_active;%[Ye: to be tuned for other systems]
            possible_contact_indices = zeros(total_possible_contact_point,1);
            for i = 1:length(active)
                possible_contact_indices(active(i)) = true(1);
            end
            
            nP = 0;
            nC = num_active;
            mC = 2;%length(D);% 2D, normally mC = 2 for 3D;[3D]
            
            % derive jacobian
            J_size = [nP + (mC+1)*nC + nL,num_q];
            
            %lb = zeros(nL+nP+(mC+2)*nC,1);
            %ub = Big*ones(nL+nP+(mC+2)*nC,1);
            JD{1} = Jx;
            JD{2} = Jy;%[3D]
            JD = vertcat(JD{:});
            
            obj.LCP_cache.data.possible_contact_indices=possible_contact_indices;
            
            % new dJ to accommodate lcp formulation
            dJ_qp = zeros(prod(J_size), num_q); % was sparse, but reshape trick for the transpose below didn't work
            possible_contact_indices_found = find(possible_contact_indices);
            n_size = [numel(possible_contact_indices), num_q];
            col_indices = 1 : num_q;
            dJz = getSubMatrixGradient(reshape(dJz, [], num_q), possible_contact_indices_found, col_indices, n_size);%change nothing
            dJ_qp = setSubMatrixGradient(dJ_qp, dJz_new, [3:3:(mC+1)*nC], col_indices, J_size);%[1 : J_size(2)]
            JD_size = size(JD);
            dJD_matrix = zeros(prod(JD_size), num_q);
            
            JD_single_size = size(Jx);
            dJD_single_matrix{1} = zeros(prod(JD_single_size), num_q);%is same as dJx_new
            dJD_single_matrix{2} = zeros(prod(JD_single_size), num_q);%is same as dJy_new
            row_start = 0;
            for i = 1 : mC
                dJD_possible_contact = getSubMatrixGradient(dJD_new{i}, possible_contact_indices_found, col_indices, n_size);
                dJD_matrix = setSubMatrixGradient(dJD_matrix, dJD_possible_contact, row_start + (1 : nC), col_indices, JD_size);
                dJD_single_matrix{i} = setSubMatrixGradient(dJD_single_matrix{i}, dJD_possible_contact, (1 : nC), col_indices, JD_single_size);
                row_start = row_start + nC;
            end
            %dJD_new = dJD_matrix;
            
            dJ_qp = setSubMatrixGradient(dJ_qp, dJD_single_matrix{1}, [1:3:(mC+1)*nC], col_indices, J_size);
            dJ_qp = setSubMatrixGradient(dJ_qp, dJD_single_matrix{2}, [2:3:(mC+1)*nC], col_indices, J_size);
            %reshape(dJ_qp(:,1),52,14) is
            % [p1q1x p1q2x ..., p1q14x
            %  p1q1y p1q2y ..., p1q14y      
            %  p1q1z p1q2z ..., p1q14z
            %  p2q1x p2q2x ..., p2q14x
            %  p2q1y p2q2y ..., p2q14y
            %  p2q1z p2q2z ..., p2q14z
            %  ...   ...   ..., ...
            %  p12q1x p12q2x ..., p12q14x
            %  p12q1y p12q2y ..., p12q14y
            %  p12q1z p12q2z ..., p12q14z
            %  zeros(16,14) for joint limit]
            dJ_qp = reshape(dJ_qp, [], num_q^2); % expand dJ in column format. For one column, stack Aidx first and then num_q
            
            %wvn = v + Hinv*tau;
            %here I add the first column Hinv*tau, which is different
            %from the LCP formulation. h is not a decision variable. In
            %fact, it does not matter. But when I compare numerical value and analytical values, the first state h is
            %perturbed as well. Thus, it is better to add Hinv*tau here.
            dwvn = [Hinv*tau, zeros(num_v,num_q),eye(num_v),zeros(num_v,obj.num_u)] + ...
                h*Hinv*dtau - [zeros(num_v,1),h*Hinv*matGradMult(dH(:,1:num_q),Hinv*tau),zeros(num_v,num_q),zeros(num_v,obj.num_u)];
            dJtranspose_qp = reshape(permute(reshape(dJ_qp,size(J,1),size(J,2),[]),[2,1,3]),numel(J),[]);% all it does is to get dJ transpose
            % double verify that this dJtranspose_qp is correct.
            %reshape(dJtranspose_qp(:,1),14,52) is
            % [p1q1x p1q1y p1q1z, p2q1x, p2q1y, p2q1z ..., p12q1x, p12q1y, p12q1z         zeros(1,16)
            %  p1q2x p1q2y p1q2z, p2q2x, p2q2y, p2q2z ..., p12q2x, p12q2y, p12q2z         zeros(1,16)
            %  ...   ...   ...     ...    ...    ...   ..., ...     ...      ...          zeros(1,16)
            %  p1q14x p1q14y p1q14z, p2q14x, p2q14y, p2q14z ..., p12q14x, p12q14y, p12q14z  zeros(1,16)]% the ending zeros are for joit limit
            
            % but we don't need it here because our original dJ is already stacked in [x1(q1);y1(q1);z1(q1);x2(q2);y2(q2);z2(q2);x3(q2);y3(q2);z3(q2);...]
            % the original dJ in LCP is stacked in [x1(q1);x2(q2);x3(q2);y1(q1);y2(q2);y3(q2);z1(q1);z2(q2);z3(q2);...]
            dMvn = [zeros(numel(Mvn),1),reshape(Hinv*reshape(dJtranspose_qp - matGradMult(dH(:,1:num_q),Hinv*J'),num_q,[]),numel(Mvn),[]),zeros(numel(Mvn),num_v+obj.num_u)];
            % Hinv*reshape(dJtranspose,num_q,[]) is the same as Hinv*vToqdot'*dJdq(:,:,i)' with i=1:num_q
            % Hinv*reshape(matGradMult(dH(:,1:num_q),Hinv*J'),num_q,[]) is same as Hinv*dHdq(:,:,i)*Hinv*vToqdot'*J' with i=1:num_q
            
            % compute dlambda/dx, dlambda/du
            % this part is a key difference from LCP solver since the
            % force vector is reformulated.
            lambda = f;
            
            b = V'*S_weighting*c;
            %A = J*vToqdot*Hinv*vToqdot'*J';
            %c = J*vToqdot*v + J*vToqdot*Hinv*tau*h;
                        
            %% partial derivative of A, b w.r.t. h, q, v and u
            dcdh = J*vToqdot*Hinv*tau;
            dbdh = V'*S_weighting*dcdh;
            % preallocate size
            dHdq = zeros(num_q,num_q,num_q);
            dCdq = zeros(num_q,1,num_q);
            dJdq = zeros(num_full_dim,num_q,num_q);
            dAdq_tmp = zeros(num_full_dim,num_full_dim,num_q);
            dAdq = zeros(num_params,num_params,num_q);
            dAdv = zeros(num_params,num_params,num_q);
            dbdq = zeros(num_params,num_q);
            dcdq = zeros(num_full_dim,num_q);
            for i=1:num_q % assume num_q == num_v
                dHdq(:,:,i) = reshape(dH(:,i),num_q,num_q);
                dCdq(:,:,i) = reshape(dC(:,i),num_q,1);
                dJx_reshape = reshape(dJx(:,i),num_q,[])';
                dJy_reshape = reshape(dJy(:,i),num_q,[])';%[3D]
                dJz_reshape = reshape(dJz(:,i),num_q,[])';
                dJdq_tmp = [];
                for j=1:nC
                    dJdq_tmp = [dJdq_tmp;dJx_reshape(j,:);dJy_reshape(j,:);dJz_reshape(j,:)];%3D
                end
                dJdq(:,:,i) = [dJdq_tmp;zeros(nL,num_q)];%[3D]
                dAdq_tmp(:,:,i) = -J*vToqdot*Hinv*dHdq(:,:,i)*Hinv*vToqdot'*J' + dJdq(:,:,i)*vToqdot*Hinv*vToqdot'*J' ...
                    +J*vToqdot*Hinv*vToqdot'*dJdq(:,:,i)';
                
                dAdq(:,:,i) = V'*S_weighting*(dAdq_tmp(:,:,i)+dRdq(:,:,i))*S_weighting*V + dV(:,:,i)'*S_weighting*(A+R)*S_weighting*V ...
                    + V'*S_weighting*(A+R)*S_weighting*dV(:,:,i);
                
                dcdq(:,i) = dJdq(:,:,i)*vToqdot*(v+Hinv*tau*h) - J*vToqdot*Hinv*dHdq(:,:,i)*Hinv*tau*h ...
                    - J*vToqdot*Hinv*dCdq(:,:,i)*h;
                dbdq(:,i) = V'*S_weighting*dcdq(:,i) + dV(:,:,i)'*S_weighting*c;
                dAdv(:,:,i) = zeros(num_params);
            end
            dCdv = dC(:,num_q+1:2*num_q);
            dcdv = J*vToqdot - J*vToqdot*Hinv*dCdv*h;
            dbdv = V'*S_weighting*dcdv;
            
            dAdu = zeros(num_params,num_params,obj.num_u);
            dcdu = zeros(num_full_dim,obj.num_u);
            dbdu = zeros(num_params,obj.num_u);
            for i=1:obj.num_u
                dAdu(:,:,i) = zeros(num_params);
                dcdu(:,i) = J*vToqdot*Hinv*B(:,i)*h;
                dbdu(:,i) = V'*S_weighting*dcdu(:,i);
            end
            
            %% partial derivative of equality constraints (active)
            num_dynamicsConstraint = num_active+nL;
            dynamicsConstraint_active_set = find(active_set<=num_active);%[Ye: to be tuned for different systems]
            jointConstraint_active_set = find(active_set>num_active & active_set<=num_dynamicsConstraint);%[Ye: to be tuned for different systems]
            boundingConstraint_active_set = find(active_set>num_dynamicsConstraint);
            num_dynamicsConstraint_active_set = length(dynamicsConstraint_active_set);
            num_jointConstraint_active_set = length(jointConstraint_active_set);
            num_boundingConstraint_active_set = length(boundingConstraint_active_set);
            num_constraint_active_set = num_dynamicsConstraint_active_set + num_jointConstraint_active_set + num_boundingConstraint_active_set;
            
            % partial derivative of b_tilde w.r.t h
            %handle dynamics constraints
            db_tildedh_dyn = zeros(num_dynamicsConstraint_active_set,1);
            for j=1:num_dynamicsConstraint_active_set
                row = active_set(j);
                idx = (row-1)*dim + (1:dim);
                db_tildedh_dyn(j) = -(phi(row)/h^2-normal(:,row)'*dcdh(idx));%phi(row)/h^2 shows up because of v_min, it has a negative sign because
                % the first componint of bin_fqp is -bin, the sign is flipped
            end
            %handle joint constraints
            db_tildedh_joint = zeros(num_jointConstraint_active_set,1);
            for j=1:num_jointConstraint_active_set
                row = active_set(j+num_dynamicsConstraint_active_set);
                idx = num_active*dim+(row-num_active);
                db_tildedh_joint(j) = -(phi(row)/h^2 - dcdh(idx));%phi(row)/h^2 shows up because of v_min, it has a negative sign because
                % the first componint of bin_fqp is -bin, the sign is flipped
            end
            
            if ~isempty(dynamicsConstraint_active_set) && isempty(jointConstraint_active_set)
                db_tildedh = [db_tildedh_dyn;zeros(num_boundingConstraint_active_set, 1)];
            elseif ~isempty(dynamicsConstraint_active_set) && ~isempty(jointConstraint_active_set)
                db_tildedh = [db_tildedh_dyn;db_tildedh_joint;zeros(num_boundingConstraint_active_set, 1)];
            elseif isempty(dynamicsConstraint_active_set) && isempty(jointConstraint_active_set)
                db_tildedh = [zeros(num_boundingConstraint_active_set, 1)];
            elseif isempty(dynamicsConstraint_active_set) && ~isempty(jointConstraint_active_set)
                db_tildedh = [db_tildedh_joint;zeros(num_boundingConstraint_active_set, 1)];
            end
            
            dA_tildedq_dyn = zeros(num_dynamicsConstraint_active_set,num_params,num_q);
            db_tildedq_dyn = zeros(num_dynamicsConstraint_active_set,1,num_q);
            db_tildedv_dyn = zeros(num_dynamicsConstraint_active_set,1,num_q);
            dA_tildedq_joint = zeros(num_jointConstraint_active_set,num_params,num_q);
            db_tildedq_joint = zeros(num_jointConstraint_active_set,1,num_q);
            db_tildedv_joint = zeros(num_jointConstraint_active_set,1,num_q);
            dA_tildedq = zeros(num_constraint_active_set,num_params,num_q);
            dA_tildedv = zeros(num_constraint_active_set,num_params,num_q);
            db_tildedq = zeros(num_constraint_active_set,1,num_q);
            db_tildedv = zeros(num_constraint_active_set,1,num_q);
            for i=1:num_q % assume num_q == num_v
                % partial derivative of A_tilde and b_tilde w.r.t q and v
                %handle dynamics constraints
                for j=1:num_dynamicsConstraint_active_set
                    row = active_set(j);
                    idx = (row-1)*dim + (1:dim);
                    dA_tildedq_dyn(j,:,i) = -normal(:,row)'*(dAdq_tmp(idx,:,i)*V+A(idx,:)*dV(:,:,i));
                    db_tildedq_dyn(j,:,i) = n(row,i)/h + normal(:,row)'*dcdq(idx,i);% switch the sign to correct one 
                    db_tildedv_dyn(j,:,i) = normal(:,row)'*dcdv(idx,i);
                end
                %handle joint constraints
                for j=1:num_jointConstraint_active_set
                    row = active_set(j+num_dynamicsConstraint_active_set);
                    idx = num_active*dim+(row-num_active);
                    dA_tildedq_joint(j,:,i) = -(dAdq_tmp(idx,:,i)*V+A(idx,:)*dV(:,:,i));
                    db_tildedq_joint(j,:,i) = JL(row-num_active,i)/h + dcdq(idx,i);% switch the sign to correct one
                    db_tildedv_joint(j,:,i) = dcdv(idx,i);
                end
                
                if ~isempty(dynamicsConstraint_active_set) && isempty(jointConstraint_active_set)
                    dA_tildedq(:,:,i) = [dA_tildedq_dyn(:,:,i);zeros(num_boundingConstraint_active_set, num_params)];
                    db_tildedq(:,:,i) = [db_tildedq_dyn(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                    db_tildedv(:,:,i) = [db_tildedv_dyn(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                elseif ~isempty(dynamicsConstraint_active_set) && ~isempty(jointConstraint_active_set)
                    dA_tildedq(:,:,i) = [dA_tildedq_dyn(:,:,i);dA_tildedq_joint(:,:,i);zeros(num_boundingConstraint_active_set, num_params)];
                    db_tildedq(:,:,i) = [db_tildedq_dyn(:,:,i);db_tildedq_joint(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                    db_tildedv(:,:,i) = [db_tildedv_dyn(:,:,i);db_tildedv_joint(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                elseif isempty(dynamicsConstraint_active_set) && isempty(jointConstraint_active_set)
                    dA_tildedq(:,:,i) = [zeros(num_boundingConstraint_active_set, num_params)];
                    db_tildedq(:,:,i) = [zeros(num_boundingConstraint_active_set, 1)];
                    db_tildedv(:,:,i) = [zeros(num_boundingConstraint_active_set, 1)];
                elseif isempty(dynamicsConstraint_active_set) && ~isempty(jointConstraint_active_set)
                    dA_tildedq(:,:,i) = [dA_tildedq_joint(:,:,i);zeros(num_boundingConstraint_active_set, num_params)];
                    db_tildedq(:,:,i) = [db_tildedq_joint(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                    db_tildedv(:,:,i) = [db_tildedv_joint(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                end
                assert(length(dA_tildedq(:,1,i)) == length(active_set));
            end
             
            % partial derivative of b_tilde w.r.t u
            db_tildedu_dyn = zeros(num_dynamicsConstraint_active_set,1,obj.num_u);
            db_tildedu_joint = zeros(num_jointConstraint_active_set,1,obj.num_u);
            db_tildedu = zeros(num_constraint_active_set,1,obj.num_u);
            for i=1:obj.num_u
                %handle dynamics constraints
                for j=1:num_dynamicsConstraint_active_set
                    row = active_set(j);
                    idx = (row-1)*dim + (1:dim);
                    db_tildedu_dyn(j,:,i) = normal(:,row)'*dcdu(idx,i);
                end
                %handle joint constraints
                for j=1:num_jointConstraint_active_set
                    row = active_set(j+num_dynamicsConstraint_active_set);
                    idx = num_active*dim+(row-num_active);
                    db_tildedu_joint(j,:,i) = dcdu(idx,i);
                end
                
                if ~isempty(dynamicsConstraint_active_set) && isempty(jointConstraint_active_set)
                    db_tildedu(:,:,i) = [db_tildedu_dyn(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                elseif ~isempty(dynamicsConstraint_active_set) && ~isempty(jointConstraint_active_set)
                    db_tildedu(:,:,i) = [db_tildedu_dyn(:,:,i);db_tildedu_joint(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                elseif isempty(dynamicsConstraint_active_set) && isempty(jointConstraint_active_set)
                    db_tildedu(:,:,i) = [zeros(num_boundingConstraint_active_set, 1)];
                elseif isempty(dynamicsConstraint_active_set) && ~isempty(jointConstraint_active_set)
                    db_tildedu(:,:,i) = [db_tildedu_joint(:,:,i);zeros(num_boundingConstraint_active_set, 1)];
                end
            end
             
            %% partial derivative of KKT matrix blocks
            dGdq = zeros(num_params,num_params,num_q);
            dEdq = zeros(num_params,num_constraint_active_set,num_q);
            dFdq = zeros(num_constraint_active_set,num_constraint_active_set,num_q);
            for i=1:num_q
                % partial dervative w.r.t q
                M = Qinv*dAdq(:,:,i)*Qinv;%Qinv*reshape(- matGradMult(reshape(dAdq(:,:,i),[],1),Qinv),size(Qinv,1),[])
                dQinvdq(:,:,i) = -M;
                N = pinv(Ain_fqp_active*Qinv*Ain_fqp_active');
                dGdq(:,:,i) = -M + M*Ain_fqp_active'*N*Ain_fqp_active*Qinv + Qinv*Ain_fqp_active'*N*Ain_fqp_active*M ...
                    -Qinv*Ain_fqp_active'*N*Ain_fqp_active*M*Ain_fqp_active'*N*Ain_fqp_active*Qinv ...
                    -Qinv*dA_tildedq(:,:,i)'*N*Ain_fqp_active*Qinv - Qinv*Ain_fqp_active'*N*dA_tildedq(:,:,i)*Qinv ...
                    +Qinv*Ain_fqp_active'*N*dA_tildedq(:,:,i)*Qinv*Ain_fqp_active'*N*Ain_fqp_active*Qinv ...
                    +Qinv*Ain_fqp_active'*N*Ain_fqp_active*Qinv*dA_tildedq(:,:,i)'*N*Ain_fqp_active*Qinv;
                dEdq(:,:,i) = - M*Ain_fqp_active'*N + Qinv*Ain_fqp_active'*N*Ain_fqp_active*M*Ain_fqp_active'*N ...
                    + Qinv*dA_tildedq(:,:,i)'*N - Qinv*Ain_fqp_active'*N*dA_tildedq(:,:,i)*Qinv*Ain_fqp_active'*N ...
                    - Qinv*Ain_fqp_active'*N*Ain_fqp_active*Qinv*dA_tildedq(:,:,i)'*N;
                dEdq(:,:,i) = - dEdq(:,:,i);
                % again, Ain in Ain_fqp_active is negative, thus the sign of E is negative since there are three
                % Ain_fqp_active multiplied in E.
                
                % not used, can be used for slack variable analytical solution
                dFdq(:,:,i) = N*Ain_fqp_active*M*Ain_fqp_active'*N - N*dA_tildedq(:,:,i)*Qinv*Ain_fqp_active'*N ...
                    -N*Ain_fqp_active*Qinv*dA_tildedq(:,:,i)'*N;
                
                % partial dervative w.r.t v
                % dGdv = 0, dEdv = 0, dFdv = 0.
                % partial dervative w.r.t u
                % dGdu = 0, dEdu = 0, dFdu = 0.
            end
            
            %note that lambda = f = S_weighting*V*result_qp (there is coefficient S_weighting*V)
            %% partial derivative of lambda w.r.t. h, q, v, and u
            dlambdadh = S_weighting*V*(- G*dbdh - E*db_tildedh);
            
            dlambdadq = zeros(num_full_dim,num_q);
            dlambdadv = zeros(num_full_dim,num_q);
            for i=1:num_q
                dlambdadq(:,i) =  S_weighting*V*(- dGdq(:,:,i)*b - G*dbdq(:,i) - dEdq(:,:,i)*bin_fqp_active - E*db_tildedq(:,:,i)) ...
                    + S_weighting*dV(:,:,i)*result_analytical;
                dlambdadv(:,i) =  S_weighting*V*(- G*dbdv(:,i) - E*db_tildedv(:,:,i));
            end
            
            dlambdadu = zeros(num_full_dim,obj.num_u);
            for i=1:obj.num_u
                dlambdadu(:,i) = S_weighting*V*(- G*dbdu(:,i) - E*db_tildedu(:,:,i));
            end
            
            dlambda = [dlambdadh, dlambdadq, dlambdadv, dlambdadu];
            
            % Find quaternion indices
            quat_bodies = obj.manip.body([obj.manip.body.floating] == 2);
            quat_positions = [quat_bodies.position_num];
            for i=1:size(quat_positions,2)
                quat_dot = qdn(quat_positions(4:7,i));
                if norm(quat_dot) > 0
                    % Update quaternion by following geodesic
                    qn(quat_positions(4:7,i)) = q(quat_positions(4:7,i)) + quat_dot/norm(quat_dot)*tan(norm(h*quat_dot));
                    qn(quat_positions(4:7,i)) = qn(quat_positions(4:7,i))/norm(qn(quat_positions(4:7,i)));
                end
            end
            
            xdn = [qn;vn];
            
            % compute gradients
            if (nargout>1)
                if isempty(lambda)
                    dqdn = dwvn;
                else
                    dqdn = matGradMult(dMvn,lambda) + Mvn*dlambda + dwvn;
                end
                df = [ [qdn, eye(num_q), zeros(num_q,num_q+obj.num_u)]+h*dqdn; dqdn ];%[Ye: +h*dqdn part miss a vToqdot matrix]
            end
             
            %xdn = reshape(Mvn,[],1);
            %df = dMvn;
            %xdn = reshape(H,[],1);%correct
            %df = [zeros(size(dH,1),1), dH, zeros(size(dH,1),obj.num_u)];%correct
            %xdn = reshape(Hinv,[],1);%correct
            %df = [zeros(numel(Hinv),1),reshape(Hinv*reshape(- matGradMult(dH(:,1:num_q),Hinv),num_q,[]),numel(Hinv),[]),zeros(numel(Hinv),num_v+obj.num_u)];%correct
            %xdn = reshape(J_contact,[],1);
            %xdn = reshape(J',[],1);%correct
            %df = [zeros(numel(J),1), dJtranspose_qp, zeros(numel(J),num_v+obj.num_u)];%correct

            %A_length = size(A,1)*size(A,2);%correct at 1e-4 level
            %xdn = reshape(A, A_length,1);%correct
            %df = [zeros(A_length,1), reshape(dAdq_tmp,A_length,14), zeros(A_length,14), zeros(A_length,8)];%correct

            %A_length = size(V'*(A+R)*V,1)*size(V'*(A+R)*V,2);%correct at 1e-4 level
            %xdn = reshape(V'*(A+R)*V, A_length,1);%correct
            %df = [zeros(A_length,1), reshape(dAdq,A_length,14), zeros(A_length,14), zeros(A_length,8)];%correct
            
            %V_length = size(V,1)*size(V,2);%correct
            %xdn = reshape(V, V_length,1);%correct
            %df = [zeros(V_length,1), reshape(dV,V_length,14), zeros(V_length,14), zeros(V_length,8)];%correct

            %reshape(Hinv*reshape(dJtranspose_qp - matGradMult(dH(:,1:num_q),Hinv*J'),num_q,[]),numel(Mvn),[]),zeros(numel(Mvn),num_v+obj.num_u)];

            %xdn = reshape(J',[],1);
            %df = [zeros(length(xdn),1), dJtranspose_qp, zeros(length(xdn),14+8)];
            
            %xdn = lambda;%correct with a medium-level difference
            %df = dlambda;%correct with a medium-level difference
            %xdn = b;% correct
            %df = [dbdh, dbdq, dbdv, dbdu];% correct
            
            %Ain_length = size(Ain_fqp_active,1)*size(Ain_fqp_active,2);%correct at 1e-6 level
            %xdn = reshape(Ain_fqp_active,Ain_length,[]);%correct
            %df = [zeros(Ain_length,1), reshape(dA_tildedq,Ain_length,14), reshape(dA_tildedv,Ain_length,14), zeros(Ain_length,8)];%correct
            
            %bin_length = size(bin_fqp_active,1);%correct, only a few elements have large deviations.
            %xdn = bin_fqp_active;%correct, only a few elements have large deviations.
            %df = [db_tildedh, reshape(db_tildedq,bin_length,num_q), reshape(db_tildedv,bin_length,num_q), reshape(db_tildedu,bin_length,obj.num_u)];%correct
            
            %Q_length = size(Q,1)*size(Q,2);%correct with constant R value, sensitive to nonlinear R
            %xdn = reshape(Q,[],1);
            %df = [zeros(Q_length,1), reshape(0.5*dAdq,Q_length,14), zeros(Q_length,14), zeros(Q_length,8)];
            
            %Qinv_length = size(Qinv,1)*size(Qinv,2);%correct, if using linear scaling for R, accuracy at the 1e-4 level if 1e0 constant in Q. (accuracy at the 1e-1 level if 1e-2 constant in Q.)
            % if using nonlinear scaling for R (with nonlinear coeff 20), accuracy at the 1e-6 level if 1e0 constant in Q. (accuracy at the 1e-2 level if 1e-2 constant in Q.)
            %xdn = reshape(Qinv,[],1);
            %df = [zeros(Qinv_length,1), reshape(dQinvdq,Qinv_length,14), zeros(Qinv_length,14), zeros(Qinv_length,8)];
            %df = [zeros(numel(Qinv),1),reshape(Qinv*reshape(- matGradMult(reshape(dAdq(:,:,:),[],num_q),Qinv),size(Qinv,1),[]),numel(Qinv),[]),zeros(numel(Qinv),num_v+obj.num_u)];%correct

            %E_length = size(E,1)*size(E,2);%accuracy depends on constant term in Q matrix
            %xdn = reshape(E,E_length,[]);
            %df = [zeros(E_length,1), reshape(dEdq,E_length,14), zeros(E_length,14), zeros(E_length,8)];
            
            %F_length = size(F,1)*size(F,2);%accuracy depends on constant term in Q matrix
            %xdn = reshape(F,F_length,[]);
            %df = [zeros(F_length,1), reshape(dFdq,F_length,14), zeros(F_length,14), zeros(F_length,8)];
            
            %G_length = size(G,1)*size(G,2);%accuracy depends on constant term in Q matrix
            %xdn = reshape(G,G_length,[]);
            %df = [zeros(G_length,1), reshape(dGdq,G_length,14), zeros(G_length,14), zeros(G_length,8)];
            
            %xdn = wvn;%correct
            %df = dwvn;%correct
            %xdn = tau; 
            %df = dtau;
            %xdn = C;
            %df = [zeros(14,1), dC, zeros(14,8)];
            
            for i=1:length(obj.sensor)
                if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
                    if (nargout>1)
                        [obj,xdn_sensor,df_sensor] = update(obj.sensor{i},obj,t,x,u);
                    else
                        [obj,xdn_sensor] = update(obj.sensor{i},obj,t,x,u);
                    end
                    xdn = [xdn;xdn_sensor];
                    if (nargout>1)
                        df = [df; df_sensor];
                    end
                end
            end
            
            %Final comment of gradient debugging
            %dMvn and dwvn are correct
            %dbdh, dbdq, dbdv, dbdu are correct
            %Ain_fqp_active and dA_tildedq match
            %bin_fqp_active and db_tildedq match
            %dlambdadv and dlambdadu are correct.
            %The only numerical issue comes from dGdq(9:14), sometimes
            %the deviation from the numerical gradient is slightly large, but
            %overall reasonable, each individual value is close to each
            %other and makes sense. dGdq(1:8) is correct. No error there.
            %The numerical difference is slightly magnified in xdn and df.
            %By checking dimension by dimension, the values match pretty good.
            %A special part to pay attention to: the first dimension is
            %gradient w.r.t. h, in this case, the gradient is non-zero.
            %However, in the update() function, the gradient is w.r.t.
            %t. There its graident component is zero. This is the
            %difference between our QP solver and LCP solver.
            %Using gevl() to compaure the values, the values match quite well
            %(all the differences are numerically tolerable). Only a
            %few large deviations (the numerical values are very large, does
            %not make sense), overall, the analytical gradient are more
            %reliable.
        end
        
        function [xdn,df] = update(obj,t,x,u) 
            X0 = [t;x;u];%note that the first input is sigma point index
            [xdn,df] = solveQP(obj,X0);
            return;
            
            %t = index;
            
            %set the gradient w.r.t index to be zero.
            %df(:,1) = zeros(size(df,1),1);
            %xdn_numeric = xdn;
            %df_numeric = df;
            % fun = @(X0) solveQP(obj,X0);
            % DerivCheck(fun, X0)
            %
            % [xdn,df] = solveQP(obj,X0);
            %
            %[xdn_numeric,df_numeric] = geval(@(X0) solveQP(obj,X0),X0,struct('grad_method','numerical'));
            % valuecheck(xdn,xdn_numeric,1e-5);
            % valuecheck(df,df_numeric,1e-5);
            % [xdn_QP,df_QP] = solveQP(obj,X0);
            %
            % % gradient check lambda and dlambda
            % X0 = [t;x;u];
            % X0 = X0 + randn(size(X0))*0.1;
            
            function DerivCheck(funptr, X0, ~, varargin)
                
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %`
                %  Checks the analytic gradient of a function 'funptr' at a point X0, and
                %  compares to numerical gradient.  Useful for checking gradients computed
                %  for fminunc and fmincon.
                %
                %  Call with same arguments as you would call for optimization (fminunc).
                %
                % $id$
                
                [~, JJ] = feval(funptr, X0, varargin{:});  % Evaluate function at X0
                
                % Pick a random small vector in parameter space
                tol = 1e-6;  % Size of numerical step to take
                rr = sqrt(eps(X0));%randn(length(X0),1)*tol;  % Generate small random-direction vector
                
                % Evaluate at symmetric points around X0
                [f1, JJ1] = feval(funptr, X0-rr/2, varargin{:});  % Evaluate function at X0
                [f2, JJ2] = feval(funptr, X0+rr/2, varargin{:});  % Evaluate function at X0
                
                % Print results
                fprintf('Derivs: Analytic vs. Finite Diff = [%.12e, %.12e]\n', sum(sum(JJ*rr)), sum(sum(f2-f1)));
                %dd =  dot(rr, JJ)-f2+f1
                dd =  sum(sum(JJ*rr))-sum(sum(f2-f1))
            end
            
            if (nargout>1)
                [obj,z,Mvn,wvn,dz,dMvn,dwvn] = solveLCP(obj,t,x,u);
            else
                [obj,z,Mvn,wvn] = solveLCP(obj,t,x,u);
            end
            
            num_q = obj.manip.num_positions;
            q=x(1:num_q); v=x((num_q+1):end);
            h = obj.timestep;
            
            if isempty(z)
                vn = wvn;
            else
                vn = Mvn*z + wvn;
            end
            
            kinsol = obj.manip.doKinematics(q);
            vToqdot = obj.manip.vToqdot(kinsol);
            qdn = vToqdot*vn;
            qn = q+ h*qdn;
            % Find quaternion indices
            quat_bodies = obj.manip.body([obj.manip.body.floating] == 2);
            quat_positions = [quat_bodies.position_num];
            for i=1:size(quat_positions,2)
                quat_dot = qdn(quat_positions(4:7,i));
                if norm(quat_dot) > 0
                    % Update quaternion by following geodesic
                    qn(quat_positions(4:7,i)) = q(quat_positions(4:7,i)) + quat_dot/norm(quat_dot)*tan(norm(h*quat_dot));
                    qn(quat_positions(4:7,i)) = qn(quat_positions(4:7,i))/norm(qn(quat_positions(4:7,i)));
                end
            end
            xdn = [qn;vn];
            
            if (nargout>1)  % compute gradients
                if isempty(z)
                    dqdn = dwvn;
                else
                    dqdn = matGradMult(dMvn,z) + Mvn*dz + dwvn;
                end
                df = [ [zeros(num_q,1), eye(num_q), zeros(num_q,num_q+obj.num_u)]+h*dqdn; dqdn ];
            end
            
            for i=1:length(obj.sensor)
                if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
                    if (nargout>1)
                        [obj,xdn_sensor,df_sensor] = update(obj.sensor{i},obj,t,x,u);
                    else
                        [obj,xdn_sensor] = update(obj.sensor{i},obj,t,x,u);
                    end
                    xdn = [xdn;xdn_sensor];
                    if (nargout>1)
                        df = [df; df_sensor];
                    end
                end
            end
        end
        
        function [xdn,df] = solveLCP_new(obj,X0)
            
            t = X0(1);
            x = X0(2:29);
            u = X0(30:37);
            global timestep_updated
            
            if (nargout>1)
                [obj,z,Mvn,wvn,dz,dMvn,dwvn] = solveLCP(obj,t,x,u);
            else
                [obj,z,Mvn,wvn] = solveLCP(obj,t,x,u);
            end
            
            num_q = obj.manip.num_positions;
            q=x(1:num_q); v=x((num_q+1):end);
            
            h = timestep_updated;
            
            if isempty(z)
                vn = wvn;
            else
                vn = Mvn*z + wvn;
            end
            
            kinsol = obj.manip.doKinematics(q);
            vToqdot = obj.manip.vToqdot(kinsol);
            qdn = vToqdot*vn;
            qn = q+ h*qdn;
            % Find quaternion indices
            quat_bodies = obj.manip.body([obj.manip.body.floating] == 2);
            quat_positions = [quat_bodies.position_num];
            for i=1:size(quat_positions,2)
                quat_dot = qdn(quat_positions(4:7,i));
                if norm(quat_dot) > 0
                    % Update quaternion by following geodesic
                    qn(quat_positions(4:7,i)) = q(quat_positions(4:7,i)) + quat_dot/norm(quat_dot)*tan(norm(h*quat_dot));
                    qn(quat_positions(4:7,i)) = qn(quat_positions(4:7,i))/norm(qn(quat_positions(4:7,i)));
                end
            end
            xdn = [qn;vn];
            
            if (nargout>1)  % compute gradients
                if isempty(z)
                    dqdn = dwvn;
                else
                    dqdn = matGradMult(dMvn,z) + Mvn*dz + dwvn;
                end
                df = [ [zeros(num_q,1), eye(num_q), zeros(num_q,num_q+obj.num_u)]+h*dqdn; dqdn];%[Ye: +h*dqdn part miss a vToqdot matrix]
            end
            
            for i=1:length(obj.sensor)
                if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
                    if (nargout>1)
                        [obj,xdn_sensor,df_sensor] = update(obj.sensor{i},obj,t,x,u);
                    else
                        [obj,xdn_sensor] = update(obj.sensor{i},obj,t,x,u);
                    end
                    xdn = [xdn;xdn_sensor];
                    if (nargout>1)
                        df = [df; df_sensor];
                    end
                end
            end
        end
        
        function [xdn,df] = update2(obj,t,x,u)
            %if obj.update_convex && nargout>1
            %t
            global timestep_updated
            global x_initial
            %timestep_updated = 5e-4;
            % this is the key part.
            %if t == 0
            %    x = x_initial;
            %end
            X0 = [t;x;u];
            
            [xdn,df] = solveLCP_new(obj,X0);
            %tElapsed = toc(tStart);
            %xdn_QP_vec = [xdn_QP_vec,xdn];
            
            %return;
            %disp('finish solveQP QP')
            
            %% add gradient check
            %
            if X0(1) > 0.05
                fun = @(X0) solveLCP_new(obj,X0);
                DerivCheck(fun, X0)
                
                [xdn_numeric,df_numeric] = geval(@(X0) solveLCP_new(obj,X0),X0,struct('grad_method','numerical'));
                valuecheck(xdn,xdn_numeric,1e-5);
                valuecheck(df,df_numeric,1e-5);
            end
            
            function DerivCheck(funptr, X0, ~, varargin)
                
                % DerivCheck(funptr, X0, opts, arg1, arg2, arg3, ....);
                %`
                %  Checks the analytic gradient of a function 'funptr' at a point X0, and
                %  compares to numerical gradient.  Useful for checking gradients computed
                %  for fminunc and fmincon.
                %
                %  Call with same arguments as you would call for optimization (fminunc).
                %
                % $id$
                
                [~, JJ] = feval(funptr, X0, varargin{:});  % Evaluate function at X0
                
                % Pick a random small vector in parameter space
                tol = 1e-6;  % Size of numerical step to take
                rr = sqrt(eps(X0));%randn(length(X0),1)*tol;  % Generate small random-direction vector
                
                % Evaluate at symmetric points around X0
                [f1, JJ1] = feval(funptr, X0-rr/2, varargin{:});  % Evaluate function at X0
                [f2, JJ2] = feval(funptr, X0+rr/2, varargin{:});  % Evaluate function at X0
                
                % Print results
                fprintf('Derivs: Analytic vs. Finite Diff = [%.12e, %.12e]\n', sum(sum(JJ*rr)), sum(sum(f2-f1)));
                %dd =  dot(rr, JJ)-f2+f1
                dd =  sum(sum(JJ*rr))-sum(sum(f2-f1))
            end
        end
        
        function hit = cacheHit(obj,t,x,u,num_args_out)
            hit = (t==obj.LCP_cache.data.t && all(x==obj.LCP_cache.data.x) && ...
                all(u==obj.LCP_cache.data.u) && num_args_out <= obj.LCP_cache.data.nargout);
        end
        
        function [obj, z, Mqdn, wqdn] = solveMexLCP(obj, t, x, u)
            num_q = obj.manip.num_positions;
            q=x(1:num_q);
            v=x(num_q+(1:obj.manip.num_velocities));
            kinsol = doKinematics(obj, q, v);
            [H,C,B] = manipulatorDynamics(obj.manip, q, v);
            [phiC,normal,d,xA,xB,idxA,idxB,mu,n,D] = obj.manip.contactConstraints(kinsol, obj.multiple_contacts);
            [z, Mqdn, wqdn, possible_contact_indices, possible_jointlimit_indices] = solveLCPmex(obj.manip.mex_model_ptr, kinsol.mex_ptr, u, phiC, n, D, obj.timestep, obj.z_inactive_guess_tol, obj.LCP_cache.data.z, H, C, B, obj.enable_fastqp);
            possible_contact_indices = logical(possible_contact_indices);
            contact_data.normal = normal(:,possible_contact_indices);
            
            for i=1:length(d)
                contact_data.d{i} = d{i}(:,possible_contact_indices);
            end
            
            contact_data.xA = xA(:,possible_contact_indices);
            contact_data.xB = xB(:,possible_contact_indices);
            contact_data.idxA = idxA(possible_contact_indices);
            contact_data.idxB = idxB(possible_contact_indices);
            obj.LCP_cache.data.contact_data = contact_data;
            obj.LCP_cache.data.z = z;
            obj.LCP_cache.data.possible_limit_indices = logical(possible_jointlimit_indices)';
        end
        
        function [obj,z,Mvn,wvn,dz,dMvn,dwvn] = solveLCP(obj,t,x,u)
            % if (nargout<5 && obj.gurobi_present && obj.manip.only_loops && obj.manip.mex_model_ptr~=0 && ~obj.position_control)
            %     [obj,z,Mvn,wvn] = solveMexLCP(obj,t,x,u);
            %     return;
            % end
             
            %gravity compensation for Kuka arm
            Nq = obj.getNumPositions();
            Nq_arm = 8;
            Nu = obj.getNumInputs();
            Nv = obj.getNumVelocities();
            Nx = Nq+Nv;
            [H,C,B] = manipulatorDynamics(obj,x(1:Nq),zeros(Nv,1));
            u = B(1:Nq_arm,:)\C(1:Nq_arm);
            % hard coding open-loop torque modifications
            u(2) = u(2) + 7;
            u(4) = u(4) - 7;
            u(8) = -40;
            
            % global active_set_fail_count
            % do LCP time-stepping
            % todo: implement some basic caching here
            if cacheHit(obj,t,x,u,nargout)
                z = obj.LCP_cache.data.z;
                Mvn = obj.LCP_cache.data.Mqdn;
                wvn = obj.LCP_cache.data.wqdn;
                if nargout > 4
                    dz = obj.LCP_cache.data.dz;
                    dMvn = obj.LCP_cache.data.dMqdn;
                    dwvn = obj.LCP_cache.data.dwqdn;
                end
            else
                
                obj.LCP_cache.data.t = t;
                obj.LCP_cache.data.x = x;
                obj.LCP_cache.data.u = u;
                obj.LCP_cache.data.nargout = nargout;
                
                num_q = obj.manip.getNumPositions;
                num_v = obj.manip.getNumVelocities;
                q=x(1:num_q); v=x(num_q+(1:num_v));
                
                kinematics_options.compute_gradients = nargout > 4;
                kinsol = doKinematics(obj, q, [], kinematics_options);
                vToqdot = obj.manip.vToqdot(kinsol);
                qd = vToqdot*v;
                h = obj.timestep;
                
                if (nargout<5)
                    [H,C,B] = manipulatorDynamics(obj.manip,q,v);
                    if (obj.num_u>0 && ~obj.position_control)
                        tau = B*u - C;
                    else
                        tau = -C;
                    end
                else
                    [H,C,B,dH,dC,dB] = manipulatorDynamics(obj.manip,q,v);
                    if (obj.num_u>0 && ~obj.position_control)
                        tau = B*u - C;
                        dtau = [zeros(num_v,1), matGradMult(dB,u) - dC, B];
                    else
                        tau = -C;
                        dtau = [zeros(num_v,1), -dC, zeros(size(B))];
                    end
                end
                
                if (obj.position_control)
                    pos_control_index = getActuatedJoints(obj.manip);
                    nL = 2*length(pos_control_index);
                else
                    nL = sum([obj.manip.joint_limit_min~=-inf;obj.manip.joint_limit_max~=inf]); % number of joint limits
                end
                nContactPairs = obj.manip.getNumContactPairs;
                nP = obj.manip.num_position_constraints;  % number of position constraints
                nV = obj.manip.num_velocity_constraints;
                Big = 1e20;
                
                % Set up the LCP:
                % z >= 0, Mz + w >= 0, z'*(Mz + w) = 0
                % for documentation below, use slack vars: s = Mz + w >= 0
                %
                % use qn = q + h*qdn
                % where H(q)*(qdn - qd)/h = B*u - C(q) + J(q)'*z
                %  or qdn = qd + H\(h*tau + J'*z)
                %  with z = [h*cL; h*cP; h*cN; h*beta{1}; ...; h*beta{mC}; lambda]
                %
                % and implement equation (7) from Anitescu97, by collecting
                %   J = [JL; JP; n; D{1}; ...; D{mC}; zeros(nC,num_q)]
                
                % a "possible" index is one that is close enough to contact or the
                % limit to have a chance at being active.  now we don't even
                % construct the M and w matrices for contacts/limits that seem
                % impossible to reach within this timestep.
                possible_limit_indices = [];
                possible_contact_indices = [];
                possible_indices_changed = false;
                
                while (1)
                    if (nL > 0)
                        if (obj.position_control)
                            phiL = q(pos_control_index) - u;
                            JL = sparse(1:obj.manip.num_u,pos_control_index,1,obj.manip.num_u,obj.manip.num_positions);
                            phiL = [phiL;-phiL]; JL = [JL;-JL];
                            % dJ = 0 by default, which is correct here
                            dJL = zeros(length(phiL),num_q^2);
                        else
                            if (nargout>4)
                                [phiL,JL,dJL] = obj.manip.jointLimitConstraints(q);
                            else
                                [phiL,JL] = obj.manip.jointLimitConstraints(q);
                            end
                            if isempty(possible_limit_indices)
                                possible_limit_indices = (phiL + h*JL*qd) < obj.z_inactive_guess_tol;
                            end
                            nL = sum(possible_limit_indices);
                            
                            % phi_check and J_check are the "impossible indices"
                            % which get checked at the end of the method (to make sure
                            % that they did not somehow make contact or hit the limit)
                            phi_check = phiL(~possible_limit_indices);
                            J_check = JL(~possible_limit_indices,:);
                            phiL = phiL(possible_limit_indices);
                            JL = JL(possible_limit_indices,:);
                            if (nargout>4)
                                dJL = dJL(possible_limit_indices,:);
                            end
                            if isempty(obj.LCP_cache.data.possible_limit_indices) || any(obj.LCP_cache.data.possible_limit_indices~=possible_limit_indices)
                                possible_indices_changed = true;
                            end
                            obj.LCP_cache.data.possible_limit_indices=possible_limit_indices;
                        end
                    else
                        phi_check = zeros(0,1);
                        J_check = zeros(0,num_q);
                        JL = zeros(0,num_q^2);
                    end
                    
                    has_contacts = (nContactPairs > 0);
                    
                    if has_contacts
                        if (nargout>4)
                            [phiC,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = obj.contactConstraints_manual(kinsol, obj.multiple_contacts);
                        else
                            [phiC,normal,d,xA,xB,idxA,idxB,mu,n,D] = obj.contactConstraints_manual(kinsol, obj.multiple_contacts);
                        end
                        if ~isempty(phiC)
                            if isempty(possible_contact_indices)
                                possible_contact_indices = (phiC+h*n*qd) < obj.z_inactive_guess_tol;
                            end
                            
                            nC = sum(possible_contact_indices);
                            mC = length(D);
                            
                            J_size = [nL + nP + (mC+2)*nC,num_q];
                            J = zeros(J_size)*q(1); % *q(1) is for taylorvar
                            lb = zeros(nL+nP+(mC+2)*nC,1);
                            ub = Big*ones(nL+nP+(mC+2)*nC,1);
                            D = vertcat(D{:});
                            % just keep the likely contacts (and store data to check the unlikely):
                            phi_check = [phi_check;phiC(~possible_contact_indices)];
                            J_check = [J_check; n(~possible_contact_indices,:)];
                            phiC = phiC(possible_contact_indices);
                            n = n(possible_contact_indices,:);
                            D = D(repmat(possible_contact_indices,mC,1),:);
                            mu = mu(possible_contact_indices,:);
                            
                            if isempty(obj.LCP_cache.data.possible_contact_indices) || ...
                                    numel(obj.LCP_cache.data.possible_contact_indices)~= numel(possible_contact_indices) || ...
                                    any(obj.LCP_cache.data.possible_contact_indices~=possible_contact_indices)
                                possible_indices_changed = true;
                            end
                            
                            obj.LCP_cache.data.possible_contact_indices=possible_contact_indices;
                            
                            J(nL+nP+(1:nC),:) = n;
                            J(nL+nP+nC+(1:mC*nC),:) = D;
                            
                            if nargout>4
                                dJ = zeros(prod(J_size), num_q); % was sparse, but reshape trick for the transpose below didn't work
                                possible_contact_indices_found = find(possible_contact_indices);
                                n_size = [numel(possible_contact_indices), num_q];
                                col_indices = 1 : num_q;
                                dn = getSubMatrixGradient(reshape(dn, [], num_q), possible_contact_indices_found, col_indices, n_size);
                                dJ = setSubMatrixGradient(dJ, dn, nL+nP+(1:nC), 1 : J_size(2), J_size);
                                D_size = size(D);
                                dD_matrix = zeros(prod(D_size), num_q);
                                row_start = 0;
                                for i = 1 : mC
                                    dD_possible_contact = getSubMatrixGradient(dD{i}, possible_contact_indices_found, col_indices, n_size);
                                    dD_matrix = setSubMatrixGradient(dD_matrix, dD_possible_contact, row_start + (1 : nC), col_indices, D_size);
                                    row_start = row_start + nC;
                                end
                                dD = dD_matrix;
                                dJ = setSubMatrixGradient(dJ, dD, nL+nP+nC + (1 : mC * nC), col_indices, J_size);
                                dJ = reshape(dJ, [], num_q^2);
                            end
                            
                            contact_data.normal = normal(:,possible_contact_indices);
                            for i=1:length(d)
                                contact_data.d{i} = d{i}(:,possible_contact_indices);
                            end
                            contact_data.xA = xA(:,possible_contact_indices);
                            contact_data.xB = xB(:,possible_contact_indices);
                            contact_data.idxA = idxA(possible_contact_indices);
                            contact_data.idxB = idxB(possible_contact_indices);
                        else
                            has_contacts = false;
                        end
                    else
                        has_contacts = false;
                    end
                    
                    if ~has_contacts
                        mC=0;
                        nC=0;
                        J = zeros(nL+nP,num_q);
                        lb = zeros(nL+nP,1);
                        ub = Big*ones(nL+nP,1);
                        if (nargout>4)
                            dJ = zeros(nL+nP,num_q^2);
                        end
                        contact_data = struct();
                        phi_check = zeros(0,1);
                        J_check = zeros(0,num_q);
                    end
                    
                    obj.LCP_cache.data.contact_data = contact_data;
                    
                    if (nC+nL+nP+nV==0)  % if there are no possible contacts
                        z = [];
                        Mvn = [];
                        wvn = v + h*(H\tau);
                        obj.LCP_cache.data.z = z;
                        obj.LCP_cache.data.Mqdn = Mvn;
                        obj.LCP_cache.data.wqdn = wvn;
                        if (nargout>4)
                            Hinv = inv(H);
                            dz = zeros(0,1+obj.num_x+obj.num_u);
                            dwvn = [zeros(num_v,1+num_q),eye(num_v),zeros(num_v,obj.num_u)] + ...
                                h*Hinv*dtau - [zeros(num_v,1),h*Hinv*matGradMult(dH(:,1:num_q),Hinv*tau),zeros(num_v,num_q),zeros(num_v,obj.num_u)];
                            dMvn = zeros(0,1+obj.num_x+obj.num_u);
                            obj.LCP_cache.data.dz = dz;
                            obj.LCP_cache.data.dMqdn = dMvn;
                            obj.LCP_cache.data.dwqdn = dwvn;
                        end
                        return;
                    end
                    
                    % now that J is allocated (since it's size is known), apply JL
                    if (nL>0)
                        J(1:nL,:) = JL;
                        if nargout>=5, dJL(1:nL,:) = dJL; end
                    end
                    
                    %% Bilateral position constraints
                    if nP > 0
                        % write as
                        %   phiP + h*JP*qdn >= 0 && -phiP - h*JP*qdn >= 0
                        if (nargout<5)
                            [phiP,JP] = obj.manip.positionConstraints(q);
                        else
                            [phiP,JP,dJP] = obj.manip.positionConstraints(q);
                        end
                        J(nL+(1:nP),:) = JP;
                        lb(nL+(1:nP),1) = -Big;
                    end
                    
                    %% Bilateral velocity constraints
                    if nV > 0
                        error('not implemented yet');  % but shouldn't be hard
                    end
                    
                    M = zeros(nL+nP+(mC+2)*nC)*q(1);
                    w = zeros(nL+nP+(mC+2)*nC,1)*q(1);
                    
                    Hinv = inv(H);
                    wvn = v + h*Hinv*tau;
                    Mvn = Hinv*vToqdot'*J';
                    
                    if (nargout>4)
                        dM = zeros(size(M,1),size(M,2),1+num_q+num_v+obj.num_u);
                        dw = zeros(size(w,1),1+num_q+num_v+obj.num_u);
                        dwvn = [zeros(num_v,1+num_q),eye(num_v),zeros(num_v,obj.num_u)] + ...
                            h*Hinv*dtau - [zeros(num_v,1),h*Hinv*matGradMult(dH(:,1:num_q),Hinv*tau),zeros(num_v,num_q),zeros(num_v,obj.num_u)];
                        dJtranspose = reshape(permute(reshape(dJ,size(J,1),size(J,2),[]),[2,1,3]),numel(J),[]);
                        dMvn = [zeros(numel(Mvn),1),reshape(Hinv*reshape(dJtranspose - matGradMult(dH(:,1:num_q),Hinv*J'),num_q,[]),numel(Mvn),[]),zeros(numel(Mvn),num_v+obj.num_u)];
                    end
                    
                    % check gradients
                    %      xdn = Mqdn;
                    %      if (nargout>1)
                    %        df = dMqdn;
                    %        df = [zeros(prod(size(xdn)),1),reshape(dJ,prod(size(xdn)),[]),zeros(prod(size(xdn)),num_q+obj.num_u)];
                    %      end
                    %      return;
                    
                    %% Joint Limits:
                    % phiL(qn) is distance from each limit (in joint space)
                    % phiL_i(qn) >= 0, cL_i >=0, phiL_i(qn) * cL_I = 0
                    % z(1:nL) = cL (nL includes min AND max; 0<=nL<=2*num_q)
                    % s(1:nL) = phiL(qn) approx phiL + h*JL*qdn
                    if (nL > 0)
                        w(1:nL) = phiL + h*JL*vToqdot*wvn;
                        M(1:nL,:) = h*JL*vToqdot*Mvn;
                        if (nargout>4)
                            dJL = [zeros(numel(JL),1),reshape(dJL,numel(JL),[]),zeros(numel(JL),num_q+obj.num_u)];
                            if (obj.position_control)
                                dw(1:nL,:) = [zeros(size(JL,1),1),JL,zeros(size(JL,1),num_q),...
                                    [-1*ones(length(pos_control_index),obj.num_u);1*ones(length(pos_control_index),obj.num_u)]] + h*matGradMultMat(JL,wvn,dJL,dwvn);
                            else
                                dw(1:nL,:) = [zeros(size(JL,1),1),JL,zeros(size(JL,1),num_q+obj.num_u)] + h*matGradMultMat(JL,wvn,dJL,dwvn);
                            end
                            dM(1:nL,1:size(Mvn,2),:) = reshape(h*matGradMultMat(JL,Mvn,dJL,dMvn),nL,size(Mvn,2),[]);
                        end
                    end
                    
                    %% Bilateral Position Constraints:
                    % enforcing eq7, line 2
                    if (nP > 0)
                        w(nL+(1:nP)) = phiP + h*JP*vToqdot*wvn;
                        M(nL+(1:nP),:) = h*JP*vToqdot*Mvn;
                        if (nargout>4)
                            dJP = [zeros(numel(JP),1),reshape(dJP,numel(JP),[]),zeros(numel(JP),num_q+obj.num_u)];
                            dw(nL+(1:nP),:) = [zeros(size(JP,1),1),JP,zeros(size(JP,1),num_q+obj.num_u)] + h*matGradMultMat(JP,wvn,dJP,dwvn);
                            dM(nL+(1:nP),1:size(Mvn,2),:) = reshape(h*matGradMultMat(JP,Mvn,dJP,qMqdn),nP,size(Mvn,2),[]);
                        end
                    end
                    
                    %% Contact Forces:
                    % s(nL+nP+(1:nC)) = phiC+h*n*qdn  (modified (fixed?) from eq7, line 3)
                    % z(nL+nP+(1:nC)) = cN
                    % s(nL+nP+nC+(1:mC*nC)) = repmat(lambda,mC,1) + D*qdn  (eq7, line 4)
                    % z(nL+nP+nC+(1:mC*nC)) = [beta_1;...;beta_mC]
                    % s(nL+nP+(mC+1)*nC+(1:nC)) = mu*cn - sum_mC beta_mC (eq7, line 5)
                    % z(nL+nP+(mC+1)*nC+(1:nC)) = lambda
                    %
                    % The second set of conditions gives:
                    %   lambda_i >= the largest projection of the velocity vector
                    %   onto the d vectors (since lambda_i >= -(D*qdn)_i for all i,
                    % and by construction of d always having the mirror vectors,
                    %   lambda_i >= (D_qdn)_i
                    %
                    % The last eqs give
                    %  lambda_i > 0 iff (sum beta)_i = mu_i*cn_i
                    % where i is for the ith contact.
                    % Assume for a moment that mu_i*cn_i = 1, then (sum_beta)_i = 1
                    % is a constraint ensuring that sum_beta_j D_j is a convex
                    % combination of the D vectors (since beta_j is also > 0)
                    % So lambda_i >0 if forces for the ith contact are on the boundary of
                    % the friction cone (lambda_i could also be > 0 if some of the beta_j
                    % D_j's are internally canceling each other out)
                    %
                    % So one solution is
                    %  v_i = 0,
                    %  beta_i >= 0
                    %  lambda_i = 0,
                    %  sum_d beta_i < mu*cn_i
                    % and another solution is
                    %  v_i > 0  (sliding)
                    %  lambda_i = max_d (v_i)
                    %  beta_i = mu*cn_i only in the direction of the largest negative velocity
                    %  beta_i = 0 otherwise
                    % By virtue of the eqs of motion connecting v_i and beta_i, only one
                    % of these two can exist. (the first is actually many solutions, with
                    % beta_i in opposite directions canceling each other out).
                    if (nC > 0)
                        w(nL+nP+(1:nC)) = phiC+h*n*vToqdot*wvn;
                        M(nL+nP+(1:nC),:) = h*n*vToqdot*Mvn;
                        
                        w(nL+nP+nC+(1:mC*nC)) = D*vToqdot*wvn;
                        M(nL+nP+nC+(1:mC*nC),:) = D*vToqdot*Mvn;
                        M(nL+nP+nC+(1:mC*nC),nL+nP+(1+mC)*nC+(1:nC)) = repmat(eye(nC),mC,1);
                        
                        M(nL+nP+(mC+1)*nC+(1:nC),nL+nP+(1:(mC+1)*nC)) = [diag(mu), repmat(-eye(nC),1,mC)];
                        
                        if (nargout>4)
                            % n, dn, and dD were only w/ respect to q.  filled them out for [t,x,u]
                            dn = [zeros(size(dn,1),1),dn,zeros(size(dn,1),num_q+obj.num_u)];
                            dD = [zeros(numel(D),1),reshape(dD,numel(D),[]),zeros(numel(D),num_q+obj.num_u)];
                            
                            dw(nL+nP+(1:nC),:) = [zeros(size(n,1),1),n,zeros(size(n,1),num_q+obj.num_u)]+h*matGradMultMat(n,wvn,dn,dwvn);
                            dM(nL+nP+(1:nC),1:size(Mvn,2),:) = reshape(h*matGradMultMat(n,Mvn,dn,dMvn),nC,size(Mvn,2),[]);
                            
                            dw(nL+nP+nC+(1:mC*nC),:) = matGradMultMat(D,wvn,dD,dwvn);
                            dM(nL+nP+nC+(1:mC*nC),1:size(Mvn,2),:) = reshape(matGradMultMat(D,Mvn,dD,dMvn),mC*nC,size(Mvn,2),[]);
                        end
                    end
                    
                    QP_FAILED = true;
                    if ~possible_indices_changed && obj.enable_fastqp
                        z = zeros(nL+nP+(mC+2)*nC,1);
                        if isempty(obj.LCP_cache.data.z) || numel(obj.LCP_cache.data.z) ~= numel(lb)
                            z_inactive = true(nL+nP+(mC+2)*nC,1);
                            obj.LCP_cache.data.z = z;
                        else
                            z_inactive = obj.LCP_cache.data.z>lb+1e-8;
                            % use conservative guess of M_active to avoid occasional numerical issues
                            % when M*z_inactive + w > 1e-8 by a small amount
                        end
                        n_z_inactive = sum(z_inactive);
                        if n_z_inactive > 0
                            Aeq = M(z_inactive,z_inactive);
                            beq = -w(z_inactive);
                            Ain_fqp = -eye(n_z_inactive);
                            bin_fqp = -lb(z_inactive);
                            QblkDiag = {eye(n_z_inactive)};
                            fqp = zeros(n_z_inactive, 1);
                            [z_,info_fqp] = fastQPmex(QblkDiag,fqp,Ain_fqp,bin_fqp,Aeq,beq,[]);
                            QP_FAILED = info_fqp<0;
                            if ~QP_FAILED
                                z(z_inactive) = z_;
                                obj.LCP_cache.data.fastqp_active_set = find(abs(Ain_fqp*z_ - bin_fqp)<1e-6);
                                % we know:
                                % z(z_inactive) >= lb
                                % M(M_active, z_inactive)*z(z_inactive) + w(M_active) = 0
                                % M(M_inactive, z_inactive)*z(z_inactive) + w(M_inactive) >= 0
                                %
                                % check:
                                % z(z_inactive)'*(M(z_inactive, z_inactive)*z(z_inactive) + w(z_inactive)) == 0  % since it is not checked by QP
                                % M(z_active,z_inactive)*z(z_inactive)+w(z_active) >= 0  % since we're assuming that z(z_active) == 0
                                z_active = ~z_inactive(1:(nL+nP+nC));  % only look at joint limit, position, and normal variables since if cn_i = 0,
                                % then that's a solution and we don't care about the relative velocity \lambda_i
                                
                                QP_FAILED = (~isempty(w(z_active)) && any(M(z_active,z_inactive)*z(z_inactive)+w(z_active) < 0)) || ...
                                    any(((z(z_inactive)>lb(z_inactive)+1e-8) & (M(z_inactive, z_inactive)*z(z_inactive)+w(z_inactive)>1e-8))) || ...
                                    any(abs(z_'*(Aeq*z_ - beq)) > 1e-6);
                            end
                        end
                    end
                    
                    if QP_FAILED
                        % then the active set has changed, call pathlcp
                        try
                            z = pathlcp(M,w,lb,ub);
                        catch
                            %z = obj.LCP_cache.data.z;
                            %Mvn = obj.LCP_cache.data.Mqdn;
                            %wvn = obj.LCP_cache.data.wqdn;
                            %dz = obj.LCP_cache.data.dz;
                            %dMvn = obj.LCP_cache.data.dMqdn;
                            %dwvn = obj.LCP_cache.data.dwqdn;
                            disp('use previous iteration solution');
                            keyboard
                            return
                        end
                         
                        obj.LCP_cache.data.fastqp_active_set = [];
                    end
                    % for debugging
                    %cN = z(nL+nP+(1:nC))
                    %beta1 = z(nL+nP+nC+(1:nC))
                    %beta2 = z(nL+nP+2*nC+(1:nC))
                    %lambda = z(nL+nP+3*nC+(1:nC))
                    % end debugging
                    % more debugging
                    %        path_convergence_tolerance = 1e-6; % default according to http://pages.cs.wisc.edu/~ferris/path/options.pdf
                    %        assert(all(z>=0));
                    %        assert(all(M*z+w>=-path_convergence_tolerance));
                    %        valuecheck(z'*(M*z+w),0,path_convergence_tolerance);
                    % end more debugging
                    
                    if obj.lcmgl_contact_forces_scale>0
                        cN = z(nL+nP+(1:nC));
                        cartesian_force = repmat(cN',3,1).*contact_data.normal;
                        for i=1:mC/2  % d is mirrored in contactConstraints
                            beta = z(nL+nP+i*nC+(1:nC));
                            cartesian_force = cartesian_force + repmat(beta',3,1).*contact_data.d{i};
                            beta = z(nL+nP+(mC/2+i)*nC+(1:nC));
                            cartesian_force = cartesian_force - repmat(beta',3,1).*contact_data.d{i};
                        end
                        
                        lcmgl = drake.matlab.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton,'LCP contact forces');
                        for j=1:nC
                            point = forwardKin(obj.manip,kinsol,contact_data.idxA(j),contact_data.xA(:,j));
                            lcmgl.glColor3f(.4,.2,.4);
                            lcmgl.drawVector3d(point,-obj.lcmgl_contact_forces_scale*cartesian_force(:,j));
                            lcmgl.glColor3f(.2,.4,.2);
                            lcmgl.drawVector3d(point,obj.lcmgl_contact_forces_scale*cartesian_force(:,j));
                        end
                        lcmgl.switchBuffers();
                    end
                    
                    obj.LCP_cache.data.z = z;
                    obj.LCP_cache.data.Mqdn = Mvn;
                    obj.LCP_cache.data.wqdn = wvn;
                    
                    if (nargout>4)
                        % Quick derivation:
                        % The LCP solves for z given that:
                        % M(a)*z + q(a) >= 0
                        % z >= 0
                        % z'*(M(a)*z + q(a)) = 0
                        % where the vector inequalities are element-wise, and 'a' is a vector of  parameters (here the state x and control input u).
                        %
                        % Our goal is to solve for the gradients dz/da.
                        %
                        % First we solve the LCP to obtain z.
                        %
                        % Then, for all i where z_i = 0, then dz_i / da = 0.
                        % Call the remaining active constraints (where z_i >0)  Mbar(a), zbar, and  qbar(a).  then we have
                        % Mbar(a) * zbar + qbar(a) = 0
                        %
                        % and the remaining gradients are given by
                        % for all j, dMbar/da_j * zbar + Mbar * dzbar / da_j + dqbar / da_j = 0
                        % or
                        %
                        % dzbar / da_j =  -pinv(Mbar)*(dMbar/da_j * zbar + dqbar / da_j)
                        %
                        % Note that there may be multiple solutions to the above equation
                        %    so we use the pseudoinverse to select the min norm solution
                        
                        dz = zeros(size(z,1),1+obj.num_x+obj.num_u);
                        zposind = find(z>0);
                        if ~isempty(zposind)
                            Mbar = M(zposind,zposind);
                            dMbar = reshape(dM(zposind,zposind,:),numel(Mbar),[]);
                            zbar = z(zposind);
                            dwbar = dw(zposind,:);
                            dz(zposind,:) = -pinv(Mbar)*(matGradMult(dMbar,zbar) + dwbar);
                        end
                        obj.LCP_cache.data.dz = dz;
                        obj.LCP_cache.data.dMqdn = dMvn;
                        obj.LCP_cache.data.dwqdn = dwvn;
                    else
                        obj.LCP_cache.data.dz = [];
                        obj.LCP_cache.data.dMqdn = [];
                        obj.LCP_cache.data.dwqdn = [];
                    end
                     
                    penetration = phi_check + h*J_check*vToqdot*(Mvn*z + wvn) < 0;
                    if any(penetration)
                        % then add the constraint and run the entire loop again!
                        try
                            limits = sum(~possible_limit_indices);
                            possible_limit_indices(~possible_limit_indices) = penetration(1:limits);
                            possible_contact_indices(~possible_contact_indices) = penetration(limits+1:end);
                            obj.warning_manager.warnOnce('Drake:TimeSteppingRigidBodyManipulator_Kuka:ResolvingLCP','This timestep violated our assumptions about which contacts could possibly become active in one timestep.  Consider reducing your dt.  If it seems to happen a lot, then please let us know about it.');
                        catch
                            disp('penetration')
                            break;
                        end
                    else
                        break;
                    end
                end
            end
        end
        
        function obj = addSensor(obj,sensor)
            if isa(sensor,'RigidBodySensor')
                obj.manip = obj.manip.addSensor(sensor);
            else
                typecheck(sensor,'TimeSteppingRigidBodySensor');
                obj.sensor{end+1} = sensor;
            end
        end
        
        function [obj,frame_id] = addFrame(obj,frame)
            [obj.manip,frame_id] = obj.manip.addFrame(frame);
        end
        
        function varargout = pdcontrol(sys,Kp,Kd,index)
            if nargin<4, index=[]; end
            [pdff,pdfb] = pdcontrol(sys.manip,Kp,Kd,index);
            pdfb = setInputFrame(pdfb,sys.manip.getStateFrame());
            pdfb = setOutputFrame(pdfb,sys.getInputFrame());
            pdff = setOutputFrame(pdff,sys.getInputFrame());
            if nargout>1
                varargout{1} = pdff;
                varargout{2} = pdfb;
            else
                % note: design the PD controller with the (non time-stepping
                % manipulator), but build the closed loop system with the
                % time-stepping manipulator:
                varargout{1} = cascade(pdff,feedback(sys,pdfb));
            end
        end
        
    end
    
    methods  % pass through methods (to the manipulator)
        function B = getB(obj)
            B = getB(obj.manip);
        end
        
        function g = getGravity(obj)
            g = getGravity(obj.manip);
        end
        
        function num_q = getNumPositions(obj)
            num_q = obj.manip.num_positions;
        end
        
        function num_v = getNumVelocities(obj)
            num_v = obj.manip.getNumVelocities();
        end
        
        function num_JL = getNumJointLimitConstraints(obj)
            num_JL = obj.manip.getNumJointLimitConstraints();
        end
        
        function obj = setStateFrame(obj,fr)
            obj = setStateFrame@DrakeSystem(obj,fr);
            
            % make sure there is a transform defined to and from the
            % manipulator state frame.  (the trivial transform is the correct
            % one)
            if ~isempty(obj.manip) % this also gets called on the initial constructor
                mfr = getStateFrame(obj.manip);
                if isempty(findTransform(fr,mfr))
                    addTransform(fr,AffineTransform(fr,mfr,eye(obj.manip.num_x,obj.num_x),zeros(obj.manip.num_x,1)));
                end
                if isempty(findTransform(mfr,fr))
                    addTransform(mfr,AffineTransform(mfr,fr,eye(obj.num_x,obj.manip.num_x),zeros(obj.num_x,1)));
                end
            end
        end
        
        function obj = setTerrain(obj,varargin)
            obj.manip = setTerrain(obj.manip,varargin{:});
        end
        
        function terrain = getTerrain(obj)
            terrain = obj.manip.terrain;
        end
        
        function varargout = getTerrainHeight(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = getTerrainHeight(obj.manip,varargin{:});
        end
        
        function obj = setJointLimits(obj,varargin)
            obj.manip = setJointLimits(obj.manip,varargin{:});
        end
        
        function obj=addRobotFromURDF(obj,varargin)
            if obj.twoD
                w = warning('off','Drake:PlanarRigidBodyManipulator:UnsupportedContactPoints');
                warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            else
                w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            end
            obj.manip=addRobotFromURDF(obj.manip,varargin{:});
            obj=compile(obj);  % note: compiles the manip twice, but it's ok.
            warning(w);
        end
        
        function obj=addRobotFromSDF(obj,varargin)
            obj.manip=addRobotFromSDF(obj.manip,varargin{:});
            obj=compile(obj);  % note: compiles the manip twice, but it's ok.
        end
        
        function varargout = doKinematics(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=doKinematics(obj.manip,varargin{:});
        end
        
        function varargout = forwardKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=forwardKin(obj.manip,varargin{:});
        end
        
        function varargout = bodyKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=bodyKin(obj.manip,varargin{:});
        end
        
        function varargout = approximateIK(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=approximateIK(obj.manip,varargin{:});
        end
        
        function varargout = inverseKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=inverseKin(obj.manip,varargin{:});
        end
        
        function varargout = inverseKinPointwise(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinPointwise(obj.manip,varargin{:});
        end
        
        function varargout = inverseKinTraj(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinTraj(obj.manip,varargin{:});
        end
        
        function varargout = inverseKinWrapup(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinWrapup(obj.manip,varargin{:});
        end
        
        function varargout = findFixedPoint(obj,x0,varargin)
            varargout=cell(1,nargout);
            if isnumeric(x0)
                x0 = Point(obj.getStateFrame(),x0);
            end
            [varargout{:}]=findFixedPoint(obj.manip,x0,varargin{:});
            varargout{1} = varargout{1}.inFrame(obj.getStateFrame());
        end
        
        function varargout = collisionDetect(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=collisionDetect(obj.manip,varargin{:});
        end
        
        function varargout = collisionDetectTerrain(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=collisionDetectTerrain(obj.manip,varargin{:});
        end
        
        function [obj,id] = addStateConstraint(obj,con)
            % keep two copies of the constraints around ... :(
            % todo: re-evaluate whether that is really necessary
            [obj,id] = addStateConstraint@DrakeSystem(obj,con);
            [obj.manip,manip_id] = obj.manip.addStateConstraint(obj,con);
            assert(id==manip_id);
        end
        
        function obj = updateStateConstraint(obj,id,con)
            obj = updateStateConstraint@DrakeSystem(obj,id,con);
            obj.manip = updateStateConstraint(obj.manip,id,con);
        end
        
        function obj = removeAllStateConstraints(obj)
            obj = removeAllStateConstraints@DrakeSystem(obj);
            obj.manip = removeAllStateConstraints(obj.manip);
        end
        
        function varargout = positionConstraints(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = positionConstraints(obj.manip,varargin{:});
        end
        
        function varargout = velocityConstraints(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = velocityConstraints(obj.manip,varargin{:});
        end
        
        function varargout = manipulatorDynamics(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = manipulatorDynamics(obj.manip,varargin{:});
        end
        
        function varargout = contactConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraints(obj.manip,varargin{:});
        end
        
        function varargout = contactConstraintDerivatives(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraintDerivatives(obj.manip,varargin{:});
        end
        
        function varargout = contactConstraintsBV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraintsBV(obj.manip,varargin{:});
        end
        
        function varargout = pairwiseContactConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = pairwiseContactConstraints(obj.manip,varargin{:});
        end
        
        function varargout = pairwiseContactConstraintsBV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = pairwiseContactConstraintsBV(obj.manip,varargin{:});
        end
        
        function varargout = resolveConstraints(obj,x0,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = resolveConstraints(obj.manip,x0,varargin{:});
            varargout{1} = varargout{1}.inFrame(obj.getStateFrame());
        end
        
        function varargout = getMass(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getMass(obj.manip,varargin{:});
        end
        
        function varargout = getCOM(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getCOM(obj.manip,varargin{:});
        end
        
        function varargout = centerOfMassJacobianDotTimesV(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = centerOfMassJacobianDotTimesV(obj.manip,varargin{:});
        end
        
        function varargout = centroidalMomentumMatrixDotTimesV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = centroidalMomentumMatrixDotTimesV(obj.manip,varargin{:});
        end
        
        function varargout = centroidalMomentumMatrix(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = centroidalMomentumMatrix(obj.manip,varargin{:});
        end
        
        function varargout = parseBodyOrFrameID(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = parseBodyOrFrameID(obj.manip,varargin{:});
        end
        
        function joint_ind = findJointId(model,varargin)
            joint_ind = findJointId(model.manip,varargin{:});
        end
        
        function body_ind = findLinkId(model,varargin)
            body_ind = findLinkId(model.manip,varargin{:});
        end
        
        function indices = findPositionIndices(model, varargin)
            indices = findPositionIndices(model.manip,varargin{:});
        end
        
        function body = findLink(model,varargin)
            body = findLink(model.manip,varargin{:});
        end
        
        function frame_id = findFrameId(model,varargin)
            frame_id = findFrameId(model.manip,varargin{:});
        end
        
        function ancestor_bodies = findAncestorBodies(obj, body_index)
            ancestor_bodies = obj.manip.findAncestorBodies(body_index);
        end
        
        function [body_path, joint_path, signs] = findKinematicPath(obj, start_body, end_body)
            [body_path, joint_path, signs] = obj.manip.findKinematicPath(start_body, end_body);
        end
        
        function obj = weldJoint(obj,body_ind_or_joint_name,robot)
            if nargin>2
                obj.manip = weldJoint(obj.manip,body_ind_or_joint_name,robot);
            else
                obj.manip = weldJoint(obj.manip,body_ind_or_joint_name);
            end
            obj.dirty = true;
        end
        
        function body = getBody(model,varargin)
            body = getBody(model.manip,varargin{:});
        end
        
        function frame = getFrame(model,varargin)
            frame = getFrame(model.manip,varargin{:});
        end
        
        function str = getBodyOrFrameName(obj,varargin)
            str = obj.manip.getBodyOrFrameName(varargin{:});
        end
        
        function model = setBody(model,varargin)
            model.manip = setBody(model.manip,varargin{:});
            model.dirty = true;
        end
        
        function v = constructVisualizer(obj,varargin)
            v = constructVisualizer(obj.manip,varargin{:});
        end
        
        function getNumContacts(~)
            error('getNumContacts is no longer supported, in anticipation of alowing multiple contacts per body pair. Use getNumContactPairs for cases where the number of contacts is fixed');
        end
        
        function n=getNumContactPairs(obj)
            n = obj.manip.getNumContactPairs;
        end
        
        function c = getBodyContacts(obj,body_idx)
            c = obj.manip.body(body_idx).collision_geometry;
        end
        
        function addContactShapeToBody(varargin)
            errorDeprecatedFunction('addCollisionGeometryToBody');
        end
        
        function obj = addCollisionGeometryToBody(obj,varargin)
            obj.manip = addCollisionGeometryToBody(obj.manip,varargin{:});
        end
        
        function addVisualShapeToBody(varargin)
            errorDeprecatedFunction('addVisualGeometryToBody');
        end
        
        function obj = addVisualGeometryToBody(obj,varargin)
            obj.manip = addVisualGeometryToBody(obj.manip,varargin{:});
        end
        
        function addShapeToBody(varargin)
            errorDeprecatedFunction('addGeometryToBody');
        end
        
        function obj = addGeometryToBody(obj,varargin)
            obj.manip = addGeometryToBody(obj.manip,varargin{:});
        end
        
        function replaceContactShapesWithCHull(varargin)
            errorDeprecatedFunction('replaceCollisionGeometryWithConvexHull');
        end
        
        function obj = replaceCollisionGeometryWithConvexHull(obj,body_indices,varargin)
            obj.manip = replaceCollisionGeometryWithConvexHull(obj.manip,body_indices,varargin{:});
        end
        
        function getContactShapeGroupNames(varargin)
            errorDeprecatedFunction('getCollisionGeometryGroupNames');
        end
        
        function groups = getCollisionGeometryGroupNames(obj)
            groups = getCollisionGeometryGroupNames(obj.manip);
        end
        
        function f_friction = computeFrictionForce(obj,qd)
            f_friction = computeFrictionForce(obj.manip,qd);
        end
        
        function obj = removeCollisionGroups(obj,contact_groups)
            obj.manip = removeCollisionGroups(obj.manip,contact_groups);
        end
        
        function obj = removeCollisionGroupsExcept(obj,varargin)
            obj.manip = removeCollisionGroupsExcept(obj.manip,varargin{:});
        end
        
        function str = getLinkName(obj,body_ind)
            str = obj.manip.getLinkName(body_ind);
        end
        
        function link_names = getLinkNames(obj)
            link_names =  {obj.manip.body.linkname}';
        end
        
        function joint_names = getJointNames(obj)
            joint_names =  {obj.manip.body.jointname}';
        end
        
        function num_bodies = getNumBodies(obj)
            num_bodies = length(obj.manip.body);
        end
        
        function [jl_min, jl_max] = getJointLimits(obj)
            jl_min = obj.manip.joint_limit_min;
            jl_max = obj.manip.joint_limit_max;
        end
        
        function varargout = jointLimitConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = jointLimitConstraints(obj.manip,varargin{:});
        end
        
        function varargout = getZeroConfiguration(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getZeroConfiguration(obj.manip,varargin{:});
        end
        
        function index = getActuatedJoints(obj)
            index = getActuatedJoints(obj.manip);
        end
        
        function ptr = getMexModelPtr(obj)
            ptr = getMexModelPtr(obj.manip);
        end
        
        function [phi,Jphi] = closestDistance(obj,varargin)
            [phi,Jphi] = closestDistance(obj.manip,varargin{:});
        end
        
        function obj = addLinksToCollisionFilterGroup(obj,linknames,collision_fg_name,robotnums)
            obj.manip = addLinksToCollisionFilterGroup(obj.manip,linknames,collision_fg_name,robotnums);
        end
        
        function out = name(obj)
            out = obj.manip.name;
        end
        
        function fr = getParamFrame(model)
            fr = getParamFrame(model.manip);
        end
        
        function model = setParams(model,p)
            model.manip = setParams(model.manip,p);
        end
         
        function terrain_contact_point_struct = getTerrainContactPoints(obj,varargin)
            terrain_contact_point_struct = getTerrainContactPoints(obj.manip,varargin{:});
        end
        
        function varargout = terrainContactPositions(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = terrainContactPositions(obj.manip,varargin{:});
        end
        
        function varargout = terrainContactJacobianDotTimesV(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = terrainContactJacobianDotTimesV(obj.manip,varargin{:});
        end
        
        function distance = collisionRaycast(obj, kinsol, origin, point_on_ray, use_margins)
            if nargin < 5
                use_margins = true;
            end
            distance = collisionRaycast(obj.manip, kinsol, origin, point_on_ray, use_margins);
        end
    end
end