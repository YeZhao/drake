classdef KukaArm < TimeSteppingRigidBodyManipulator_Kuka
    
    properties
        hand_name = 'iiwa_link_ee';
        disturbance_type = 1; % 1 ee-force, 2-state error, 3-torque
        cylinder_id
        left_finger_id
        right_finger_id
        iiwa_link_7_id
        uncertainty_source
        uncertainty_source_default
        object_initial_position
        friction_coeff
        time_step
        uncertain_mu
        uncertain_mu_mean
        uncertain_mu_set
        uncertain_phi
        uncertain_ori
        uncertain_position_mean
        uncertain_position_set
        uncertain_orientation_mean
        uncertain_orientation_set
    end
    
    methods
        function obj = KukaArm(options)
            if nargin < 1
                options = struct();s
            end
            
            if ~isfield(options,'floating')
                options.floating = false;
            end
            if ~isfield(options,'urdf')
                options.urdf = 'urdf/iiwa14_fixed_gripper.urdf';
            end
            if ~isfield(options,'with_weight')
                options.with_weight = false;
            end
            if ~isfield(options,'with_box')
                options.with_box = false;
            end
            if ~isfield(options,'with_shelf')
                options.with_shelf = false;
            end
            if ~isfield(options,'with_shelf_and_boxes')
                options.with_shelf_and_boxes = false;
            end
            if options.with_weight
                options.urdf = 'urdf/iiwa14_fixed_gripper.urdf';
            end
            
            warning('off','Drake:RigidBodyManipulator:UnsupportedVelocityLimits');
            warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            warning('off','Drake:RigidBodyManipulator:WeldedLinkInd');
            warning('off','Drake:RigidBodyManipulator:UnsupportedJointLimits');
            obj = obj@TimeSteppingRigidBodyManipulator_Kuka(options.urdf,0.001,options);
            
            if options.with_weight
                options_hand.weld_to_link = findLinkId(obj,obj.hand_name);
                obj = obj.addRobotFromURDF('urdf/robotiq_simple.urdf', [0;0;0.099], [pi/2;0;0], options_hand);
            end
            if options.with_box
                obj = obj.addRobotFromURDF('urdf/box.urdf', [0.6;0;1.4], [0;0;0]);
            end
            if options.with_shelf
                obj = obj.addRobotFromURDF('urdf/shelf.urdf', [0.6;0;0.88], [0;0;0]);
                obj = obj.addRobotFromURDF('urdf/shelf.urdf', [0.6;0.2;1.08], [pi/2;0;0]);
                obj = obj.addRobotFromURDF('urdf/shelf.urdf', [0.8;0.0;1.08], [0;pi/2;0]);
            end
            if options.with_shelf_and_boxes
                %obj = obj.addRobotFromURDF('urdf/box.urdf', [-0.3;-.9;.9], [0;0;0]);
                obj = obj.addRobotFromURDF('urdf/box.urdf', [-0.3;0;1.5], [0;0;0]);
                obj = obj.addRobotFromURDF('urdf/shelf.urdf', [0.6;0;0.88], [0;0;0]);
                obj = obj.addRobotFromURDF('urdf/shelf.urdf', [0.6;0.2;1.08], [pi/2;0;0]);
                obj = obj.addRobotFromURDF('urdf/shelf.urdf', [0.8;0.0;1.08], [0;pi/2;0]);
            end
            
            obj = obj.removeCollisionGroupsExcept({'manip'});
            options.floating = true;
            %obj = obj.addRobotFromURDF('urdf/cylinder.urdf',[],[],options);
            %obj = obj.addRobotFromURDF('urdf/cylinder_small.urdf',[],[],options);
            obj = obj.addRobotFromURDF('urdf/cylinder_small_with_inner_sliding_pendulum.urdf',[],[],options);
            
            obj.cylinder_id = obj.findLinkId('cylinder');
            obj.left_finger_id = obj.findLinkId('left_finger');
            obj.right_finger_id = obj.findLinkId('iiwa_link_7+iiwa_link_ee+base_link+right_finger+iiwa_link_ee_kuka');
            obj.iiwa_link_7_id = obj.findLinkId('iiwa_link_7');
            
            obj = compile(obj);
            
        end
        
        function nw = getNumDisturbances(obj)
            switch obj.disturbance_type
                case 1
                    nw = 3;
                case 2
                    nw = obj.getNumContStates();
                case 3
                    nw = obj.getNumInputs();
                otherwise
                    error('Unknown disturbance type');
            end
        end
        
        function [f,df,d2f] = dynamics_w(obj,t,x,u,w)
            
            [f,df] = dynamics_w_ee(obj,t,x,u,w);
            
            if nargout == 3
                %Finite diff to get 2nd derivatives
                nx = length(x);
                nu = length(u);
                nw = length(w);
                
                Dx = 1e-6*eye(nx);
                Du = 1e-6*eye(nu);
                Dw = 1e-6*eye(nw);
                
                d2f = zeros(nx, 1+nx+nu+nw, 1+nx+nu+nw);
                for k = 1:nx
                    [~,df_p] = dynamics_w_ee(obj,t,x+Dx(:,k),u,w);
                    d2f(:,:,1+k) = df_p-df;
                end
                for k = 1:nu
                    [~,df_p] = dynamics_w_ee(obj,t,x,u+Du(:,k),w);
                    d2f(:,:,1+nx+k) = df_p-df;
                end
                for k = 1:nw
                    [~,df_p] = dynamics_w_ee(obj,t,x,u,w+Dw(:,k));
                    d2f(:,:,1+nx+nu+k) = df_p-df;
                end
                
                d2f = reshape(d2f,nx,(1+nx+nu+nw)*(1+nx+nu+nw));
                
            end
            
        end
        
        function [f,df] = dynamics_w_ee(obj,t,x,u,w)
            % w should be a 3x1 force vector in world coordinates acting at the
            % robot's end effector
            
            nq = obj.getNumPositions;
            q=x(1:nq);
            
            if (nargout>1)
                nx = obj.getNumStates;
                nu = obj.getNumInputs;
                
                kinsol = doKinematics(obj, q, [], struct('compute_gradients', true));
                [~,J,dJ] = obj.forwardKin(kinsol,findLinkId(obj,obj.hand_name),[0;0;0]);
                uw = J'*w;
                
                dJtw = zeros(nq,nq);
                for i=1:nq
                    dJtw(:,i) = dJ(:,(i-1)*nq+(1:nq))'*w;
                end
                
                [f,df] = dynamics(obj,t,x,u+uw);
                df_du = df(:,1+nx+(1:nu));
                df_dq = df(:,1+(1:nq)) + df_du*dJtw;
                df_dqd = df(:,1+nq+(1:nq));
                df_dw = df_du*J';
                
                df = [df(:,1),df_dq,df_dqd,df_du,df_dw];
            else
                kinsol = doKinematics(obj, q, []);
                [~,J] = obj.forwardKin(kinsol,findLinkId(obj,obj.hand_name),[0;0;0]);
                uw = J'*w;
                
                [f,df] = dynamics(obj,t,x,u+uw);
            end
        end
        
        function [f,df] = dynamics_w_x(obj,t,x,u,w)
            % w is a state error vector
            [f,df] = dynamics(obj,t,x+w,u);
            df = [df,df(:,1+(1:obj.getNumStates))];
        end
        
        %         function c = getZeroConfiguration(obj)
        %             c=zeros(obj.getNumStates,1);
        %         end
        
        function [f,df] = dynamics_w_u(obj,t,x,u,w)
            % u is a state error vector
            [f,df] = dynamics(obj,t,x,u+w);
            df = [df,df(:,1+obj.getNumStates+(1:obj.getNumInputs))];
        end
        
        function u0 = findTrim(obj,q0)
            Nq = obj.getNumPositions();
            Nq_arm = 8;
            Nu = obj.getNumInputs();
            Nv = obj.getNumVelocities();
            Nx = Nq+Nv;
            
            [H,C,B] = manipulatorDynamics(obj,q0,zeros(Nv,1));
            
            u0 = B(1:Nq_arm,:)\C(1:Nq_arm);
        end
        
        function [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = contactConstraints_old(obj,kinsol,allow_multiple_contacts,active_collision_options)
            
            % @retval phi (m x 1) Vector of gap function values (typically contact distance), for m possible contacts
            % @retval normal (3 x m) Contact normal vector in world coordinates, points from B to A
            % @retval d {k} (3 x m) Contact friction basis vectors in world coordinates, points from B to A
            % @retval xA (3 x m) The closest point on body A to contact with body B, relative to body A origin and in body A frame
            % @retval xB (3 x m) The closest point on body B to contact with body A, relative to body B origin and in body B frame
            % @retval idxA (m x 1) The index of body A. 0 is the special case for the environment/terrain
            % @retval idxB (m x 1) The index of body B. 0 is the special case for the environment/terrain
            % @retval mu (m x 1) Coefficients of friction
            % @retval n (m x n) normal vector in joint coordinates, state vector length n
            % @retval D {2k}(m x n) friction cone basis in joint coordinates, for k directions
            % @retval dn (mn x n) dn/dq derivative
            % @retval dD {2k}(mn x n) dD/dq derivative
            
            compute_first_derivative = nargout > 8;
            compute_kinematics_gradients = nargout > 10;
            
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kin_options = struct('compute_gradients', compute_kinematics_gradients);
                kinsol = doKinematics(obj, kinsol, [], kin_options);
            end
            
            % Scott's implementation
            brick_size = 0.06;
            finger_contact_left = [0;0;.04];
            finger_contact_right = [0;  0.0400;  0.1225];%-0.0001
            
            box_pose = obj.forwardKin(kinsol,obj.cylinder_id,[0;0;0],1);
            left_finger_tip = obj.forwardKin(kinsol,obj.left_finger_id,finger_contact_left,0);
            right_finger_tip = obj.forwardKin(kinsol,obj.right_finger_id,finger_contact_right,0);
            
            xA = [finger_contact_right, finger_contact_left];
            idxA = [obj.right_finger_id; obj.left_finger_id];
            idxB = [obj.cylinder_id; obj.cylinder_id];
            mu = 1.0;
            
            R_box = rpy2rotmat(box_pose(4:6));
            right_normal = R_box*[1;0;0];
            left_normal = R_box*[-1;0;0];
            
            normal = [right_normal, left_normal];
            d = cell(1,2);
            d{1} = [R_box*[0;1;0],R_box*[0;1;0]];
            d{2} = [R_box*[0;0;1],R_box*[0;0;1]];
            
            right_finger_in_B = right_finger_tip-box_pose(1:3);
            left_finger_in_B = left_finger_tip-box_pose(1:3);
            
            phi_right = right_finger_in_B'*right_normal - brick_size/2.0;
            phi_left = left_finger_in_B'*left_normal - brick_size/2.0;
            xB = [right_finger_in_B - right_finger_in_B'*right_normal *right_normal + [brick_size/2;0;0], ...
                left_finger_in_B - left_finger_in_B'*left_normal * left_normal - [brick_size/2;0;0]];
            
            phi = [phi_right;phi_left];
            
            if compute_kinematics_gradients
                [n, D, dn, dD] = contactConstraintDerivatives(obj, normal, kinsol, idxA, idxB, xA, xB, d);
            elseif compute_first_derivative
                [n, D] = contactConstraintDerivatives(obj, normal, kinsol, idxA, idxB, xA, xB, d);
            end
        end
        
        function [phi,normal,d,xA,xB,idxA,idxB,mu] = worldContactConstraints(obj,kinsol)
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            
            %%
            normal_ground = repmat([0;0;1],1,4);
            d_ground{1} = repmat([1;0;0],1,4);
            d_ground{2} = repmat([0;-1;0],1,4);
            
            obj_pos = kinsol.q(9:11);
            obj_ori = kinsol.q(12:14);
            cylinder_radius = 0.03;
            cylinder_height_half = 0.09;
            R_world_to_obj = rpy2rotmat(obj_ori);
            
            xA_ground = [cylinder_radius, -cylinder_radius, 0, 0;
                0, 0, cylinder_radius, -cylinder_radius;
                -cylinder_height_half, -cylinder_height_half, -cylinder_height_half, -cylinder_height_half];
            
            xB_ground(:,1) = obj_pos + R_world_to_obj*[cylinder_radius;0;-cylinder_height_half];
            xB_ground(:,2) = obj_pos + R_world_to_obj*[-cylinder_radius;0;-cylinder_height_half];
            xB_ground(:,3) = obj_pos + R_world_to_obj*[0;cylinder_radius;-cylinder_height_half];
            xB_ground(:,4) = obj_pos + R_world_to_obj*[0;-cylinder_radius;-cylinder_height_half];
            
            phi_ground = xB_ground(3,:)';
            xB_ground(3,:) = zeros(1,4);
            
            %% pieces of code above compute the ground contact constraint manually
            % terrain_options=struct();
            % terrain_options.active_collision_options.terrain_only = true;
            % [phi_ground,normal_ground,d_ground,xA_ground,xB_ground] = obj.contactConstraints(kinsol.q,false,terrain_options.active_collision_options);
            %%
            
            n_ground_contact_point = 4;
            % note that, here A and B are inverted
            phi_ground = phi_ground(1:n_ground_contact_point);
            normal_ground = -normal_ground(:,1:n_ground_contact_point);
            d_ground{1} = -d_ground{1}(:,1:n_ground_contact_point);
            d_ground{2} = -d_ground{2}(:,1:n_ground_contact_point);
            xA_ground = xA_ground(:,1:n_ground_contact_point);
            xB_ground = xB_ground(:,1:n_ground_contact_point);
            xB_ground_tmp = xB_ground;
            xB_ground = xA_ground;
            xA_ground = xB_ground_tmp;
            
            % modified object and four contact points on each finger stick
            finger_contact_delta = 0.01;
            right_finger_y_shift = 0.04;
            finger_contact_left1 = [finger_contact_delta;0;.04];
            finger_contact_left2 = [finger_contact_delta;0;.02];
            finger_contact_left3 = [-finger_contact_delta;0;.04];
            finger_contact_left4 = [-finger_contact_delta;0;.02];
            finger_contact_right1 = [finger_contact_delta;right_finger_y_shift;0.1225];
            finger_contact_right2 = [finger_contact_delta;right_finger_y_shift;0.1025];
            finger_contact_right3 = [-finger_contact_delta;right_finger_y_shift;0.1225];
            finger_contact_right4 = [-finger_contact_delta;right_finger_y_shift;0.1025];
            
            b = kinsol.q(9:14);%obj.forwardKin(kinsol,obj.cylinder_id,[0;0;0],1);
            R_world_to_B = rpy2rotmat(b(4:6));
            
            %% D-H parameters of kuka arm
            theta = kinsol.q(1:7);
            % d1 = 0.4200;
            % d2 = 0.4000;
            % a = zeros(1,7);
            % alpha = [pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,0];
            % d = [0,0,d1,0,d2,0,0];
            % T = eye(4);
            % for i=1:7
            %     A(:,:,i) = [cos(theta(i)),-sin(theta(i))*cos(alpha(i)),sin(theta(i))*sin(alpha(i)),a(i)*cos(theta(i));
            %                 sin(theta(i)),cos(theta(i))*cos(alpha(i)),-cos(theta(i))*sin(alpha(i)),a(i)*sin(theta(i));
            %                 0,sin(alpha(i)),cos(alpha(i)),d(i);
            %                 0,0,0,1];
            %     T = T*A(:,:,i);
            % end
            
            % pieces of code above are embedded into the analytical solution below
            T = [sin(theta(7))*(sin(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) - cos(theta(7))*(sin(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3))))),   cos(theta(7))*(sin(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) + sin(theta(7))*(sin(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3))))), sin(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) - cos(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))), - (21*cos(theta(1))*sin(theta(2)))/50 - (2*sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))))/5 - (2*cos(theta(1))*cos(theta(4))*sin(theta(2)))/5;
                cos(theta(7))*(sin(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3))))) - sin(theta(7))*(sin(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))), - cos(theta(7))*(sin(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))) - sin(theta(7))*(sin(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3))))), cos(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) - sin(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))),   (2*sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))))/5 - (21*sin(theta(1))*sin(theta(2)))/50 - (2*cos(theta(4))*sin(theta(1))*sin(theta(2)))/5;
                sin(theta(7))*(sin(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) - cos(theta(5))*sin(theta(2))*sin(theta(3))) - cos(theta(7))*(cos(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) - sin(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4)))),                                                                                                                                                                                                                                                                 cos(theta(7))*(sin(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) - cos(theta(5))*sin(theta(2))*sin(theta(3))) + sin(theta(7))*(cos(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) - sin(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4)))),                                                                                                                                                            sin(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) + cos(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4))),                                                                        (21*cos(theta(2)))/50 + (2*cos(theta(2))*cos(theta(4)))/5 + (2*cos(theta(3))*sin(theta(2))*sin(theta(4)))/5;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                            1];
            
            R_DHbase_world = [-1,0,0;0,-1,0;0,0,1];% rotation from DH base coordinate to world coordinate (rotation around z axis)
            pos_ee = R_DHbase_world*[T(1,4);T(2,4);T(3,4)+0.36];%0.36 is the base vertical position w.r.t. world coordinate
            %R_ee = R_DHbase_world*T(1:3,1:3)*R_DHbase_world;
            %R_ee_fwdkin = rpy2rotmat(iiwa_link_7_init(4:6));
            
            fr1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
            fr2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
            fr3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
            fr4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
            
            fl1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1225];
            fl2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1025];
            fl3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1225];
            fl4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1025];
            %%
            
            fl1 = R_world_to_B'*fl1;
            fl2 = R_world_to_B'*fl2;
            fl3 = R_world_to_B'*fl3;
            fl4 = R_world_to_B'*fl4;
            
            fr1 = R_world_to_B'*fr1;
            fr2 = R_world_to_B'*fr2;
            fr3 = R_world_to_B'*fr3;
            fr4 = R_world_to_B'*fr4;
            
            b_local = R_world_to_B'*b(1:3);
            phi = [phi_ground; ...
                norm(fr1(1:2)-b_local(1:2))-cylinder_radius; norm(fr2(1:2)-b_local(1:2))-cylinder_radius; norm(fr3(1:2)-b_local(1:2))-cylinder_radius; ...
                norm(fr4(1:2)-b_local(1:2))-cylinder_radius; ...
                norm(fl1(1:2)-b_local(1:2))-cylinder_radius; norm(fl2(1:2)-b_local(1:2))-cylinder_radius; norm(fl3(1:2)-b_local(1:2))-cylinder_radius; ...
                norm(fl4(1:2)-b_local(1:2))-cylinder_radius];
            cylinder_normal = R_world_to_B'*[zeros(2,n_ground_contact_point);-ones(1,n_ground_contact_point)];%cylinder normal expressed in cylinder coordinate
            right_normal1 = [fr1(1:2) - b_local(1:2);0];
            right_normal1 = right_normal1./sqrt(right_normal1'*right_normal1);
            right_normal2 = [fr2(1:2) - b_local(1:2);0];
            right_normal2 = right_normal2./sqrt(right_normal2'*right_normal2);
            right_normal3 = [fr3(1:2) - b_local(1:2);0];
            right_normal3 = right_normal3./sqrt(right_normal3'*right_normal3);
            right_normal4 = [fr4(1:2) - b_local(1:2);0];
            right_normal4 = right_normal4./sqrt(right_normal4'*right_normal4);
            left_normal1 = [fl1(1:2) - b_local(1:2);0];
            left_normal1 = left_normal1./sqrt(left_normal1'*left_normal1);
            left_normal2 = [fl2(1:2) - b_local(1:2);0];
            left_normal2 = left_normal2./sqrt(left_normal2'*left_normal2);
            left_normal3 = [fl3(1:2) - b_local(1:2);0];
            left_normal3 = left_normal3./sqrt(left_normal3'*left_normal3);
            left_normal4 = [fl4(1:2) - b_local(1:2);0];
            left_normal4 = left_normal4./sqrt(left_normal4'*left_normal4);
            normal = [cylinder_normal, right_normal1, right_normal2, ...
                right_normal3, right_normal4, left_normal1, left_normal2, ...
                left_normal3, left_normal4];
            
            d = cell(1,2);
            Tr11 = [right_normal1(2), -right_normal1(1) ,0]'/sqrt(right_normal1(1)^2 + right_normal1(2)^2);
            Tr12 = [0,0,-1]';
            
            Tr21 = [right_normal2(2), -right_normal2(1) ,0]'/sqrt(right_normal2(1)^2 + right_normal2(2)^2);
            Tr22 = [0,0,-1]';
            
            Tr31 = [right_normal3(2), -right_normal3(1) ,0]'/sqrt(right_normal3(1)^2 + right_normal3(2)^2);
            Tr32 = [0,0,-1]';
            
            Tr41 = [right_normal4(2), -right_normal4(1) ,0]'/sqrt(right_normal4(1)^2 + right_normal4(2)^2);
            Tr42 = [0,0,-1]';
            
            Tl11 = [left_normal1(2), -left_normal1(1) ,0]'/sqrt(left_normal1(1)^2 + left_normal1(2)^2);
            Tl12 = [0,0,-1]';
            
            Tl21 = [left_normal2(2), -left_normal2(1) ,0]'/sqrt(left_normal2(1)^2 + left_normal2(2)^2);
            Tl22 = [0,0,-1]';
            
            Tl31 = [left_normal3(2), -left_normal3(1) ,0]'/sqrt(left_normal3(1)^2 + left_normal3(2)^2);
            Tl32 = [0,0,-1]';
            
            Tl41 = [left_normal4(2), -left_normal4(1) ,0]'/sqrt(left_normal4(1)^2 + left_normal4(2)^2);
            Tl42 = [0,0,-1]';
            
            d{1} = [R_world_to_B'*[-ones(1,n_ground_contact_point);zeros(2,n_ground_contact_point)],Tr11,Tr21,Tr31,Tr41,Tl11,Tl21,Tl31,Tl41];
            d{2} = [R_world_to_B'*[zeros(1,n_ground_contact_point);ones(1,n_ground_contact_point);zeros(1,n_ground_contact_point)],Tr12,Tr22,Tr32,Tr42,Tl12,Tl22,Tl32,Tl42];
            
            d{1} = R_world_to_B*d{1};
            d{2} = R_world_to_B*d{2};
            
            xA = [xA_ground, finger_contact_right1, finger_contact_right2, ...
                finger_contact_right3, finger_contact_right4, ...
                finger_contact_left1, finger_contact_left2, ...
                finger_contact_left3, finger_contact_left4];
            
            %define horizontal 2D position on the cylinder surface
            xB = cylinder_radius*normal;
            %define vertical heights of closest point on the cylinder w.r.t cylinder coordinate
            %cylinder_local = R_world_to_B'*[0;0;cylinder_height/2];
            %xB(3,1) = - cylinder_local(3);
            xB(:,1:n_ground_contact_point) = xB_ground;
            % x and y direction is not accurate, currently assume the central point. It should be a point on the edge
            xB(3,n_ground_contact_point+1) = fr1(3) - b_local(3);
            xB(3,n_ground_contact_point+2) = fr2(3) - b_local(3);
            xB(3,n_ground_contact_point+3) = fr3(3) - b_local(3);
            xB(3,n_ground_contact_point+4) = fr4(3) - b_local(3);
            
            xB(3,n_ground_contact_point+5) = fl1(3) - b_local(3);
            xB(3,n_ground_contact_point+6) = fl2(3) - b_local(3);
            xB(3,n_ground_contact_point+7) = fl3(3) - b_local(3);
            xB(3,n_ground_contact_point+8) = fl4(3) - b_local(3);
            
            normal = R_world_to_B*normal;
            
            % todo: when the object is not betwen two finger tips.
            % todo: add caps at the end of cylinders
            nC = 8+n_ground_contact_point;
            idxA = [ones(n_ground_contact_point,1); obj.right_finger_id; obj.right_finger_id; ...
                obj.right_finger_id; obj.right_finger_id; ...
                obj.left_finger_id; obj.left_finger_id; ...
                obj.left_finger_id; obj.left_finger_id];
            idxB = obj.cylinder_id*ones(nC,1);
            
            if strcmp(obj.uncertainty_source, 'friction_coeff') || strcmp(obj.uncertainty_source, 'friction_coeff+object_initial_position')
                if isempty(obj.uncertain_mu)
                    mu = obj.uncertain_mu_mean*ones(nC,1);
                else
                    mu = obj.uncertain_mu*ones(nC,1);
                end
            else
                mu = 4.0*ones(nC,1);
            end
        end
        
        function [normal,d,mu] = worldContactConstraints_v2(obj,kinsol)
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            
            obj_pos = kinsol.q(9:11);
            obj_ori = kinsol.q(12:14);
            cylinder_radius = 0.03;
            cylinder_height_half = 0.09;
            R_world_to_obj = rpy2rotmat(obj_ori);
            
            % xA_ground = [cylinder_radius, -cylinder_radius, 0, 0;
            %     0, 0, cylinder_radius, -cylinder_radius;
            %     -cylinder_height_half, -cylinder_height_half, -cylinder_height_half, -cylinder_height_half];
            
            xB_ground(:,1) = obj_pos + R_world_to_obj*[cylinder_radius;0;-cylinder_height_half];
            xB_ground(:,2) = obj_pos + R_world_to_obj*[-cylinder_radius;0;-cylinder_height_half];
            xB_ground(:,3) = obj_pos + R_world_to_obj*[0;cylinder_radius;-cylinder_height_half];
            xB_ground(:,4) = obj_pos + R_world_to_obj*[0;-cylinder_radius;-cylinder_height_half];
            
            % phi_ground = xB_ground(3,:)';
            % xB_ground(3,:) = zeros(1,4);
            
            %% pieces of code above compute the ground contact constraint manually
            % terrain_options=struct();
            % terrain_options.active_collision_options.terrain_only = true;
            % [phi_ground,normal_ground,d_ground,xA_ground,xB_ground] = obj.contactConstraints(kinsol.q,false,terrain_options.active_collision_options);
            %%
            
            n_ground_contact_point = 4;
            % note that, here A and B are inverted
            % phi_ground = phi_ground(1:n_ground_contact_point);
            % normal_ground = -normal_ground(:,1:n_ground_contact_point);
            % d_ground{1} = -d_ground{1}(:,1:n_ground_contact_point);
            % d_ground{2} = -d_ground{2}(:,1:n_ground_contact_point);
            % xA_ground = xA_ground(:,1:n_ground_contact_point);
            % xB_ground = xB_ground(:,1:n_ground_contact_point);
            % xB_ground_tmp = xB_ground;
            % xB_ground = xA_ground;
            % xA_ground = xB_ground_tmp;
            %
            % modified object and four contact points on each finger stick
            finger_contact_delta = 0.01;
            right_finger_y_shift = 0.04;
            
            b = kinsol.q(9:14);%obj.forwardKin(kinsol,obj.cylinder_id,[0;0;0],1);
            R_world_to_B = rpy2rotmat(b(4:6));
            
            %% D-H parameters of kuka arm
            theta = kinsol.q(1:7);
            % d1 = 0.4200;
            % d2 = 0.4000;
            % a = zeros(1,7);
            % alpha = [pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,0];
            % d = [0,0,d1,0,d2,0,0];
            % T = eye(4);
            % for i=1:7
            %     A(:,:,i) = [cos(theta(i)),-sin(theta(i))*cos(alpha(i)),sin(theta(i))*sin(alpha(i)),a(i)*cos(theta(i));
            %                 sin(theta(i)),cos(theta(i))*cos(alpha(i)),-cos(theta(i))*sin(alpha(i)),a(i)*sin(theta(i));
            %                 0,sin(alpha(i)),cos(alpha(i)),d(i);
            %                 0,0,0,1];
            %     T = T*A(:,:,i);
            % end
            
            % pieces of code above are embedded into the analytical solution below
            T = [sin(theta(7))*(sin(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) - cos(theta(7))*(sin(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3))))),   cos(theta(7))*(sin(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) + sin(theta(7))*(sin(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3))))), sin(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) - cos(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))), - (21*cos(theta(1))*sin(theta(2)))/50 - (2*sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))))/5 - (2*cos(theta(1))*cos(theta(4))*sin(theta(2)))/5;
                cos(theta(7))*(sin(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3))))) - sin(theta(7))*(sin(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))), - cos(theta(7))*(sin(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))) - sin(theta(7))*(sin(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3))))), cos(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) - sin(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))),   (2*sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))))/5 - (21*sin(theta(1))*sin(theta(2)))/50 - (2*cos(theta(4))*sin(theta(1))*sin(theta(2)))/5;
                sin(theta(7))*(sin(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) - cos(theta(5))*sin(theta(2))*sin(theta(3))) - cos(theta(7))*(cos(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) - sin(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4)))),                                                                                                                                                                                                                                                                 cos(theta(7))*(sin(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) - cos(theta(5))*sin(theta(2))*sin(theta(3))) + sin(theta(7))*(cos(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) - sin(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4)))),                                                                                                                                                            sin(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) + cos(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4))),                                                                        (21*cos(theta(2)))/50 + (2*cos(theta(2))*cos(theta(4)))/5 + (2*cos(theta(3))*sin(theta(2))*sin(theta(4)))/5;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                            1];
            
            R_DHbase_world = [-1,0,0;0,-1,0;0,0,1];% rotation from DH base coordinate to world coordinate (rotation around z axis)
            pos_ee = R_DHbase_world*[T(1,4);T(2,4);T(3,4)+0.36];%0.36 is the base vertical position w.r.t. world coordinate
            %R_ee = R_DHbase_world*T(1:3,1:3)*R_DHbase_world;
            %R_ee_fwdkin = rpy2rotmat(iiwa_link_7_init(4:6));
            
            fr1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
            fr2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
            fr3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
            fr4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
            
            fl1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1225];
            fl2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1025];
            fl3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1225];
            fl4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1025];
            
            %%
            fl1 = R_world_to_B'*fl1;
            fl2 = R_world_to_B'*fl2;
            fl3 = R_world_to_B'*fl3;
            fl4 = R_world_to_B'*fl4;
            
            fr1 = R_world_to_B'*fr1;
            fr2 = R_world_to_B'*fr2;
            fr3 = R_world_to_B'*fr3;
            fr4 = R_world_to_B'*fr4;
            
            b_local = R_world_to_B'*b(1:3);
            cylinder_normal = R_world_to_B'*[zeros(2,n_ground_contact_point);-ones(1,n_ground_contact_point)];%cylinder normal expressed in cylinder coordinate
            right_normal1 = [fr1(1:2) - b_local(1:2);0];
            right_normal1 = right_normal1./sqrt(right_normal1'*right_normal1);
            right_normal2 = [fr2(1:2) - b_local(1:2);0];
            right_normal2 = right_normal2./sqrt(right_normal2'*right_normal2);
            right_normal3 = [fr3(1:2) - b_local(1:2);0];
            right_normal3 = right_normal3./sqrt(right_normal3'*right_normal3);
            right_normal4 = [fr4(1:2) - b_local(1:2);0];
            right_normal4 = right_normal4./sqrt(right_normal4'*right_normal4);
            left_normal1 = [fl1(1:2) - b_local(1:2);0];
            left_normal1 = left_normal1./sqrt(left_normal1'*left_normal1);
            left_normal2 = [fl2(1:2) - b_local(1:2);0];
            left_normal2 = left_normal2./sqrt(left_normal2'*left_normal2);
            left_normal3 = [fl3(1:2) - b_local(1:2);0];
            left_normal3 = left_normal3./sqrt(left_normal3'*left_normal3);
            left_normal4 = [fl4(1:2) - b_local(1:2);0];
            left_normal4 = left_normal4./sqrt(left_normal4'*left_normal4);
            normal = [cylinder_normal, right_normal1, right_normal2, ...
                right_normal3, right_normal4, left_normal1, left_normal2, ...
                left_normal3, left_normal4];
            
            d = cell(1,2);
            % Tr11 = cross(right_normal1,[0;0;1]);
            % Tr11 = Tr11/norm(Tr11);
            % Tr12 = cross(right_normal1,Tr11);
            %
            % Tr21 = cross(right_normal2,[0;0;1]);
            % Tr21 = Tr21/norm(Tr21);
            % Tr22 = cross(right_normal2,Tr21);
            %
            % Tr31 = cross(right_normal3,[0;0;1]);
            % Tr31 = Tr31/norm(Tr31);
            % Tr32 = cross(right_normal3,Tr31);
            %
            % Tr41 = cross(right_normal4,[0;0;1]);
            % Tr41 = Tr41/norm(Tr41);
            % Tr42 = cross(right_normal4,Tr41);
            %
            % Tl11 = cross(left_normal1,[0;0;1]);
            % Tl11 = Tl11/norm(Tl11);
            % Tl12 = cross(left_normal1,Tl11);
            %
            % Tl21 = cross(left_normal2,[0;0;1]);
            % Tl21 = Tl21/norm(Tl21);
            % Tl22 = cross(left_normal2,Tl21);
            %
            % Tl31 = cross(left_normal3,[0;0;1]);
            % Tl31 = Tl31/norm(Tl31);
            % Tl32 = cross(left_normal3,Tl31);
            %
            % Tl41 = cross(left_normal4,[0;0;1]);
            % Tl41 = Tl41/norm(Tl41);
            % Tl42 = cross(left_normal4,Tl41);
            
            Tr11 = [right_normal1(2), -right_normal1(1) ,0]'/sqrt(right_normal1(1)^2 + right_normal1(2)^2);
            Tr12 = [0,0,-1]';
            
            Tr21 = [right_normal2(2), -right_normal2(1) ,0]'/sqrt(right_normal2(1)^2 + right_normal2(2)^2);
            Tr22 = [0,0,-1]';
            
            Tr31 = [right_normal3(2), -right_normal3(1) ,0]'/sqrt(right_normal3(1)^2 + right_normal3(2)^2);
            Tr32 = [0,0,-1]';
            
            Tr41 = [right_normal4(2), -right_normal4(1) ,0]'/sqrt(right_normal4(1)^2 + right_normal4(2)^2);
            Tr42 = [0,0,-1]';
            
            Tl11 = [left_normal1(2), -left_normal1(1) ,0]'/sqrt(left_normal1(1)^2 + left_normal1(2)^2);
            Tl12 = [0,0,-1]';
            
            Tl21 = [left_normal2(2), -left_normal2(1) ,0]'/sqrt(left_normal2(1)^2 + left_normal2(2)^2);
            Tl22 = [0,0,-1]';
            
            Tl31 = [left_normal3(2), -left_normal3(1) ,0]'/sqrt(left_normal3(1)^2 + left_normal3(2)^2);
            Tl32 = [0,0,-1]';
            
            Tl41 = [left_normal4(2), -left_normal4(1) ,0]'/sqrt(left_normal4(1)^2 + left_normal4(2)^2);
            Tl42 = [0,0,-1]';
            
            d{1} = [R_world_to_B'*[-ones(1,n_ground_contact_point);zeros(2,n_ground_contact_point)],Tr11,Tr21,Tr31,Tr41,Tl11,Tl21,Tl31,Tl41];
            d{2} = [R_world_to_B'*[zeros(1,n_ground_contact_point);ones(1,n_ground_contact_point);zeros(1,n_ground_contact_point)],Tr12,Tr22,Tr32,Tr42,Tl12,Tl22,Tl32,Tl42];
            
            d{1} = R_world_to_B*d{1};
            d{2} = R_world_to_B*d{2};
            
            % xA = [xA_ground, finger_contact_right1, finger_contact_right2, ...
            %     finger_contact_right3, finger_contact_right4, ...
            %     finger_contact_left1, finger_contact_left2, ...
            %     finger_contact_left3, finger_contact_left4];
            %
            % %define horizontal 2D position on the cylinder surface
            % xB = cylinder_radius*normal;
            % %define vertical heights of closest point on the cylinder w.r.t cylinder coordinate
            % %cylinder_local = R_world_to_B'*[0;0;cylinder_height/2];
            % %xB(3,1) = - cylinder_local(3);
            % xB(:,1:n_ground_contact_point) = xB_ground;
            % % x and y direction is not accurate, currently assume the central point. It should be a point on the edge
            % xB(3,n_ground_contact_point+1) = fr1(3) - b_local(3);
            % xB(3,n_ground_contact_point+2) = fr2(3) - b_local(3);
            % xB(3,n_ground_contact_point+3) = fr3(3) - b_local(3);
            % xB(3,n_ground_contact_point+4) = fr4(3) - b_local(3);
            %
            % xB(3,n_ground_contact_point+5) = fl1(3) - b_local(3);
            % xB(3,n_ground_contact_point+6) = fl2(3) - b_local(3);
            % xB(3,n_ground_contact_point+7) = fl3(3) - b_local(3);
            % xB(3,n_ground_contact_point+8) = fl4(3) - b_local(3);
            
            normal = R_world_to_B*normal;
            
            % todo: when the object is not betwen two finger tips.
            % todo: add caps at the end of cylinders
            nC = 8+n_ground_contact_point;
            
            if strcmp(obj.uncertainty_source, 'friction_coeff') || strcmp(obj.uncertainty_source, 'friction_coeff+object_initial_position')
                mu = obj.uncertain_mu*ones(nC,1);
            else
                mu = 1.0*ones(nC,1);
            end
        end
        
        function [xB] = computexB(obj,kinsol)
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            
            obj_pos = kinsol.q(9:11);
            obj_ori = kinsol.q(12:14);
            cylinder_radius = 0.03;
            cylinder_height_half = 0.09;
            R_world_to_obj = rpy2rotmat(obj_ori);
            
            xA_ground = [cylinder_radius, -cylinder_radius, 0, 0;
                0, 0, cylinder_radius, -cylinder_radius;
                -cylinder_height_half, -cylinder_height_half, -cylinder_height_half, -cylinder_height_half];
            
            %% pieces of code above compute the ground contact constraint manually
            % terrain_options=struct();
            % terrain_options.active_collision_options.terrain_only = true;
            % [phi_ground,normal_ground,d_ground,xA_ground,xB_ground] = obj.contactConstraints(kinsol.q,false,terrain_options.active_collision_options);
            %%
            
            n_ground_contact_point = 4;
            % note that, here A and B are inverted
            xA_ground = xA_ground(:,1:n_ground_contact_point);
            xB_ground = xA_ground;
            
            % modified object and four contact points on each finger stick
            finger_contact_delta = 0.01;
            right_finger_y_shift = 0.04;
            
            b = kinsol.q(9:14);%obj.forwardKin(kinsol,obj.cylinder_id,[0;0;0],1);
            R_world_to_B = rpy2rotmat(b(4:6));
            
            %% D-H parameters of kuka arm
            theta = kinsol.q(1:7);
            % d1 = 0.4200;
            % d2 = 0.4000;
            % a = zeros(1,7);
            % alpha = [pi/2,-pi/2,-pi/2,pi/2,pi/2,-pi/2,0];
            % d = [0,0,d1,0,d2,0,0];
            % T = eye(4);
            % for i=1:7
            %     A(:,:,i) = [cos(theta(i)),-sin(theta(i))*cos(alpha(i)),sin(theta(i))*sin(alpha(i)),a(i)*cos(theta(i));
            %                 sin(theta(i)),cos(theta(i))*cos(alpha(i)),-cos(theta(i))*sin(alpha(i)),a(i)*sin(theta(i));
            %                 0,sin(alpha(i)),cos(alpha(i)),d(i);
            %                 0,0,0,1];
            %     T = T*A(:,:,i);
            % end
            
            % pieces of code above are embedded into the analytical solution below
            T = [sin(theta(7))*(sin(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) - cos(theta(7))*(sin(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3))))),   cos(theta(7))*(sin(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) + sin(theta(7))*(sin(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3))))), sin(theta(6))*(cos(theta(5))*(cos(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) - cos(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(3))*sin(theta(1)) + cos(theta(1))*cos(theta(2))*sin(theta(3)))) - cos(theta(6))*(sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))) + cos(theta(1))*cos(theta(4))*sin(theta(2))), - (21*cos(theta(1))*sin(theta(2)))/50 - (2*sin(theta(4))*(sin(theta(1))*sin(theta(3)) - cos(theta(1))*cos(theta(2))*cos(theta(3))))/5 - (2*cos(theta(1))*cos(theta(4))*sin(theta(2)))/5;
                cos(theta(7))*(sin(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3))))) - sin(theta(7))*(sin(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))), - cos(theta(7))*(sin(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) - cos(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))) - sin(theta(7))*(sin(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) + cos(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3))))), cos(theta(6))*(sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) - cos(theta(4))*sin(theta(1))*sin(theta(2))) - sin(theta(6))*(cos(theta(5))*(cos(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))) + sin(theta(1))*sin(theta(2))*sin(theta(4))) + sin(theta(5))*(cos(theta(1))*cos(theta(3)) - cos(theta(2))*sin(theta(1))*sin(theta(3)))),   (2*sin(theta(4))*(cos(theta(1))*sin(theta(3)) + cos(theta(2))*cos(theta(3))*sin(theta(1))))/5 - (21*sin(theta(1))*sin(theta(2)))/50 - (2*cos(theta(4))*sin(theta(1))*sin(theta(2)))/5;
                sin(theta(7))*(sin(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) - cos(theta(5))*sin(theta(2))*sin(theta(3))) - cos(theta(7))*(cos(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) - sin(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4)))),                                                                                                                                                                                                                                                                 cos(theta(7))*(sin(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) - cos(theta(5))*sin(theta(2))*sin(theta(3))) + sin(theta(7))*(cos(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) - sin(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4)))),                                                                                                                                                            sin(theta(6))*(cos(theta(5))*(cos(theta(2))*sin(theta(4)) - cos(theta(3))*cos(theta(4))*sin(theta(2))) + sin(theta(2))*sin(theta(3))*sin(theta(5))) + cos(theta(6))*(cos(theta(2))*cos(theta(4)) + cos(theta(3))*sin(theta(2))*sin(theta(4))),                                                                        (21*cos(theta(2)))/50 + (2*cos(theta(2))*cos(theta(4)))/5 + (2*cos(theta(3))*sin(theta(2))*sin(theta(4)))/5;
                0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                                                                                                                                                                                                                                        0,                                                                                                                                                                            1];
            
            R_DHbase_world = [-1,0,0;0,-1,0;0,0,1];% rotation from DH base coordinate to world coordinate (rotation around z axis)
            pos_ee = R_DHbase_world*[T(1,4);T(2,4);T(3,4)+0.36];%0.36 is the base vertical position w.r.t. world coordinate
            %R_ee = R_DHbase_world*T(1:3,1:3)*R_DHbase_world;
            %R_ee_fwdkin = rpy2rotmat(iiwa_link_7_init(4:6));
            
            fr1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
            fr2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
            fr3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1225];
            fr4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;-right_finger_y_shift;0.081+0.1025];
            
            fl1 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1225];
            fl2 = pos_ee + R_DHbase_world*T(1:3,1:3)*[-finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1025];
            fl3 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1225];
            fl4 = pos_ee + R_DHbase_world*T(1:3,1:3)*[finger_contact_delta;(kinsol.q(8)-0.04);0.081+0.1025];
            %%
            
            fl1 = R_world_to_B'*fl1;
            fl2 = R_world_to_B'*fl2;
            fl3 = R_world_to_B'*fl3;
            fl4 = R_world_to_B'*fl4;
            
            fr1 = R_world_to_B'*fr1;
            fr2 = R_world_to_B'*fr2;
            fr3 = R_world_to_B'*fr3;
            fr4 = R_world_to_B'*fr4;
            
            b_local = R_world_to_B'*b(1:3);
            %cylinder_normal = R_world_to_B'*[zeros(2,n_ground_contact_point);-ones(1,n_ground_contact_point)];%cylinder normal expressed in cylinder coordinate
            right_normal1 = [fr1(1:2) - b_local(1:2);0];
            right_normal1 = right_normal1./sqrt(right_normal1'*right_normal1);
            right_normal2 = [fr2(1:2) - b_local(1:2);0];
            right_normal2 = right_normal2./sqrt(right_normal2'*right_normal2);
            right_normal3 = [fr3(1:2) - b_local(1:2);0];
            right_normal3 = right_normal3./sqrt(right_normal3'*right_normal3);
            right_normal4 = [fr4(1:2) - b_local(1:2);0];
            right_normal4 = right_normal4./sqrt(right_normal4'*right_normal4);
            left_normal1 = [fl1(1:2) - b_local(1:2);0];
            left_normal1 = left_normal1./sqrt(left_normal1'*left_normal1);
            left_normal2 = [fl2(1:2) - b_local(1:2);0];
            left_normal2 = left_normal2./sqrt(left_normal2'*left_normal2);
            left_normal3 = [fl3(1:2) - b_local(1:2);0];
            left_normal3 = left_normal3./sqrt(left_normal3'*left_normal3);
            left_normal4 = [fl4(1:2) - b_local(1:2);0];
            left_normal4 = left_normal4./sqrt(left_normal4'*left_normal4);
            normal = [zeros(3,n_ground_contact_point), right_normal1, right_normal2, ...
                right_normal3, right_normal4, left_normal1, left_normal2, ...
                left_normal3, left_normal4];
            
            %define horizontal 2D position on the cylinder surface
            xB = cylinder_radius*normal;
            %define vertical heights of closest point on the cylinder w.r.t cylinder coordinate
            %cylinder_local = R_world_to_B'*[0;0;cylinder_height/2];
            %xB(3,1) = - cylinder_local(3);
            xB(:,1:n_ground_contact_point) = xB_ground;
            % x and y direction is not accurate, currently assume the central point. It should be a point on the edge
            xB(3,n_ground_contact_point+1) = fr1(3) - b_local(3);
            xB(3,n_ground_contact_point+2) = fr2(3) - b_local(3);
            xB(3,n_ground_contact_point+3) = fr3(3) - b_local(3);
            xB(3,n_ground_contact_point+4) = fr4(3) - b_local(3);
            
            xB(3,n_ground_contact_point+5) = fl1(3) - b_local(3);
            xB(3,n_ground_contact_point+6) = fl2(3) - b_local(3);
            xB(3,n_ground_contact_point+7) = fl3(3) - b_local(3);
            xB(3,n_ground_contact_point+8) = fl4(3) - b_local(3);
        end
        
        function [n,D] = jointContactJacobians(obj,kinsol)
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            
            [~,normal,d,xA,xB,idxA,idxB,mu] = worldContactConstraints(obj,kinsol);
            
            [n, D] = contactConstraintDerivatives(obj, normal, kinsol, idxA, idxB, xA, xB, d);
        end
        
        function [phi,normal,d,xA,xB,idxA,idxB,mu,n,D,dn,dD] = contactConstraints_manual(obj,kinsol,allow_multiple_contacts,active_collision_options)
            % @retval phi (m x 1) Vector of gap function values (typically contact distance), for m possible contacts
            % @retval normal (3 x m) Contact normal vector in world coordinates, points from B to A
            % @retval d {k} (3 x m) Contact friction basis vectors in world coordinates, points from B to A
            % @retval xA (3 x m) The closest point on body A to contact with body B, relative to body A origin and in body A frame
            % @retval xB (3 x m) The closest point on body B to contact with body A, relative to body B origin and in body B frame
            % @retval idxA (m x 1) The index of body A. 0 is the special case for the environment/terrain
            % @retval idxB (m x 1) The index of body B. 0 is the special case for the environment/terrain
            % @retval mu (m x 1) Coefficients of friction
            % @retval n (m x n) normal vector in joint coordinates, state vector length n
            % @retval D {2k}(m x n) friction cone basis in joint coordinates, for k directions
            % @retval dn (mn x n) dn/dq derivative
            % @retval dD {2k}(mn x n) dD/dq derivative
            
            compute_first_derivative = nargout > 8;
            compute_kinematics_gradients = nargout > 10;
            
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            
            [phi,normal,d,xA,xB,idxA,idxB,mu] = worldContactConstraints(obj,kinsol);
            
            if compute_first_derivative
                if ~isstruct(kinsol)
                    % treat input as contactPositions(obj,q)
                    kinsol = doKinematics(obj, kinsol, []);
                end
                [n, D] = contactConstraintDerivatives(obj, normal, kinsol, idxA, idxB, xA, xB, d);
            end
            
            if compute_kinematics_gradients
                %finite diff since we don't trust this
                %[~, ~, dn, dD] = contactConstraintDerivatives(obj, normal, kinsol, idxA, idxB, xA, xB, d);
                
                dn = zeros(numel(n),size(n,2));
                dD = cell(1,4);
                dD{1} = dn;
                dD{2} = dn;
                dD{3} = dn;
                dD{4} = dn;
                dq = diag(sqrt(eps(kinsol.q)))*1e-2;
                for k = 1:14
                    [np,Dp] = jointContactJacobians(obj,kinsol.q+dq(:,k));
                    [nm,Dm] = jointContactJacobians(obj,kinsol.q-dq(:,k));
                    dn(:,k) = vec(np-nm)/(2*dq(k,k));%1e-6;% [Ye: why 1e-6]
                    dD{1}(:,k) = vec(Dp{1}-Dm{1})/(2*dq(k,k));
                    dD{2}(:,k) = vec(Dp{2}-Dm{2})/(2*dq(k,k));
                    dD{3}(:,k) = vec(Dp{3}-Dm{3})/(2*dq(k,k));
                    dD{4}(:,k) = vec(Dp{4}-Dm{4})/(2*dq(k,k));
                end
            else
                dn = [];
                dD = [];
            end
        end
        
        function [normal,d,mu] = contactConstraints_manual_v2(obj,kinsol)
            % @retval normal (3 x m) Contact normal vector in world coordinates, points from B to A
            % @retval d {k} (3 x m) Contact friction basis vectors in world coordinates, points from B to A
            % @retval mu (m x 1) Coefficients of friction
            
            if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            [normal,d,mu] = worldContactConstraints_v2(obj,kinsol);
        end
    end
end