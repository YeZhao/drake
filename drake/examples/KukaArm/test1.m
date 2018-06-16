if ~isstruct(kinsol)
                % treat input as contactPositions(obj,q)
                kinsol = doKinematics(obj, kinsol, []);
            end
            
            obj_pos = kinsol.q(9:11);
            obj_ori = kinsol.q(12:14);
            cylinder_radius = 0.03;
            cylinder_height_half = 0.09;
            R_world_to_obj = rpy2rotmat(obj_ori);
            
            %             xA_ground = [cylinder_radius, -cylinder_radius, 0, 0;
            %                              0, 0, cylinder_radius, -cylinder_radius;
            %                              -cylinder_height_half, -cylinder_height_half, -cylinder_height_half, -cylinder_height_half];
            
            xB_ground(:,1) = obj_pos + R_world_to_obj*[cylinder_radius;0;-cylinder_height_half];
            xB_ground(:,2) = obj_pos + R_world_to_obj*[-cylinder_radius;0;-cylinder_height_half];
            xB_ground(:,3) = obj_pos + R_world_to_obj*[0;cylinder_radius;-cylinder_height_half];
            xB_ground(:,4) = obj_pos + R_world_to_obj*[0;-cylinder_radius;-cylinder_height_half];
            
            %             phi_ground = xB_ground(3,:)';
            %             xB_ground(3,:) = zeros(1,4);
            
            %% pieces of code above compute the ground contact constraint manually
            % terrain_options=struct();
            % terrain_options.active_collision_options.terrain_only = true;
            % [phi_ground,normal_ground,d_ground,xA_ground,xB_ground] = obj.contactConstraints(kinsol.q,false,terrain_options.active_collision_options);
            %%
            
            n_ground_contact_point = 4;
            % note that, here A and B are inverted
            %             phi_ground = phi_ground(1:n_ground_contact_point);
            %             normal_ground = -normal_ground(:,1:n_ground_contact_point);
            %             d_ground{1} = -d_ground{1}(:,1:n_ground_contact_point);
            %             d_ground{2} = -d_ground{2}(:,1:n_ground_contact_point);
            %             xA_ground = xA_ground(:,1:n_ground_contact_point);
            %             xB_ground = xB_ground(:,1:n_ground_contact_point);
            %             xB_ground_tmp = xB_ground;
            %             xB_ground = xA_ground;
            %             xA_ground = xB_ground_tmp;
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
                        Tr11 = cross(right_normal1,[0;0;1]);
                        Tr11 = Tr11/norm(Tr11);
                        Tr12 = cross(right_normal1,Tr11);
            
                        Tr21 = cross(right_normal2,[0;0;1]);
                        Tr21 = Tr21/norm(Tr21);
                        Tr22 = cross(right_normal2,Tr21);
            
                        Tr31 = cross(right_normal3,[0;0;1]);
                        Tr31 = Tr31/norm(Tr31);
                        Tr32 = cross(right_normal3,Tr31);
            
                        Tr41 = cross(right_normal4,[0;0;1]);
                        Tr41 = Tr41/norm(Tr41);
                        Tr42 = cross(right_normal4,Tr41);
            
                        Tl11 = cross(left_normal1,[0;0;1]);
                        Tl11 = Tl11/norm(Tl11);
                        Tl12 = cross(left_normal1,Tl11);
            
                        Tl21 = cross(left_normal2,[0;0;1]);
                        Tl21 = Tl21/norm(Tl21);
                        Tl22 = cross(left_normal2,Tl21);
            
                        Tl31 = cross(left_normal3,[0;0;1]);
                        Tl31 = Tl31/norm(Tl31);
                        Tl32 = cross(left_normal3,Tl31);
            
                        Tl41 = cross(left_normal4,[0;0;1]);
                        Tl41 = Tl41/norm(Tl41);
                        Tl42 = cross(left_normal4,Tl41);
            
                        if (Tr12(3) ~= -1)
                            keyboard
                        end
                        
                        if (Tr22(3) ~= -1)
                            keyboard
                        end
                        
                        if (Tr32(3) ~= -1)
                            keyboard
                        end
                        
                        if (Tr42(3) ~= -1)
                            keyboard
                        end
                            
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
            
            %             if ((abs(sum(sum(d_test{1} - d{1}))) > 1e-7) | (abs(sum(sum(d_test{2} - d{2}))) > 1e-7))
            %                 keyboard
            %             end
            
            %             xA = [xA_ground, finger_contact_right1, finger_contact_right2, ...
            %                 finger_contact_right3, finger_contact_right4, ...
            %                 finger_contact_left1, finger_contact_left2, ...
            %                 finger_contact_left3, finger_contact_left4];
            %
            %             %define horizontal 2D position on the cylinder surface
            %             xB = cylinder_radius*normal;
            %             %define vertical heights of closest point on the cylinder w.r.t cylinder coordinate
            %             %cylinder_local = R_world_to_B'*[0;0;cylinder_height/2];
            %             %xB(3,1) = - cylinder_local(3);
            %             xB(:,1:n_ground_contact_point) = xB_ground;
            %             % x and y direction is not accurate, currently assume the central point. It should be a point on the edge
            %             xB(3,n_ground_contact_point+1) = fr1(3) - b_local(3);
            %             xB(3,n_ground_contact_point+2) = fr2(3) - b_local(3);
            %             xB(3,n_ground_contact_point+3) = fr3(3) - b_local(3);
            %             xB(3,n_ground_contact_point+4) = fr4(3) - b_local(3);
            %
            %             xB(3,n_ground_contact_point+5) = fl1(3) - b_local(3);
            %             xB(3,n_ground_contact_point+6) = fl2(3) - b_local(3);
            %             xB(3,n_ground_contact_point+7) = fl3(3) - b_local(3);
            %             xB(3,n_ground_contact_point+8) = fl4(3) - b_local(3);
            
            normal = R_world_to_B*normal;
            
            % todo: when the object is not betwen two finger tips.
            % todo: add caps at the end of cylinders
            nC = 8+n_ground_contact_point;
            
            if strcmp(obj.uncertainty_source, 'friction_coeff') || strcmp(obj.uncertainty_source, 'friction_coeff+object_initial_position')
                mu = obj.uncertain_mu*ones(nC,1);
            else
                mu = 1.0*ones(nC,1);
            end