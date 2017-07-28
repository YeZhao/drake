classdef RigidBodyStepTerrainMultipleSteps < RigidBodyTerrain
    %  This class provides an implementation of RigidBodyTerrain with z=0
    %  everywhere, except on a series of boxes of center, height specified
    %  on creation.
    
    % Boxes specified as an Nx5 matrix with columns:
    %  X center, Y center, Width (x), Length (y), Height (z)
    methods
        function obj = RigidBodyStepTerrainMultipleSteps(boxes)
            %height = 0.05*randn(1,15);
            % temporarily fix the rough terrain height
            height = [0.0052, -0.0001, -0.0116, -0.0001, -0.0069, -0.0067, ... 
                      0.0086, 0.0011, 0.0040, 0.0088, 0.0018, 0.0055, 0.0068, 0.0117, 0.0048];
            height = zeros(1,15);
            
            obj.height = height;
            if nargin < 2 
                boxes(1,:) = [0.4, 0.0, 0.4, 2, 0.5];
                for i=2:length(height)
                    boxes(i,:) = [0.4*i, 0.0, 0.4, 2, 0.5];
                end
            end
            obj.centers_x = boxes(:, 1);
            obj.centers_y = boxes(:, 2);
            obj.boxs_x = boxes(:, 3);
            obj.boxs_y = boxes(:, 4);
            obj.boxs_z = boxes(:, 5);
            obj.num_boxes = size(boxes, 1);
            obj.geom = constructGeometry(obj);
        end
        
        function [z,normal] = getHeight(obj,xy)
            [~, n] = size(xy);
            x = repmat(xy(1, :), [obj.num_boxes, 1]);
            y = repmat(xy(2, :), [obj.num_boxes, 1]);
            
            x_above_boxes = x >= repmat(obj.centers_x - obj.boxs_x/2, [1, size(xy, 2)]);
            x_below_boxes = x <= repmat(obj.centers_x+ obj.boxs_x/2, [1, size(xy, 2)]);
            x_in_boxes = x_above_boxes & x_below_boxes;
            
            y_above_boxes = y >= repmat(obj.centers_y - obj.boxs_y/2, [1, size(xy, 2)]);
            y_below_boxes = y <= repmat(obj.centers_y+ obj.boxs_y/2, [1, size(xy, 2)]);
            y_in_boxes = y_above_boxes & y_below_boxes;
            
            in_boxes = x_in_boxes & y_in_boxes;            
            % Give the height of the heighest box that each point falls in
            adjusted_heights = repmat(obj.height', [1, size(xy, 2)]) .* in_boxes; % [Ye: modified the terrain height]
            z = max(adjusted_heights, [], 1);
            if ~any(z > 0)% detect terrain height when it has a negative value
                z = min(adjusted_heights, [], 1);
            end
            
            % Normals always straight up
            normal = repmat([0;0;1],1,n);
        end
        
        function geom = getCollisionGeometry(obj)
            geom = obj.geom;
        end
        
        function geom = getVisualGeometry(obj)
            geom = obj.geom;
        end
        
        function geom = constructGeometry(obj)
            ground_width = 1.5;
            ground_depth = 0.5;
            geom_ground = RigidBodyBox([ground_width;ground_width;ground_depth]);
            geom_ground.T(1,4) = 0 - ground_width/2 + 0.2;%assume 0.2 is slightly smaller than one step length
            geom_ground.T(3,4) = 0 - ground_depth/2;
            geom_ground.c = hex2dec({'ee','cb','ad'})'/256;  % something a little brighter (peach puff 2 from http://www.tayloredmktg.com/rgb/)
            geom_ground.name = 'terrain';
            %      geom.c = hex2dec({'cd','af','95'})'/256;
            
            geom = cell(obj.num_boxes + 1, 1);
            geom{1} = geom_ground;
            for i=1:obj.num_boxes
                geom_box = RigidBodyBox([obj.boxs_x(i); obj.boxs_y(i); obj.boxs_z(i)]);
                geom_box.T(1:3, 4) = [obj.centers_x(i); obj.centers_y(i); obj.height(i)-obj.boxs_z(i)/2];
                geom_box.c = hex2dec({'ee','cb','ad'})'/256;
                geom_box.name = 'stair';
                geom{i+1} = geom_box;
            end
        end
    end
    
    properties
        geom;
        centers_x;
        centers_y;
        boxs_x;
        boxs_y;
        boxs_z;
        num_boxes;
        height;
    end
end