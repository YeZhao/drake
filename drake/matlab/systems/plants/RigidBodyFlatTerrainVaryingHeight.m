classdef RigidBodyFlatTerrainVaryingHeight < RigidBodyTerrain
%  This class provides an implementation of RigidBodyTerrain with varying z
%  everywhere
  
  methods 
    function obj = RigidBodyFlatTerrainVaryingHeight(height)
      % Construct a flat terrain map at a given z height (varying parameter)
      
      obj.z = 0;
      obj.geom = constructRigidBodyGeometry(obj,height);
    end
    
    function [z,normal] = getHeight(obj,xy)
      [m n] = size(xy);
      z = repmat(obj.z,1,n);
      normal = repmat([0;0;1],1,n);
    end

    function geom = getCollisionGeometry(obj)
      geom = obj.geom;
    end

    function geom = getVisualGeometry(obj)
      geom = obj.geom;
    end

    function obj = setGeometryColor(obj, color)
      geom = obj.getVisualGeometry();
      geom.c = reshape(color, 3, 1);
      obj.geom = geom;
    end
    
    function geom = constructRigidBodyGeometry(obj,height)
      box_width = 1.5;
      box_depth = 0.5;
      geom{1} = RigidBodyBox([box_width;box_width;box_depth]);
      geom{1}.T(1,4) = 0 - box_width/2 + 0.2;%assume 0.2 is slightly smaller than one step length
      geom{1}.T(3,4) = obj.z - box_depth/2;
      geom{1}.c = hex2dec({'ee','cb','ad'})'/256;  % something a little brighter (peach puff 2 from http://www.tayloredmktg.com/rgb/);
      geom{1}.name = 'terrain';
      
      % add a next-step rough terrain
      box_width2 = 2;
      box_depth2 = 0.5;
      geom{2} = RigidBodyBox([box_width2;box_width2;box_depth2]);
      geom{2}.T(1,4) = box_width2 - 0.8;
      geom{2}.T(3,4) = height - box_depth2/2;
      geom{2}.c = hex2dec({'ee','cb','ad'})'/256;  % something a little brighter (peach puff 2 from http://www.tayloredmktg.com/rgb/);
      geom{2}.name = 'terrain';
      
      % add a negative 1-step rough terrain
      box_width3 = 1.4;
      box_depth3 = 0.5;
      geom{3} = RigidBodyBox([box_width3;box_width3;box_depth3]);
      geom{3}.T(1,4) = -2;
      geom{3}.T(3,4) = obj.z - box_depth3/2 - 0.2;%0.2 is perturbed height
      geom{3}.c = hex2dec({'ee','cb','ad'})'/256;  % something a little brighter (peach puff 2 from http://www.tayloredmktg.com/rgb/);
      geom{3}.name = 'terrain';
      
%      geom.c = hex2dec({'cd','af','95'})'/256;
    end
  end
  
  properties
    geom;
    z;
  end
end
