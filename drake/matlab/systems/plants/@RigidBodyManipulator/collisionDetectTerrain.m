function [phi,normal,xB,idxB] = collisionDetectTerrain(obj, xA_in_world)
%   if sum(strcmp(fieldnames(obj), 'terrain_index'))%check whether there is a field named as terrain_index
%       [z,normal] = getHeight(obj.terrain_sample{obj.terrain_index},xA_in_world(1:2,:));
%   else
      [z,normal] = getHeight(obj.terrain,xA_in_world(1:2,:));
%   end
  xB = [xA_in_world(1:2,:);z];
  idxB = ones(1,size(xA_in_world,2));
  %phi = sqrt(sum((xA_in_world-xB).^2))';
  phi = (xA_in_world(3,:)-z)';
end
