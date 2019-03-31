function [edge_midpoints, edge_parallels, edge_normals] = getedges(model)
edgenodes = findNodes(model.Mesh,'region','Edge',[1 2 3 4]);
edgenodes = [edgenodes, edgenodes(1)];
edgecoords = model.Mesh.Nodes(:,edgenodes);


edge_midpoints = (edgecoords(:,1:end-1) + edgecoords(:,2:end))/2;
edge_parallels = edgecoords(:,2:end) - edgecoords(:,1:end-1);
edgeparanorm = sqrt(edge_parallels(1,:).^2 + edge_parallels(2,:).^2);
edge_parallels = edge_parallels./edgeparanorm;
edge_normals = [edge_parallels(2,:); -edge_parallels(1,:)];
edgenormalnorm = sqrt(edge_normals(1,:).^2 + edge_normals(2,:).^2);
edge_normals = edge_normals./edgenormalnorm;
end