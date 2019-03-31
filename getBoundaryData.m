function [N,A,H1,H2] = getBoundaryData(model,results)
% Extract edge information from mesh - gives lists of the edges' midpoints,
% parallel vectors and normals
[edge_midpoints, edge_parallels, edge_normals] = getedges(model); 

% Evaluate gradient of solution and conductivity at edge midpoints
[gradx,grady] = evaluateGradient(results,edge_midpoints(1,:),edge_midpoints(2,:));
grad = [gradx grady];
c_vals = (c(edge_midpoints(1,:),edge_midpoints(2,:)))';

figure;
plot(c_vals); title('Exact solution')

% Construct Neumann/Dirichlet data
N = c_vals.*sum(grad.*edge_normals',2); % Add direct noise
A = abs(sum(grad.*edge_parallels',2)); % Add indirect noise by adding noise to results before gradient calculation
H1 = c_vals.*sqrt(sum(grad.*grad,2));
H2 = c_vals.*sum(grad.*grad,2); % Add direct noise - Gaussian noise (L^2)
end