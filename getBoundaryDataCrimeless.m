function [N,A,H1,H2] = getBoundaryDataCrimeless(model,results,M)
% Noise levels in percentages. Percentages of what? Who knows
thetaNoise = 0.05;
NNoise = 0.05;
resultsNoise = 0.1;
H2Noise = 0.05;

disp(['Angle noise: ', num2str(thetaNoise*100), '%'])
disp(['N noise: ', num2str(NNoise*100), '%'])
disp(['H noise: ', num2str(H2Noise*100), '%'])


% Evaluate gradient of solution and conductivity at near-unit circle
% (corrected to stay inside the mesh area)
[edge_midpoints, ~, ~] = getedges(model); 

rng(12)

h = 2*pi/M;
theta = 0:h:2*pi-h;
theta = theta + randn(size(theta))*h*thetaNoise;
N = length(edge_midpoints);
s = cos(pi/N); % Noise in this?
x = s*cos(theta);
y = s*sin(theta);
%results.NodalSolution = results.NodalSolution + resultsNoise*randn(size(results.NodalSolution)).*results.NodalSolution;
[gradx,grady] = evaluateGradient(results,x,y);
grad = [gradx grady];
c_vals = (c(x,y))';

edge_parallels = [y; -x];
edge_normals = [x;y];
edgeparanorm = sqrt(edge_parallels(1,:).^2 + edge_parallels(2,:).^2);
edge_parallels = edge_parallels./edgeparanorm;
edgenormalnorm = sqrt(edge_normals(1,:).^2 + edge_normals(2,:).^2);
edge_normals = edge_normals./edgenormalnorm;

plotData(theta,c_vals,'Exact solution','\theta','\sigma','exactSolution')

% Construct Neumann/Dirichlet data
N = c_vals.*sum(grad.*edge_normals',2); % Add direct noise
N = N + NNoise*N.*randn(size(N));
A = abs(sum(grad.*edge_parallels',2)); % Add indirect noise by adding noise to results before gradient calculation
H1 = c_vals.*sqrt(sum(grad.*grad,2));
H2 = c_vals.*sum(grad.*grad,2); % Add direct noise - Gaussian noise (L^2)
H2 = H2 + H2Noise*H2.*randn(size(H2));

end