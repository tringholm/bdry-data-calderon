function [model,results] = fwdProblem(meshres)
model = createpde();
geometryFromEdges(model,@circleg); % Specify geometry, just using circular domain for now
%pdegplot(model,'EdgeLabels','on'); % Plot geometry to see if it's all right
% axis tight
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',@bdrycond,'Vectorized','on'); % Set boundary conditions, dirichlet for now
specifyCoefficients(model,'m',0,'d',0,'c',@c_coeff,'a',0,'f',1); % We need different c and f, all else is zero
generateMesh(model,'Hmax',meshres,'GeometricOrder','quadratic'); % Generate FEM mesh
results = solvepde(model); % Solve the PDE
figure(10);
u = results.NodalSolution;
pdeplot(model,'XYData',u,'ZData',u)
colormap jet
end