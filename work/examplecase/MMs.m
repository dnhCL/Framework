%==========================================================================
% Example Case Using the FVMLab Framework
%
% Purpose:
% - Demonstrates how to set up a CFD case using the FVMLab framework.
% - Includes mesh generation, material property definition, boundary 
%   conditions, iteration parameters, and calling the solver.
%
% by Frederik Rogiers
%==========================================================================


clear variables
clc
close all

%% Case parameters
Lx = 1;
Nx = 100;
Ly = 1;
Ny = 100;
u_i = [0, 0]; % Initial velocity field

casedef.boundarynames = {'WEST', 'EAST', 'SOUTH', 'NORTH'};
BCtype = {'Dirichlet', 'Dirichlet', 'Neumann', 'Neumann'};
BCval = [0, 0; 0, 0; 0, 0; 0, 0];

%% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0, 0], [Lx, 0], Nx, 1.00, 'o');
seedJ = LineSeed.lineSeedOneWayBias([0, 0], [0, Ly], Ny, 1.00, 'o');
mesh  = TwoSeedMesher.genmesh(seedI, seedJ, casedef.boundarynames);

%% Create domain from mesh
casedef.dom = newdomain(mesh, 'MyDomain');

%% Set up initial fields
U = Field(casedef.dom.allCells, 1); % Velocity [m/s] (vector)
set(U, [u_i(1) * ones(1, U.elcountzone); u_i(2) * ones(1, U.elcountzone)]);

%% Define material properties
casedef.material.mu = 1;
casedef.material.rho = 1;

casedef.vars.Uinit = U;

%% Define boundary conditions
for i = 1:length(casedef.boundarynames)
  casedef.BC{i}.zoneID = casedef.boundarynames(i);
  casedef.BC{i}.kind = BCtype(i);
  casedef.BC{i}.data.bcval = BCval(i, :);
end

%% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol = 1e-9;
casedef.iteration.dt = 50;

%% MMS Function Definition (From Class 5)
    % solution
    sol_x = @(x, y) x.^2 ;
    sol_y = @(x, y) -y.^2;
    % source
    src_x = @(x, y) ( 4 .* x.^3 - 2 .* x.^2 .* y - 2 * casedef.material.mu ) ;
    src_y = @(x, y) ( -2 .* x .* y.^2 + 4 .* y.^3 + 2 * casedef.material.mu );
  
%% Assign MMS functions to the case definition
Src = Field(casedef.dom.allCells,1);
set(Src,[src_x(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:));
         src_y(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))])
casedef.S.Src = Src;
casedef.S.src_x = src_x;
casedef.S.src_y = src_y;

Sol = Field(casedef.dom.allCells,1);
set(Sol, [sol_x(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:));
          sol_y(casedef.dom.cCoord(1,:),casedef.dom.cCoord(2,:))])
casedef.S.Sol = Sol;
casedef.S.sol_x = sol_x;
casedef.S.sol_y = sol_y;

%% Call Solver
result = mmssolver(casedef);

%% Extract Results
u = Field(casedef.dom.allCells, 0);
set(u, result.U.data(1, :));
v = Field(casedef.dom.allCells, 0);
set(v, result.U.data(2, :));

uS = Field(casedef.dom.allCells, 0);
set(uS, casedef.S.Sol.data(1, :));
vS = Field(casedef.dom.allCells, 0);
set(vS, casedef.S.Sol.data(2, :));

%% Plot Results
figure;
hold on;
axis off;
axis equal;
colormap(jet(50));
scale = 'lin';
lw = 0;
quiver = 1;
fvmplotfield(uS, scale, lw);
fvmplotvectorfield(faceInterpolate(casedef.dom,casedef.S.Sol),quiver);
title('Chosen solution');

figure;
hold on;
axis off;
axis equal;
colormap(jet(50));
scale = 'lin';
lw = 0;
quiver = 1;
fvmplotfield(u, scale, lw);
fvmplotvectorfield(faceInterpolate(casedef.dom,U),quiver);
title('Result solver');

rel_error = norm(result.U.data - casedef.S.Sol.data)/norm(casedef.S.Sol.data);
fprintf("Relative error: %10.4e\n", rel_error);

%% Face Interpolation Function
function fieldFace=faceInterpolate(dom,field)
    tempFace=zeros(field.dim,dom.nF);
    for i=1:dom.nF
      c1=dom.fNbC(2*i-1);
      c2=dom.fNbC(2*i);
      fXiLambda=dom.fXiLambda(i);
      tempFace(:,i)=( (1-fXiLambda)*field.data(:,c1) + fXiLambda*field.data(:,c2) );
    end
    
    fieldFace = Field(dom.allFaces,1);
    set(fieldFace,tempFace);
end 

