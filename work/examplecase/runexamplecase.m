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

%% 1. Initialization
clear variables  % Clear all variables in memory
clc             % Clear the command window

%% 2. Mesh Generation
% Create structured mesh using line seeds in x and y directions
seedI = LineSeed.lineSeedOneWayBias([0 0], [1 0], 400, 1.00, 'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0], [0 1], 400, 1.00, 'o');

% Define boundary names
casedef.boundarynames = {'WESTRAND', 'OOSTRAND', 'ZUIDRAND', 'NOORDRAND'};

% Generate the mesh
mesh = TwoSeedMesher.genmesh(seedI, seedJ, casedef.boundarynames);

% Create computational domain
casedef.dom = newdomain(mesh, 'MyDomain');

%% 3. Define Initial Fields
% Temperature field (scalar)
T = Field(casedef.dom.allCells, 0);
randomdata = rand(T.elsize, T.elcountzone) - 0.5; % Random initial temperature distribution
set(T, randomdata);
casedef.T = T;

% Velocity field (vector)
U = Field(casedef.dom.allCells, 1);
velocity_magnitude = 0;
x_component = velocity_magnitude;
y_component = 0;
set(U, [x_component * ones(1, U.elcountzone); y_component * ones(1, U.elcountzone)]);
casedef.U = U;

%% 4. Define Material Properties
casedef.material.k = 1;   % Thermal conductivity
casedef.material.rho = 1; % Density
casedef.material.mu = 1;  % Dynamic viscosity

%% 5. Define Boundary Conditions for Temperature
jBC = 0;

% Dirichlet boundary condition at WESTRAND (T = 10)
jBC = jBC + 1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 10;

% Dirichlet boundary condition at OOSTRAND (T = 0)
jBC = jBC + 1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;

% Dirichlet boundary condition at ZUIDRAND (T = 0)
jBC = jBC + 1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;

% Dirichlet boundary condition at NOORDRAND (T = 0)
jBC = jBC + 1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;

%% 6. Iteration Parameters
casedef.iteration.maxniter = 1000; % Maximum number of iterations
casedef.iteration.TTol = 1e-6;     % Convergence tolerance

%% 8. Define Normal, Tangent, and Xi Fields for Face Calculations
normal = Field(casedef.dom.allFaces, 1);
tangent = Field(casedef.dom.allFaces, 1);
xi = Field(casedef.dom.allFaces, 1);
set(normal, (casedef.dom.fNormal));
set(tangent, (casedef.dom.fTangent));
set(xi, (casedef.dom.fXi));

%% 9. Solve the Case Using the Solver
result = examplesolver(casedef);

%% 10. Post-Processing: Plot Temperature Field
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T, scale, 0); % Plot temperature field

%% 11. Extract and Plot Temperature Profiles

% Get total mesh size
Lx = double(seedI.displX); % Length in x-direction
Ly = double(seedJ.displY); % Length in y-direction

% Number of cells in each direction
Nx = double(seedI.nSegm);
Ny = double(seedJ.nSegm);

% Compute cell size
dx = Lx / Nx;
dy = Ly / Ny;

% Compute target locations for profile extraction
x_target = double((Lx / 2) - (dx / 2)); % Midpoint in x-direction
y_target = double((Ly / 2) - (dy / 2)); % Midpoint in y-direction
tolerance = double(dx / 2); % Define margin for selection

% Get cell centers and temperature values
cell_centers = double(casedef.dom.cCoord);
temperatures = double(result.T.data);

% Select cells near x_target (vertical temperature profile)
cells_in_line_x = abs(cell_centers(1, :) - x_target) < tolerance;
y_values = cell_centers(2, cells_in_line_x);
temperature_profile_x = temperatures(cells_in_line_x);

% Filter values within a valid range
valid_range_y = (y_values >= 0) & (y_values <= Ly);
y_values = y_values(valid_range_y);
temperature_profile_x = temperature_profile_x(valid_range_y);

% Sort data by y-coordinates
[y_values, sort_idx] = sort(y_values);
temperature_profile_x = temperature_profile_x(sort_idx);

% Plot temperature profile at x = Lx/2
figure;
plot(y_values, temperature_profile_x, '-o');
xlabel('y');
ylabel('Temperature');
title(['Temperature Profile at x = ', num2str(Lx/2), ' (0 ≤ y ≤ ', num2str(Ly), ')']);
grid on;

% ----------------------------

% Select cells near y_target (horizontal temperature profile)
cells_in_line_y = abs(cell_centers(2, :) - y_target) < tolerance;
x_values = cell_centers(1, cells_in_line_y);
temperature_profile_y = temperatures(cells_in_line_y);

% Filter values within a valid range
valid_range_x = (x_values >= 0) & (x_values <= Lx);
x_values = x_values(valid_range_x);
temperature_profile_y = temperature_profile_y(valid_range_x);

% Sort data by x-coordinates
[x_values, sort_idx] = sort(x_values);
temperature_profile_y = temperature_profile_y(sort_idx);

% Plot temperature profile at y = Ly/2
figure;
plot(x_values, temperature_profile_y, '-o');
xlabel('x');
ylabel('Temperature');
title(['Temperature Profile at y = ', num2str(Ly/2), ' (0 ≤ x ≤ ', num2str(Lx), ')']);
grid on;



