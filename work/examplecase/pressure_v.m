%==========================================================================
% Simulation of Flow in a Channel with Fixed Pressure Gradient
% Method: Finite Volume Method (FVM)
% The numerical solution is compared with the analytical solution of 
% channel flow.
%==========================================================================

%% 1. Initialization
clear variables  % Clear all variables in memory
% clc           % Commented out to avoid clearing the console
close all       % Close all open figures

%% 2. Define Initial and Boundary Conditions
uniformU = [0, 0];  % Initial velocity field (zero in both directions)

% Define pressure gradient direction
f = [1, 0] ;

% Define boundary velocity and pressure gradient
u_out = 10;     % Velocity at boundary
dp = -1000;         % Pressure difference
gradP = dp * f;  % Compute pressure gradient

% Define boundary conditions
casedef.boundarynames = { 'WEST', 'EAST', 'SOUTH', 'NORTH' };
BCtype = { 'Neumann', 'Neumann', 'Dirichlet', 'Dirichlet' };
BCval = [0 0; 0 0; 0 0; u_out 0];

%% 3. Generate Mesh
% Define grid seeds in the x and y directions with 10 uniform divisions
seedI = LineSeed.lineSeedOneWayBias([0, 0], [1, 0], 10, 1.00, 'o');
seedJ = LineSeed.lineSeedOneWayBias([0, 0], [0, 1], 10, 1.00, 'o');

% Generate mesh from the seeds
mesh = TwoSeedMesher.genmesh(seedI, seedJ, casedef.boundarynames);

% Define computational domain
casedef.dom = newdomain(mesh, 'MyDomain');

%% 4. Define Physical Fields
% Velocity Field [m/s] (Vector)
U = Field(casedef.dom.allCells, 1);
set(U, [uniformU(1) * ones(1, U.elcountzone); 
        uniformU(2) * ones(1, U.elcountzone)]);

% Pressure Field [Pa] (Scalar)
P = Field(casedef.dom.allCells, 0);
reset(P, 0);

% Define material properties
casedef.material.mu = 0.1;     % Kinematic viscosity
casedef.material.rho = 10;      % Density
casedef.material.k = 1;        % Thermal conductivity [W/(m·K)]

% Store variables in the case definition
casedef.vars.Uinit = U;
casedef.vars.P = P;  
casedef.vars.gradP = gradP;

%% 5. Apply Boundary Conditions
for i = 1:length(casedef.boundarynames)
    casedef.BC{i}.zoneID = casedef.boundarynames(i);
    casedef.BC{i}.kind   = BCtype(i);
    casedef.BC{i}.data.bcval = BCval(i,:);
end

%% 6. Set Iteration Parameters
casedef.iteration.maxniter = 1000;  % Maximum number of iterations
casedef.iteration.TTol     = 1e-9;  % Convergence tolerance
casedef.iteration.dt       = 50;    % Time step

%% 7. Solve the Problem using examplesolver2
result = pressuresolver(casedef);

%% 8. Post-processing: Extract and Plot Velocity Field
% Define velocity field for post-processing
u = Field(casedef.dom.allCells, 0);
set(u, result.U.data(1, :));  % Extract numerical solution

% Plot the velocity field
figure; hold on; axis off; axis equal; colormap(jet(50));
title("Values of U");
scale = 'lin'; lw = 0; quiver = 1;
fvmplotfield(u, scale, lw);  % Plot scalar velocity field
%% 9. Compare Numerical Solution with Analytical Solution

% Get total mesh size
Lx = double(seedI.displX); % Length in x-direction
Ly = double(seedJ.displY); % Length in y-direction

% Number of elements in each direction
Nx = double(seedI.nSegm); % Number of cells in x
Ny = double(seedJ.nSegm); % Number of cells in y

% Cell size
dx = Lx / Nx;
dy = Ly / Ny;

% Compute target location for velocity extraction
x_target = double((Lx / 2) - (dx / 2)); % Middle of the domain in x-direction
y_target = double((Ly / 2) - (dy / 2)); % Middle of the domain in y-direction
tolerance = double(dx / 2); % Use half dx as margin for selection

% Get cell centers and velocity values
cell_centers = double(casedef.dom.cCoord); % Convert to double for precision
velocity_x = double(result.U.data(1, :)); % Extract U_x component (horizontal velocity)

% Select cells close to x_target (vertical profile at x = Lx/2)
cells_in_line_x = abs(cell_centers(1, :) - x_target) < tolerance;

% Extract y-coordinates and velocities from selected cells
y_values = cell_centers(2, cells_in_line_x);
velocity_profile_x = velocity_x(cells_in_line_x);

% Filter values within valid range
valid_range_y = (y_values >= 0) & (y_values <= Ly);
y_values = y_values(valid_range_y);
velocity_profile_x = velocity_profile_x(valid_range_y);

% Sort data by y-coordinates
[y_values, sort_idx] = sort(y_values);
velocity_profile_x = velocity_profile_x(sort_idx);

% Compute analytical solution
u_analytical = channel(y_values, -norm(gradP), Ly, ...
                                   casedef.material.nu * casedef.material.rho, u_out);

% Compute relative error
rel_error = norm(velocity_profile_x - u_analytical) / norm(u_analytical);
fprintf("Relative error of numerical solution: %10.4e\n", rel_error);

% Plot the velocity profile at x = Lx/2
figure;
plot(y_values, velocity_profile_x, '-o', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
plot(y_values, u_analytical, '*', 'MarkerSize', 6);
xlabel('y');
ylabel('U_x');
title(['Velocity Profile at x = ', num2str(Lx/2), ' (0 ≤ y ≤ ', num2str(Ly), ')']);
legend("Numerical Solution", "Analytical Solution", 'Location', 'northwest');
grid on;
hold off;

% ----------------------------

% Select cells close to y_target (horizontal profile at y = Ly/2)
cells_in_line_y = abs(cell_centers(2, :) - y_target) < tolerance;

% Extract x-coordinates and velocities from selected cells
x_values = cell_centers(1, cells_in_line_y);
velocity_profile_y = velocity_x(cells_in_line_y);

% Filter values within valid range
valid_range_x = (x_values >= 0) & (x_values <= Lx);
x_values = x_values(valid_range_x);
velocity_profile_y = velocity_profile_y(valid_range_x);

% Sort data by x-coordinates
[x_values, sort_idx] = sort(x_values);
velocity_profile_y = velocity_profile_y(sort_idx);

% Plot the velocity profile at y = Ly/2
figure;
plot(x_values, velocity_profile_y, '-o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('x');
ylabel('U_x');
title(['Velocity Profile at y = ', num2str(Ly/2), ' (0 ≤ x ≤ ', num2str(Lx), ')']);
grid on;

function u = channel(ymesh,dp,H,mu,Utop)

  u = dp*H^2/(2*mu)*( ymesh.^2/H^2 - ymesh/H ) + Utop * ymesh/H;
  
end