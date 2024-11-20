%==========================================================================
%
% Example case using the FVMLab framework 4 students
%
% Purpose: Provides an example for setting up a case and calling a solver.
%          This involves creating a mesh, defining materials, defining
%          boundary conditions, defining iteration parameters, and finally
%          calling the solver.
%
% by Frederik Rogiers
%
%==========================================================================

% TIP: use "clear variables" instead of "clear all" to clear variables
%      use "clear classes" when the interface of a class has changes
%      use "close all" to close figures
%      use "clc" to clear the command window
% 
% TIP: pressing CTRL+D while the cursor is on a function opens that function
%      in the m-editor. This is the most convenient way of browsing through
%      your source code.

clear variables
clc


% Create a mesh
seedI = LineSeed.lineSeedOneWayBias([0 0],[1 0],10,1.00,'o');
seedJ = LineSeed.lineSeedOneWayBias([0 0],[0 1],10,1.00,'o');
casedef.boundarynames = {'WESTRAND','OOSTRAND','ZUIDRAND','NOORDRAND'};
mesh  = TwoSeedMesher.genmesh(seedI,seedJ,casedef.boundarynames);
% Create domain from mesh
casedef.dom = newdomain(mesh,'MyDomain');



% Set up initial fields
T = Field(casedef.dom.allCells,0);     % Temperature [K] (scalar); empty field
% reset(T,0);                          % Reset with all zeros
randomdata = rand(T.elsize,T.elcountzone)-0.5;
set(T,randomdata);                     % Set with random numbers

U = Field(casedef.dom.allCells,1);
set(U,[ones(1,U.elcountzone); zeros(1,U.elcountzone)]); % Velocidad constante en x

casedef.U = U;
%reset(U,[1;0.2]);




% Define material properties
casedef.material.k = 1;  % Thermal conductivity [W/(m K)]
%casedef.material.rho = 1;


% Define boundary conditions


jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 10;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0; 




%{
 jBC = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'WESTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 1;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'OOSTRAND';
casedef.BC{jBC}.kind   = 'Dirichlet';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'ZUIDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
jBC = jBC+1;
casedef.BC{jBC}.zoneID = 'NOORDRAND';
casedef.BC{jBC}.kind   = 'Neumann';
casedef.BC{jBC}.data.bcval = 0;
 
%}


% Set up iteration parameters
casedef.iteration.maxniter = 1000;
casedef.iteration.TTol     = 1e-6;

normal = Field(casedef.dom.allFaces,1);
tangent = Field(casedef.dom.allFaces,1);
xi = Field(casedef.dom.allFaces,1);
set(normal, (casedef.dom.fNormal));
set(tangent,(casedef.dom.fTangent));
set(xi,(casedef.dom.fXi));


% Call solver
result = examplesolver(casedef);
% Call solver
%result2 = examplesolver2(casedef);


% Plot result
figure; hold on; axis off; axis equal; colormap(jet(50));
scale = 'lin'; lw = 1;
fvmplotfield(result.T,scale,0);
%fvmplotfield(result2.T,scale,0);
% % Uoost = restrictto(U,getzone(casedef.dom,'OOSTRAND'));
% %fvmplotvectorfield(xi,lw);
fvmplotmesh(casedef.dom,lw);
%fvmplotcellnumbers(casedef.dom,8);
% fvmplotfacenumbers(casedef.dom,8);
% fvmplotvertexnumbers(casedef.dom,8);

% Obtener coordenadas de los centros de las celdas y los valores de temperatura
cell_centers = casedef.dom.cCoord; % Coordenadas de los centros de las celdas
temperatures = result.T.data;      % Temperatura en cada celda

% Seleccionar celdas cerca de x = 0.5 (con un pequeño margen de error)
tolerance = 0.05; % Ajusta este valor si necesitas mayor precisión
x_target = 0.45;
cells_in_line = abs(cell_centers(1, :) - x_target) < tolerance;

% Extraer coordenadas y temperaturas de las celdas seleccionadas
y_values = cell_centers(2, cells_in_line);   % Coordenada y de las celdas seleccionadas
temperature_profile = temperatures(cells_in_line); % Temperaturas correspondientes

% Filtrar valores de y entre 0 y 1
valid_range = (y_values >= -0) & (y_values <= 1);
y_values = y_values(valid_range);
temperature_profile = temperature_profile(valid_range);

% Ordenar por y
[y_values, sort_idx] = sort(y_values);
temperature_profile = temperature_profile(sort_idx);

% Graficar el perfil de temperaturas en x = 0.5
figure;
plot(y_values, temperature_profile, '-o');
xlabel('y');
ylabel('Temperatura');
title('Perfil de Temperatura a lo largo de x = 0.5 (0 ≤ y ≤ 1)');
grid on;


% Seleccionar celdas cerca de y = 0.5 (con un pequeño margen de error)

y_target = 0.45;
cells_in_line = abs(cell_centers(2, :) - y_target) < tolerance;

% Extraer coordenadas y temperaturas de las celdas seleccionadas
x_values = cell_centers(1, cells_in_line);   % Coordenada x de las celdas seleccionadas
temperature_profile = temperatures(cells_in_line); % Temperaturas correspondientes

% Filtrar valores de x entre 0 y 1
valid_range = (x_values >= 0) & (x_values <= 1);
x_values = x_values(valid_range);
temperature_profile = temperature_profile(valid_range);

% Ordenar por x
[x_values, sort_idx] = sort(x_values);
temperature_profile = temperature_profile(sort_idx);

% Graficar el perfil de temperaturas en y = 0.5
figure;
plot(x_values, temperature_profile, '-o');
xlabel('x');
ylabel('Temperatura');
title('Perfil de Temperatura a lo largo de y = 0.5 (0 ≤ x ≤ 1)');
grid on;
