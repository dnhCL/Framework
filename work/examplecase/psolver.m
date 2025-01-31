function result = psolver(Nx, Ny)
    %======================================================================
    % CFD Pressure Solver Case Function (Modified for Grid Convergence Study)
    %======================================================================

    %% 1. Initialization
    % clear variables; clc;

    %% 2. Define Initial and Boundary Conditions
    uniformU = [0, 0];  % Initial velocity field (zero in both directions)

    % Define pressure gradient direction
    f = [1, 0];

    % Define boundary velocity and pressure gradient
    u_out = 10;     % Velocity at boundary
    dp = -1000;     % Pressure difference
    gradP = dp * f; % Compute pressure gradient

    % Define boundary conditions
    casedef.boundarynames = {'WEST', 'EAST', 'SOUTH', 'NORTH'};
    BCtype = {'Neumann', 'Neumann', 'Dirichlet', 'Dirichlet'};
    BCval = [0 0; 0 0; 0 0; u_out 0];

    %% 3. Generate Mesh
    seedI = LineSeed.lineSeedOneWayBias([0, 0], [1, 0], Nx, 1.00, 'o');
    seedJ = LineSeed.lineSeedOneWayBias([0, 0], [0, 1], Ny, 1.00, 'o');
    mesh = TwoSeedMesher.genmesh(seedI, seedJ, casedef.boundarynames);
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
    casedef.material.mu = 0.1;
    casedef.material.rho = 1;
    casedef.material.k = 1;

    % Store variables in the case definition
    casedef.U = U;
    casedef.P = P;
    casedef.gradP = gradP;

    %% 5. Apply Boundary Conditions
    for i = 1:length(casedef.boundarynames)
        casedef.BC{i}.zoneID = casedef.boundarynames(i);
        casedef.BC{i}.kind = BCtype(i);
        casedef.BC{i}.data.bcval = BCval(i, :);
    end

    %% 6. Set Iteration Parameters
    casedef.iteration.maxniter = 1000;
    casedef.iteration.TTol = 1e-9;
    casedef.iteration.dt = 50;

    %% 7. Solve the Problem using pressuresolver
    result = pressuresolver(casedef);

    %% 8. Compute Relative Error with Analytical Solution
    % Get total mesh size
    Lx = 1.0;
    Ly = 1.0;
    dx = Lx / Nx;
    dy = Ly / Ny;

    % Compute target location for velocity extraction
    x_target = (Lx / 2) - (dx / 2);
    y_target = (Ly / 2) - (dy / 2);
    tolerance = dx / 2;

    % Get cell centers and velocity values
    cell_centers = double(casedef.dom.cCoord);
    velocity_x = double(result.U.data(1, :));

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
                                      casedef.material.mu * casedef.material.rho, u_out);

    % Compute relative error
    rel_error = norm(velocity_profile_x - u_analytical) / norm(u_analytical);

    % Store relative error in result structure
    result.dom = casedef.dom;  
    result.rel_error = rel_error;
    result.Nx = Nx;
    result.Ny = Ny;

end

function u = channel(ymesh, dp, H, mu, Utop)
    u = dp * H^2 / (2 * mu) * (ymesh.^2 / H^2 - ymesh / H) + Utop * ymesh / H;
end
