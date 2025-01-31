function result = ssolver1(Nx, Ny)
    %======================================================================
    % CFD Master Case Function (Modified for Grid Convergence Study)
    %======================================================================

    %% 1. Initialization
    %clear variables; clc;

    %% 2. Mesh Generation
    seedI = LineSeed.lineSeedOneWayBias([0 0], [1 0], Nx, 1.00, 'o');
    seedJ = LineSeed.lineSeedOneWayBias([0 0], [0 1], Ny, 1.00, 'o');
    casedef.boundarynames = {'WESTRAND', 'OOSTRAND', 'ZUIDRAND', 'NOORDRAND'};
    mesh = TwoSeedMesher.genmesh(seedI, seedJ, casedef.boundarynames);
    casedef.dom = newdomain(mesh, 'MyDomain');

    %% 3. Define Initial Fields
    T = Field(casedef.dom.allCells, 0);
    set(T, rand(T.elsize, T.elcountzone) - 0.5);
    casedef.T = T;

    % Velocity field (vector)
    U = Field(casedef.dom.allCells, 1);
    velocity_magnitude = 0;
    x_component = velocity_magnitude;
    y_component = 0;
    set(U, [x_component * ones(1, U.elcountzone); y_component * ones(1, U.elcountzone)]);
    casedef.U = U;

    %% 4. Define Material Properties
    casedef.material.k = 1;
    casedef.material.rho = 1;
    casedef.material.mu = 1;

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
    casedef.iteration.maxniter = 1000;
    casedef.iteration.TTol = 1e-6;

    %% 7. Solve the Case Using the Solver
    result = examplesolver(casedef);

    % Return results for grid convergence analysis
    result.Nx = Nx;
    result.Ny = Ny;
    result.L2Error = norm(result.T.data - casedef.T.data) / norm(casedef.T.data);
end
