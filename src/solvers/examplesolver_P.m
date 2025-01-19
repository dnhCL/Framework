function result = examplesolver_with_pressure(casedef)

    dom = casedef.dom;

    % Create field objects
    T = Field(dom.allCells, 0);      % Temperature [K] (scalar); empty field
    reset(T, 0);                     % Reset with all zeros

    P = Field(dom.allCells, 0);      % Pressure [Pa] (scalar); empty field
    % Initialize pressure with a small perturbation to avoid singularity
    reset(P, rand(size(P.data)) * 1e-5);

    % Create an equation object for holding a scalar conservation equation
    eqn = ScalarFvEqn2(dom);

    % Extract domain and parameters
    fNbC = dom.fNbC;
    fNbCLoc = dom.fNbCLoc;
    fArea = dom.fArea;
    fXiMag = dom.fXiMag;
    kappa = casedef.material.k; % Diffusion coefficient
    rho = casedef.material.rho; % Fluid density
    U = casedef.U;

    % Iteration parameters
    iterate = true;
    niter = 0;
    maxniter = casedef.iteration.maxniter;
    tol = casedef.iteration.TTol;

    while iterate

        niter = niter + 1;

        % Reset equation object
        reset(eqn);

        % Initialize coefficients
        adiag = zeros(dom.nC, 1);
        anb_internal = zeros(2 * dom.nIf, 1);
        anb_boundary = zeros(2 * dom.nBf, 1);

        % Compute gradient of pressure at cell centers
        gradP_cell = zeros(2, dom.nC);
        for c = 1:dom.nC
            gradP_cell(:, c) = compute_pressure_gradient(P, dom, c, fNbC, fNbCLoc);
        end

        % Interpolate pressure gradient to faces (Rhie-Chow)
        gradP_face = zeros(2, dom.nIf + dom.nBf);
        for i = 1:dom.nIf
            c1 = fNbC(fNbCLoc * (i - 1) + 1);
            c2 = fNbC(fNbCLoc * (i - 1) + 2);

            if c1 > 0 && c2 > 0 && c1 <= dom.nC && c2 <= dom.nC
                % Linear interpolation of gradP
                gradP_face(:, i) = 0.5 * (gradP_cell(:, c1) + gradP_cell(:, c2));
            end
        end

        % Add gradP contribution to internal faces
        for i = 1:dom.nIf
            c1 = fNbC(fNbCLoc * (i - 1) + 1);
            c2 = fNbC(fNbCLoc * (i - 1) + 2);

            if c1 > 0 && c2 > 0 && c1 <= dom.nC && c2 <= dom.nC
                normal = dom.fNormal(:, i); % Normal vector at face
                area = fArea(i);

                % Pressure force contribution
                F_pressure = -(gradP_face(:, i)' * normal) * area / rho;

                % Add to coefficients
                adiag(c1) = adiag(c1) - F_pressure;
                adiag(c2) = adiag(c2) + F_pressure;
            end
        end

        % Ensure diagonal dominance for boundary cells
        for b = 1:dom.nBf
            boundary_cell = fNbC(fNbCLoc * (dom.nIf + b - 1) + 1);
            if boundary_cell > 0 && boundary_cell <= dom.nC
                adiag(boundary_cell) = 1;
                anb_boundary(fNbCLoc * (b - 1) + 2) = 0;
            end
        end

        % Solve for velocity (predictor step)
        [A, b] = to_msparse(eqn);
        for b = 1:dom.nBf
            boundary_cell = fNbC(fNbCLoc * (dom.nIf + b - 1) + 1);
            if boundary_cell > 0 && boundary_cell <= dom.nC
                b(boundary_cell) = 0; % Set boundary conditions in RHS vector
            end
        end

        u = get(U);
        u = u';
        u = A \ b;
        set(U, u');

        % Pressure correction (SIMPLE algorithm)
        % Compute divergence of velocity for continuity
        div_u = compute_velocity_divergence(U, dom);

        % Solve for pressure correction
        dp = solve_pressure_correction(div_u, dom, rho);
        p = get(P);
        p = p + dp; % Update pressure
        set(P, p);

        % Correct velocities with updated pressure
        correct_velocities(U, dp, dom, rho);

        % Check convergence
        if norm(div_u) < tol || niter > maxniter
            iterate = false;
        end
    end

    % Save results
    result.P = P;
    result.U = U;
    result.niter = niter;
    result.converged = norm(div_u) < tol;

end

% Helper functions
function gradP = compute_pressure_gradient(P, dom, cellID, fNbC, fNbCLoc)
    % Compute the gradient of pressure for a given cell
    gradP = zeros(2, 1);
    for i = 1:dom.nIf
        % Identify neighboring cells for the face
        c1 = fNbC(fNbCLoc * (i - 1) + 1);
        c2 = fNbC(fNbCLoc * (i - 1) + 2);

        if (c1 == cellID || c2 == cellID) && c1 > 0 && c2 > 0 && c1 <= dom.nC && c2 <= dom.nC
            % Calculate pressure difference and gradient
            distance = dom.fXiMag(i);
            if c1 == cellID
                gradP = gradP + (P.data(c2) - P.data(c1)) / distance;
            else
                gradP = gradP + (P.data(c1) - P.data(c2)) / distance;
            end
        end
    end
end

function div_u = compute_velocity_divergence(U, dom)
    % Compute divergence of velocity field
    div_u = zeros(1, dom.nC);
    for c = 1:dom.nC
        for i = 1:dom.nIf
            % Identify neighboring cells for the face
            c1 = dom.fNbC(dom.fNbCLoc * (i - 1) + 1);
            c2 = dom.fNbC(dom.fNbCLoc * (i - 1) + 2);
            
            if (c1 == c || c2 == c) && c1 > 0 && c2 > 0 && c1 <= dom.nC && c2 <= dom.nC
                faceID = i;
                normal = dom.fNormal(:, faceID);
                div_u(c) = div_u(c) + dot(U.data(:, faceID), normal);
            end
        end
    end
end

function dp = solve_pressure_correction(div_u, dom, rho)
    % Solve pressure correction equation
    % Example placeholder implementation (to be replaced with actual solve)
    dp = -div_u / rho;
end

function correct_velocities(U, dp, dom, rho)
    % Correct velocities using pressure correction
    for f = 1:dom.nIf
        c1 = dom.fNbC(f * 2 - 1);
        c2 = dom.fNbC(f * 2);
        if c1 > 0 && c2 > 0 && c1 <= dom.nC && c2 <= dom.nC
            dp_face = (dp(c1) + dp(c2)) / 2;
            U.data(:, f) = U.data(:, f) - dp_face / rho;
        end
    end
end