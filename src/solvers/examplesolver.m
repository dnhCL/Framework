

% Convection-Diffusion Solver
function result = examplesolver(casedef)

    %% 1. Extract Domain and Initialize Fields
    dom = casedef.dom;
    
    % Temperature field [K] (scalar)
    T = Field(dom.allCells, 0);
    reset(T, 0); % Initialize with zeros

    % Pressure field [Pa] (scalar)
    P = Field(dom.allCells, 0);
    reset(P, 0); % Initialize with zeros

    % Create equation object for scalar conservation equation
    eqn = ScalarFvEqn2(dom);

    %% 2. Extract Domain Properties and Parameters
    fNbC = dom.fNbC;       % Face neighbor cells
    fNbCLoc = dom.fNbCLoc; % Dimension per face (e.g., 2 in 2D)
    fArea = dom.fArea;     % Area of each face
    fXiMag = dom.fXiMag;   % Distance between neighboring cell centers
    kappa = casedef.material.k; % Diffusion coefficient
    fXiLambda = dom.fXiLambda;
    U = casedef.U;         % Velocity field
    rho = casedef.material.rho; % Density

    % Initialize lists for storing face velocities
    U_face_internal = [];
    U_face_boundary = [];

    %% 3. Iteration Loop
    iterate = true;
    niter = 0;
    
    while iterate   
       
        niter = niter + 1;

        % Reset equation object for new iteration
        reset(eqn);
        
        % Initialize coefficient vectors
        adiag = zeros(dom.nC, 1);
        anb_internal = zeros(2 * dom.nIf, 1); % Internal face coefficients
        anb_boundary = zeros(2 * dom.nBf, 1); % Boundary face coefficients

        U_face_internal = zeros(1, dom.nIf); % Reset for each iteration
        U_face_boundary = zeros(1, dom.nBf); % Reset for each iteration
        
        %% 4. Assemble Coefficients for Internal Faces
        for i = 1:dom.nIf
            % Identify neighboring cells for face i
            c1 = fNbC(fNbCLoc * (i - 1) + 1);
            c2 = fNbC(fNbCLoc * (i - 1) + 2);

            normal = dom.fNormal(:, i); % Face normal vector
            area = fArea(i); % Face area

            % Compute velocity at the face (average of adjacent cells)
            U_face = 0.5 * (U.data(:, c1) + U.data(:, c2));
            
            % Compute normal component of U_face
            U_normal_face = dot(U_face, normal);
            U_face_internal(i) = U_normal_face;

            % Compute convective flux
            F_conv = U_normal_face * area;

            % Compute diffusion coefficient between c1 and c2
            D = -kappa * fArea(i) / fXiMag(i);

            % Compute Peclet number
            Pe = abs(F_conv) / D;

            % Hybrid Scheme Selection
            if Pe <= 2
                % Central Differencing Scheme (CDS)
                a_c1_c2 = D + 0.5 * F_conv;
                a_c2_c1 = D - 0.5 * F_conv;
            else
                % Upwind Scheme
                if F_conv > 0
                    a_c1_c2 = D + F_conv;
                    a_c2_c1 = D;
                else
                    a_c1_c2 = D;
                    a_c2_c1 = D - F_conv;
                end
            end
                     
            % Store coefficients in anb_internal (symmetric)
            anb_internal(fNbCLoc * (i - 1) + 1) = a_c1_c2;
            anb_internal(fNbCLoc * (i - 1) + 2) = a_c2_c1;

            % Add contributions to diagonal coefficients
            adiag(c1) = adiag(c1) - anb_internal(fNbCLoc * (i - 1) + 1);
            adiag(c2) = adiag(c2) - anb_internal(fNbCLoc * (i - 1) + 2);
        end

        %% 5. Assemble Coefficients for Boundary Faces
        for j = 1:dom.nBf
            % Identify physical and ghost cells for boundary face j
            cP = fNbC(fNbCLoc * (dom.nIf + j - 1) + 1); % Physical cell
            cGC = fNbC(fNbCLoc * (dom.nIf + j - 1) + 2); % Ghost cell    
            
            % Compute diffusion coefficient for boundary face
            D = -kappa * fArea(dom.nIf + j) / fXiMag(dom.nIf + j);
            
            % Determine boundary condition type
            kind = '';
            bcval = 0;
            for jBC = 1:length(casedef.BC)
                boundaryZone = getzone(dom, casedef.BC{jBC}.zoneID);
                if j + dom.nIf >= boundaryZone.range(1) && j + dom.nIf <= boundaryZone.range(2)
                    bcval = casedef.BC{jBC}.data.bcval;
                    kind = casedef.BC{jBC}.kind;
                    break;
                end
            end
            
            % Apply boundary conditions
            if strcmp(kind, 'Dirichlet')
                % Dirichlet Condition (φ_f = φ*)
                anb_boundary(fNbCLoc * (j - 1) + 1) = D;
                anb_boundary(fNbCLoc * (j - 1) + 2) = fXiLambda(dom.nIf + j);
                adiag(cP) = adiag(cP) - anb_boundary(fNbCLoc * (j - 1) + 1);
                adiag(cGC) = (1 - fXiLambda(dom.nIf + j));
                eqn.bdata(cGC) = eqn.bdata(cGC) + bcval;

            elseif strcmp(kind, 'Neumann')
                % Neumann Condition (Gradient Specified)
                anb_boundary(fNbCLoc * (j - 1) + 1) = D;
                anb_boundary(fNbCLoc * (j - 1) + 2) = D;
                adiag(cP) = adiag(cP) - anb_boundary(fNbCLoc * (j - 1) + 1);
                adiag(cGC) = adiag(cGC) - anb_boundary(fNbCLoc * (j - 1) + 2);
                eqn.bdata(cP) = eqn.bdata(cP) + bcval * fArea(dom.nIf + j);
            end
        end        

        % Store coefficients in equation object
        eqn.adata = [adiag; anb_internal; anb_boundary];
       
        % Convert equation system to MATLAB sparse matrix
        [A, b] = to_msparse(eqn);
        x = get(T);
        x = x';


        % Compute residuals and check convergence
        TRes = b - A * x;
        TResnorm = norm(TRes);
        
        if TResnorm < casedef.iteration.TTol
            Tconverged = true;
            iterate = false;
        elseif niter > casedef.iteration.maxniter
            Tconverged = false;
            iterate = false;
        else
            x = A \ b; % Direct solution
            set(T, x'); % Store solution in field
        end
    end

    %% 6. Store Final Results
    result.endtime = now;
    result.Tconverged = Tconverged;
    result.niter = niter;
    result.TResnorm = TResnorm;
    result.TRes = Field(dom.allCells, 0);
    set(result.TRes, TRes');
    result.T = T;
    result.P = P;
    result.U = U;

end


    










 




