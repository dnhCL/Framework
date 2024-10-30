%==========================================================================
%
% Solver for Steady-State Diffusion using the FVMLab Framework
%
% Purpose: Implements the steady-state diffusion equation in a scalar
%          conservation equation using finite volume discretization.
%          Computes the coefficients and solves the resulting system.
%
% by Frederik Rogiers (modified for steady-state diffusion)
%
%==========================================================================
function result = examplesolver(casedef)

    dom = casedef.dom;
 
    % Create field objects
    T = Field(dom.allCells, 0);      % Temperature [K] (scalar); empty field
    reset(T, 0);                     % Reset to zeros
 
    % Create an equation object to contain the scalar conservation equation
    eqn = ScalarFvEqn2(dom);
 
    % Define material properties
    k = casedef.material.k;  % Thermal conductivity [W/(m K)]
 
    % Iteration setup
    iterate = true;
    niter = 0;
    while iterate
        niter = niter + 1;
 
        % Reset all terms to zero
        reset(eqn);

        % Initialize coefficients
        adiag = zeros(dom.nC, 1);
        anb_internal = zeros(2 * dom.nIf, 1);
        anb_boundary = zeros(2 * dom.nBf, 1);
        bdata = zeros(eqn.n, 1);

        % Compute coefficients for physical cell equations and add them to the eqn object
        internal_idx = 1;
        boundary_idx = 1;
        for jF = 1:dom.nF
            % Obtain neighboring cells of the face
            c1 = dom.fNbC((jF - 1) * dom.fNbCLoc + 1);
            c2 = dom.fNbC((jF - 1) * dom.fNbCLoc + 2);
 
            % Obtain face area and distance between centroids
            Af = dom.fArea(jF);
            Lf = dom.fXiMag(jF);
 
            % Diffusion coefficient
            D = k * Af / Lf;
 
            % If c2 is a ghost cell, then we are at a boundary
            if c2 > dom.nPc
                % Determine boundary condition value using casedef.BC
                bcval = 0; % Default value
                bcType = 'Dirichlet'; % Default type
 
                % Search for applicable boundary condition for the current face
                for jBC = 1:length(casedef.BC)
                    % Obtain the boundary zone identifier
                    boundaryZone = getzone(dom, casedef.BC{jBC}.zoneID);
                    % Check if face belongs to the boundary zone using the zone range
                    if jF >= boundaryZone.range(1) && jF <= boundaryZone.range(2)
                        bcval = casedef.BC{jBC}.data.bcval;
                        if isfield(casedef.BC{jBC}.data, 'bcType')
                            bcType = casedef.BC{jBC}.data.bcType;
                        end
                        break;
                    end
                end

                if strcmp(bcType, 'Dirichlet')
                    % Apply Dirichlet boundary condition to physical cell (c1)
                    adiag(c1) = adiag(c1) + D;
                    bdata(c1) = bdata(c1) + D * bcval;
                elseif strcmp(bcType, 'Neumann')
                    % Apply Neumann boundary condition (specified flux)
                    bdata(c1) = bdata(c1) + bcval * Af;
                end
                
                anb_boundary(boundary_idx) = -D;
                anb_boundary(boundary_idx + 1) = -D;
                boundary_idx = boundary_idx + 2;
            else
                % Internal face case: add coefficients of the two neighboring cells
                adiag(c1) = adiag(c1) + D;
                adiag(c2) = adiag(c2) + D;
                anb_internal(internal_idx) = -D;
                anb_internal(internal_idx + 1) = -D;
                internal_idx = internal_idx + 2;
            end
        end

        % Ensure symmetry and complete matrix diagonal
        for i = 1:dom.nC
            if adiag(i) == 0
                adiag(i) = -sum(anb_internal(internal_idx - 2:internal_idx - 1)) - sum(anb_boundary(boundary_idx - 2:boundary_idx - 1));
            end
        end

        % Assign coefficients to the equation
        eqn.adata = [adiag; anb_internal; anb_boundary];
        eqn.bdata = bdata;
 
        % Create a sparse linear system in MATLAB from the eqn object
        [A, b] = to_msparse(eqn);
        x = get(T);
        x = x';
        spy(A);
        title('Coefficient Matrix Structure A');
        
        % Check tolerance and iteration count
        TRes = b - A * x;
        TResnorm = norm(TRes);
        if TResnorm < casedef.iteration.TTol
            Tconverged = true;
            iterate = false;
        elseif niter > casedef.iteration.maxniter
            Tconverged = false;
            iterate = false;
        else
            % Solve the linear system
            x = A \ b;
            set(T, x'); 
        end
    end
 
    % Prepare result structure
    result.endtime = now; % Call datestr(now) to display this time
    result.Tconverged = Tconverged;
    result.niter = niter;
    result.TResnorm = TResnorm;
    result.TRes = Field(dom.allCells, 0);
    set(result.TRes, TRes');
    result.T = T;
 
    % Plot the resulting temperature field (this is not important)
    figure;
    hold on;
    axis off;
    axis equal;
    colormap(jet(50));
    fvmplotfield(result.T, 'lin', 0);
    title('Temperature Distribution in the Domain');

end



 




