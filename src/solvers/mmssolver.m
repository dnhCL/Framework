%==========================================================================
% MMS Solver Using the FVMLab Framework
%
% Objective:
% - Implements a solver for a **scalar conservation equation**.
% - Uses the **Finite Volume Method (FVM)** to obtain velocity distributions.
% - Applies **boundary conditions** and iterates the solution until convergence.
%
% Author: Frederik Rogiers
%==========================================================================

function result = mmssolver(casedef)

    %% 1. Extract Domain and Initialize Variables
    dom = casedef.dom;
    
    % Obtain the initial velocity field
    U = casedef.vars.Uinit;
    
    % Extract material properties from the problem definition
    mu = casedef.material.mu;  % Dynamic viscosity
    dt = casedef.iteration.dt;  % Time step value

    % Create equation objects to store numerical discretization
    eqn_u = ScalarFvEqn2(dom);
    eqn_v = ScalarFvEqn2(dom);
    
    %% 2. Setup Iteration Loop
    iterate = true;
    niter = 0;
    
    while iterate
    
      niter = niter + 1;
    
      % Reset equation objects before updating with new values
      reset(eqn_u);
      reset(eqn_v);
    
      % Initialize coefficient matrices for the system of equations
      adiag_u = zeros(dom.nC, 1);
      bdata_u = zeros(dom.nC, 1);
      anb_internal_u = zeros(2 * dom.nIf, 1);
      anb_boundary_u = zeros(2 * dom.nBf, 1);
    
      adiag_v = zeros(dom.nC, 1);
      bdata_v = zeros(dom.nC, 1);
      anb_internal_v = zeros(2 * dom.nIf, 1);
      anb_boundary_v = zeros(2 * dom.nBf, 1);
    
      % Compute interpolated velocity at faces based on adjacent cells
      u_face = faceInterpolate(dom, U);
      u_normal = scalarProduct(u_face, dom.fNormal);
    
      %% 3. Compute Coefficients for Internal Faces
      for i = 1:dom.nF  % Traverse through all faces in the domain
        c1 = dom.fNbC(2 * i - 1); % First adjacent cell to the face
        c2 = dom.fNbC(2 * i); % Second adjacent cell to the face
        xiMag = dom.fXiMag(i); % Distance between the centers of adjacent cells
        fArea = dom.fArea(i); % Surface area of the face
        fu_normal = u_normal(i); % Normal velocity component at the face
    
        % Compute contributions from convection and diffusion
        conv = -fArea * fu_normal / 2; % Convective effect
        D = -mu * fArea / xiMag; % Diffusive effect
    
        if i <= dom.nIf
          % Modify diagonal coefficients for internal cells
          adiag_u(c1) = adiag_u(c1) - D - conv;
          adiag_u(c2) = adiag_u(c2) - D + conv;
          
          adiag_v(c1) = adiag_v(c1) - D - conv;
          adiag_v(c2) = adiag_v(c2) - D + conv;
          
          % Update off-diagonal coefficients connecting adjacent cells
          anb_internal_u(2 * i - 1) = D - conv;
          anb_internal_u(2 * i) = D + conv;
          
          anb_internal_v(2 * i - 1) = D - conv;
          anb_internal_v(2 * i) = D + conv;
        else
          % Modify coefficients for boundary cells
          adiag_u(c1) = adiag_u(c1) - D - conv;
          anb_boundary_u(2 * (i - dom.nIf) - 1) = D - conv;
          adiag_v(c1) = adiag_v(c1) - D - conv;
          anb_boundary_v(2 * (i - dom.nIf) - 1) = D - conv;
        end
      end
      
      %% 4. Apply Time-Stepping and Source Terms
      for i = 1:dom.nPc
        c1 = i;
        cVol = dom.cVol(c1);
        fTime = cVol / dt; % Contribution from time discretization
        source_u = casedef.S.Src.data(1, c1); % External force in x-direction
        source_v = casedef.S.Src.data(2, c1); % External force in y-direction
        
        % Adjust diagonal coefficients to account for time-dependent terms
        adiag_u(c1) = adiag_u(c1) + fTime;
        adiag_v(c1) = adiag_v(c1) + fTime;
        
        % Incorporate source contributions
        bdata_u(c1) = bdata_u(c1) + U.data(1, c1) * fTime + source_u * cVol;
        bdata_v(c1) = bdata_v(c1) + U.data(2, c1) * fTime + source_v * cVol;
      end
    
      %% 5. Implement Boundary Conditions
      for j = 1:length(casedef.BC)  % Traverse through all boundary regions
        b_zone = dom.getzone(casedef.BC{j}.zoneID);
        r_zone = b_zone.range(1):b_zone.range(2);
        kind = casedef.BC{j}.kind;
        bcval = casedef.BC{j}.data.bcval;
    
        for i = r_zone
          c1 = dom.fNbC(2 * i - 1);
          cGc = dom.fNbC(2 * i); % Index for the ghost cell
          fx = dom.fCoord(1, i);
          fy = dom.fCoord(2, i);
    
          if (kind == "Dirichlet")
            % Assign fixed velocity values at the boundaries
            bdata_u(cGc) = casedef.S.sol_x(fx, fy);
            bdata_v(cGc) = casedef.S.sol_y(fx, fy);
            adiag_u(cGc) = 0.5;
            adiag_v(cGc) = 0.5;
            anb_boundary_u(2 * (i - dom.nIf)) = 0.5;
            anb_boundary_v(2 * (i - dom.nIf)) = 0.5;
          elseif (kind == "Neumann")
            % Implement a zero-gradient boundary condition
            fXiMag = dom.fXiMag(i);
            bdata_u(cGc) = (casedef.S.Sol.data(1, cGc) - casedef.S.Sol.data(1, c1)) / fXiMag;
            bdata_v(cGc) = (casedef.S.Sol.data(2, cGc) - casedef.S.Sol.data(2, c1)) / fXiMag;
            adiag_u(cGc) = 1 / fXiMag;
            adiag_v(cGc) = 1 / fXiMag;
            anb_boundary_u(2 * (i - dom.nIf)) = -1 / fXiMag;
            anb_boundary_v(2 * (i - dom.nIf)) = -1 / fXiMag;
          else
            fprintf("Unrecognized boundary condition type in zone %s", b_zone.id);
          end
        end
      end
    
      %% 6. Solve the System of Equations
      eqn_u.adata = [adiag_u; anb_internal_u; anb_boundary_u];
      eqn_u.bdata = bdata_u;
      
      eqn_v.adata = [adiag_v; anb_internal_v; anb_boundary_v];
      eqn_v.bdata = bdata_v;
    
      [Au, bu] = to_msparse(eqn_u);
      [Av, bv] = to_msparse(eqn_v);
      u = U.data(1, :)';
      v = U.data(2, :)';
    
      %% 7. Check for Convergence
      uRes = bu - Au * u;
      vRes = bv - Av * v;
      UResnorm = max([norm(uRes), norm(vRes)]);
      
      if UResnorm < casedef.iteration.TTol || niter > casedef.iteration.maxniter
        iterate = false;
        Uconverged = (UResnorm < casedef.iteration.TTol);
      else
        u = Au \ bu;
        v = Av \ bv;
        set(U, [u'; v']); 
      end
    end 
    
    %% 8. Save Computed Results
    result.endtime = now;
    result.Uconverged = Uconverged;
    result.niter = niter;
    result.UResnorm = UResnorm;
    result.U = U;
    
end

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