%Pressure 
function result = examplesolver2(casedef)

  %% 1. Extract Domain and Initialize Fields
  dom = casedef.dom;
  
  % Initialize velocity field
  U = casedef.vars.Uinit;
  
  % Extract material properties
  nu = casedef.material.nu;  % Kinematic viscosity
  rho = casedef.material.rho; % Density
  dt = casedef.iteration.dt;  % Time step size
  gradP = casedef.vars.gradP; % Pressure gradient

  % Create equation objects for velocity components (u and v)
  eqn_u = ScalarFvEqn2(dom);
  eqn_v = ScalarFvEqn2(dom);

  %% 2. Iteration Loop
  iterate = true;
  niter = 0;
  
  while iterate
  
      niter = niter + 1;
  
      % Reset equation objects for a new iteration
      reset(eqn_u);
      reset(eqn_v);
  
      % Initialize coefficient vectors for the linear system
      adiag_u = zeros(dom.nC, 1);
      bdata_u = zeros(dom.nC, 1);
      anb_internal_u = zeros(2 * dom.nIf, 1);
      anb_boundary_u = zeros(2 * dom.nBf, 1);
  
      adiag_v = zeros(dom.nC, 1);
      bdata_v = zeros(dom.nC, 1);
      anb_internal_v = zeros(2 * dom.nIf, 1);
      anb_boundary_v = zeros(2 * dom.nBf, 1);
  
      % Compute face velocities from the previous iteration
      U_face = faceInterpolate(dom, U);
      U_norm = scalarProduct(U_face, dom.fNormal);
  
      %% 3. Compute Coefficients for Internal Faces
      for i = 1:dom.nF  % Loop over all faces
          c1 = dom.fNbC(2 * i - 1); % Physical cell index
          c2 = dom.fNbC(2 * i);    % Neighboring cell index
          fXiMag = dom.fXiMag(i);    % Distance between cell centers
          fArea = dom.fArea(i);     % Face area
          fU_norm = U_norm(i);        % Normal velocity component at the face
          c1_c2 = dom.cVol([c1, c2]);

          % Compute convective and Dusive terms
          conv = -fArea * fU_norm / 2; % Convective term
          D = -nu * fArea / fXiMag; % Diffsive term

          % Update diagonal coefficients for velocity equations
          if i <= dom.nIf  % Internal faces
              adiag_u(c1) = adiag_u(c1) - D - conv + (c1_c2(1) / dt) / 4;
              adiag_u(c2) = adiag_u(c2) - D + conv + (c1_c2(2) / dt) / 4;
              adiag_v(c1) = adiag_v(c1) - D - conv + (c1_c2(1) / dt) / 4;
              adiag_v(c2) = adiag_v(c2) - D + conv + (c1_c2(2) / dt) / 4;

              anb_internal_u(2 * i - 1) = D - conv;
              anb_internal_u(2 * i) = D + conv;
              anb_internal_v(2 * i - 1) = D - conv;
              anb_internal_v(2 * i) = D + conv;

              % Compute source terms for velocity equations
              bdata_u(c1) = bdata_u(c1) + (U.data(1, c1) * c1_c2(1) / dt - gradP(1) / rho * c1_c2(1)) / 4;
              bdata_v(c1) = bdata_v(c1) + (U.data(2, c1) * c1_c2(1) / dt - gradP(2) / rho * c1_c2(1)) / 4;
              bdata_u(c2) = bdata_u(c2) + (U.data(1, c2) * c1_c2(2) / dt - gradP(1) / rho * c1_c2(2)) / 4;
              bdata_v(c2) = bdata_v(c2) + (U.data(2, c2) * c1_c2(2) / dt - gradP(2) / rho * c1_c2(2)) / 4;
            else
              adiag_u(c1) = adiag_u(c1) - D - conv + (c1_c2(1)/dt)/4;
              anb_boundary_u(2*(i - dom.nIf)-1) = D - conv;
              adiag_v(c1) = adiag_v(c1) - D - conv + (c1_c2(1)/dt)/4;
              anb_boundary_v(2*(i - dom.nIf)-1) = D - conv;
              bdata_u(c1) = bdata_u(c1) + (U.data(1,c1) * c1_c2(1) / dt - gradP(1)/rho * c1_c2(1)) /4;
              bdata_v(c1) = bdata_v(c1) + (U.data(2,c1) * c1_c2(1) / dt - gradP(2)/rho * c1_c2(1)) /4;
            end
          
      end

      %% 4. Apply Boundary Conditions
      for j = 1:length(casedef.BC)
          b_zone = dom.getzone(casedef.BC{j}.zoneID);
          r_zone = b_zone.range(1):b_zone.range(2);
          kind = casedef.BC{j}.kind;
          bcval = casedef.BC{j}.data.bcval;

          for i = r_zone
              cGC = dom.fNbC(2 * i);
              bdata_u(cGC) = bcval(1);
              bdata_v(cGC) = bcval(2);

              if strcmp(kind, "Dirichlet")
                  adiag_u(cGC) = 0.5;
                  adiag_v(cGC) = 0.5;
                  anb_boundary_u(2 * (i - dom.nIf)) = 0.5;
                  anb_boundary_v(2 * (i - dom.nIf)) = 0.5;
              elseif strcmp(kind, "Neumann")
                  adiag_u(cGC) = 1 / fXiMag;
                  adiag_v(cGC) = 1 / fXiMag;
                  anb_boundary_u(2 * (i - dom.nIf)) = -1 / fXiMag;
                  anb_boundary_v(2 * (i - dom.nIf)) = -1 / fXiMag;
              else
                  fprintf("Unknown boundary condition for zone %s", b_zone.id);
              end
          end
      end
    
      %% 5. Solve the Linear System
      eqn_u.adata = [adiag_u; anb_internal_u; anb_boundary_u];
      eqn_u.bdata = bdata_u;
    
      eqn_v.adata = [adiag_v; anb_internal_v; anb_boundary_v];
      eqn_v.bdata = bdata_v;
  
      [Au, bu] = to_msparse(eqn_u);
      [Av, bv] = to_msparse(eqn_v);
  
      u = U.data(1, :)';
      v = U.data(2, :)';
  
      %% 6. Convergence Check
      uRes = bu - Au * u;
      vRes = bv - Av * v;
      UResnorm = max([norm(uRes), norm(vRes)]);

      if UResnorm < casedef.iteration.TTol
          Uconverged = true;
          iterate = false;
      elseif niter > casedef.iteration.maxniter
          Uconverged = false;
          iterate = false;
      else
          u = Au \ bu;
          v = Av \ bv;
          set(U, [u'; v']); % Store the solution in the velocity field
      end
  end % End iteration loop

  %% 7. Store Final Results
  result.endtime = now;
  result.Uconverged = Uconverged;
  result.niter = niter;
  result.UResnorm = UResnorm;
  result.uRes = Field(dom.allCells, 0); set(result.uRes, uRes');
  result.vRes = Field(dom.allCells, 0); set(result.vRes, vRes');
  result.U = U;

end

function fieldFace=faceInterpolate(dom,field)
    tempFace=zeros(field.dim,dom.nF);
    for i=1:dom.nF
      c1=dom.fNbC(2*i-1); %physical cell index
      c2=dom.fNbC(2*i); %neighbor cell
      fXiLambda=dom.fXiLambda(i);
      tempFace(:,i)=( (1-fXiLambda)*field.data(:,c1) + fXiLambda*field.data(:,c2) );
    end
    
    fieldFace = Field(dom.allFaces,1);     % Velocity [m/s] (vector);
    set(fieldFace,tempFace);
    
end