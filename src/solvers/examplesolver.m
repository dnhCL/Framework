%==========================================================================
%
% Example solver using the FVMLab framework 4 students
%
% Purpose: Provides code structure for solving a scalar conservation
%          equation using data structures provided by the framework.
%
% by Frederik Rogiers
%
%==========================================================================
function result = examplesolver(casedef)

    dom = casedef.dom;
    
    % Create field objects
    T = Field(dom.allCells, 0);      % Temperature [K] (scalar); empty field
    reset(T, 0);                     % Reset with all zeros

    P = Field(dom.allCells, 0);      % Presión [Pa] (scalar); empty field
    reset(P, 0);                     % Reset with all zeros

    
    % Create an equation object for holding a scalar conservation equation
    eqn = ScalarFvEqn2(dom);
    
    % Extraer información del dominio y parámetros
    fNbC = dom.fNbC;
    fNbCLoc = dom.fNbCLoc; % Dimensión por cara (ej. 2 en 2D)
    fArea = dom.fArea;     % Área de cada cara
    fXiMag = dom.fXiMag;   % Distancia entre centros de celdas vecinas
    kappa = casedef.material.k; % Coeficiente de difusión
    fXiLambda = dom.fXiLambda;
    U = casedef.U;
    %rho = casedef.material.rho;

    iterate = true;
    niter = 0;
    while iterate   
       
       niter = niter + 1;
       
       % Reset all terms in the equation object for the new iteration
       reset(eqn); 
       
       % Inicializar vectores para almacenar los coeficientes
       adiag = zeros(dom.nC, 1);
       anb_internal = zeros(2 * dom.nIf, 1); % Coeficientes de caras internas
       anb_boundary = zeros(2 * dom.nBf, 1); % Coeficientes de caras de borde
       
       % Ensamblaje de coeficientes en caras internas
       for i = 1:dom.nIf
            % Identificar las celdas vecinas para la cara i
            c1 = fNbC(fNbCLoc * (i - 1) + 1);
            c2 = fNbC(fNbCLoc * (i - 1) + 2);

            normal = dom.fNormal(:, i); % Vector normal a la cara
            area = fArea(i); % Área de la cara

            % Velocidad en la cara (promedio de las celdas adyacentes)
            U_face = 0.5 * (U.data(:, c1) + U.data(:, c2));

            % Flujo convectivo
            F_conv = (U_face' * normal) * area *0 ;

            % Gradiente de presión
            pressure_grad = (P.data(c2) - P.data(c1)) / fXiMag(i);
            pressure_term = pressure_grad * normal * area; % Contribución del gradiente de presión
           
            % Calcular el coeficiente de difusión entre c1 y c2
            D = -kappa * fArea(i) / fXiMag(i);

            % Coeficientes para flujo convectivo y difusión
            a_c1_c2 = D + max(0, F_conv) ; % Upwind para convección
            a_c2_c1 = D - min(0, F_conv) ; % Upwind para convección
           
            % Almacenar en anb_internal para ambos lados (simétrico)
            anb_internal(2 * (i - 1) + 1) = a_c1_c2; % (c1 -> c2)
            anb_internal(2 * (i - 1) + 2) = a_c2_c1; % (c2 -> c1)

            % Añadir contribuciones de anb_internal a la diagonal de c1 y c2
            adiag(c1) = adiag(c1) - anb_internal(2 * (i - 1) + 1) ;
            adiag(c2) = adiag(c2) - anb_internal(2 * (i - 1) + 2) ;

           
       end

       % Ensamblaje de coeficientes en caras de borde
       for j = 1:dom.nBf
           % Identificar la celda física y la celda fantasma para la cara de borde
           cP = fNbC(fNbCLoc * (dom.nIf + j - 1) + 1); % Celda física
           cGC = fNbC(fNbCLoc * (dom.nIf + j - 1) + 2); % Celda fantasma

           normal = dom.fNormal(:, dom.nIf + j); % Vector normal a la cara
           area = fArea(dom.nIf + j); % Área de la cara

           % Velocidad en la cara (promedio de las celdas adyacentes)
           U_face = U.data(:, cP);

           % Flujo convectivo
           F_conv = (U_face' * normal) * area *0 ;
           
           % Calcular el coeficiente entre la celda física y la celda fantasma
           D = -kappa * fArea(dom.nIf + j) / fXiMag(dom.nIf + j);
           
           
           % Determinar la condición de borde para cGC
           kind = '';
           bcval = 0;
            % Buscar la condición de borde aplicable a la cara actual
            for jBC = 1:length(casedef.BC)
                % Obtenemos el identificador de la zona de frontera
                boundaryZone = getzone(dom, casedef.BC{jBC}.zoneID);
                % Verificar si la cara pertenece a la zona de frontera usando el rango de la zona
                if j + dom.nIf >= boundaryZone.range(1) && j + dom.nIf  <= boundaryZone.range(2)
                    bcval = casedef.BC{jBC}.data.bcval;
                    kind = casedef.BC{jBC}.kind;  % 'Dirichlet' o 'Neumann'
                    break;
                end
            end
           
           % Aplicar la condición de borde según el tipo
            if strcmp(kind, 'Dirichlet')
                % Condición de Dirichlet: interpolar entre φ_PC y φ_GC para que cumpla φ_f = φ*
                
                    % Flujo entrante
                    anb_boundary(fNbCLoc * (j - 1) + 1) =  (D + F_conv)  ; % (cP -> cGC, PDE en celda física)
                    anb_boundary(fNbCLoc * (j - 1) + 2) =  fXiLambda(dom.nIf + j); % (cGC -> cP, Dirichlet)
                    adiag(cP) = adiag(cP) - anb_boundary(fNbCLoc * (j - 1) + 1);
                    adiag(cGC) = (1 - fXiLambda(dom.nIf + j)) ;
                
                eqn.bdata(cGC) = eqn.bdata(cGC) + bcval; 

            elseif strcmp(kind, 'Neumann')
                % Ajustar los coeficientes dependiendo del flujo

                    % Flujo entrante
                    anb_boundary(fNbCLoc * (j - 1) + 1) = D + F_conv ; % (cP -> cGC, PDE en celda física)
                    anb_boundary(fNbCLoc * (j - 1) + 2) = D - F_conv ; % (cGC -> cP, Neumann)
                    adiag(cP) = adiag(cP) - anb_boundary(fNbCLoc * (j - 1) + 1);
                    adiag(cGC) = adiag(cGC) - anb_boundary(fNbCLoc * (j - 1) + 2)  ;

                % Ajustar bdata en función del gradiente φ'
                eqn.bdata(cP) = eqn.bdata(cP) + bcval * fArea(dom.nIf + j);
            end
       end

                  
        
                 
        


       % Asignar los datos calculados a eqn
       
       eqn.adata = [adiag; anb_internal; anb_boundary];
       
       % Crear sistema disperso en MATLAB desde el objeto eqn
       [A, b] = to_msparse(eqn);
       x = get(T);
       x = x';
       spy(A);
       
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
          x = A \ b; % Solución directa
          set(T, x'); % Guardar la solución en el campo
       end
         
    end % iteración

    % Resultados finales
    result.endtime = now; % Tiempo final
    result.Tconverged = Tconverged;
    result.niter = niter;
    result.TResnorm = TResnorm;
    result.TRes = Field(dom.allCells, 0);
    set(result.TRes, TRes');
    result.T = T;
    disp(result.T.data(105));
    disp(result.T.data(115));
    disp(result.T.data(125));
    disp(result.T.data(135));

    




end



    
    










    










 




