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
    
    % Create an equation object for holding a scalar conservation equation
    eqn = ScalarFvEqn2(dom);
    
    % Extraer información del dominio y parámetros
    fNbC = dom.fNbC;
    fNbCLoc = dom.fNbCLoc; % Dimensión por cara (ej. 2 en 2D)
    fArea = dom.fArea;     % Área de cada cara
    fXiMag = dom.fXiMag;   % Distancia entre centros de celdas vecinas
    kappa = casedef.material.k; % Coeficiente de difusión
    fXiLambda = dom.fXiLambda;

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
           
           % Calcular el coeficiente de difusión entre c1 y c2
           D = -kappa * fArea(i) / fXiMag(i);
           
           % Almacenar en anb_internal para ambos lados (simétrico)
           anb_internal(2 * (i - 1) + 1) = D; % (c1 -> c2)
           anb_internal(2 * (i - 1) + 2) = D; % (c2 -> c1)
       end

       % Ensamblaje de coeficientes en caras de borde
       for j = 1:dom.nBf
           % Identificar la celda física y la celda fantasma para la cara de borde
           cP = fNbC(fNbCLoc * (dom.nIf + j - 1) + 1); % Celda física
           cGC = fNbC(fNbCLoc * (dom.nIf + j - 1) + 2); % Celda fantasma
           
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
                % Ajustar el término de fuente y el valor de φ_GC según la interpolación
                
                % Expresar φ_GC en términos de φ_PC y φ*
                anb_boundary(fNbCLoc * (j - 1) + 1) = D  ; % (cP -> cGC, PDE en celda física)
                anb_boundary(fNbCLoc * (j - 1) + 2) = D ; % (cGC -> cP, ajuste por Dirichlet)
                adiag(cP) = adiag(cP) - anb_boundary(fNbCLoc * (j - 1) + 1) ; % Diagonal de la celda física
                adiag(cGC) = adiag(cGC) - anb_boundary(fNbCLoc * (j - 1) + 2) / (1 - fXiLambda(dom.nIf + j)) ;
                % Ajustar bdata en función de φ* según la interpolación de Dirichlet
                eqn.bdata(cP) = eqn.bdata(cP) - D *bcval * (1.05 - fXiLambda(dom.nIf + j)) ; 

            elseif strcmp(kind, 'Neumann')
                % Condición de Neumann: imponer un gradiente en la frontera
                anb_boundary(fNbCLoc * (j - 1) + 1) = D; % (cP -> cGC, PDE en celda física)
                anb_boundary(fNbCLoc * (j - 1) + 2) = D; % (cGC -> cP, Condición Neumann)
                adiag(cP) = adiag(cP) - anb_boundary(fNbCLoc * (j - 1) + 1) ; % Diagonal de la celda física
                adiag(cGC) = adiag(cGC) - anb_boundary(fNbCLoc * (j - 1) + 2) ; % Diagonal de la celda fantasma
                % Ajustar bdata en función del gradiente φ'
                eqn.bdata(cP) = eqn.bdata(cP) + bcval * fArea(dom.nIf + j);
            end
       end

       % Calcular adiag después de ensamblar anb_internal y anb_boundary
        

        % Recorrer las caras internas y añadir contribuciones a la diagonal
        for i = 1:dom.nIf
            % Identificar las celdas vecinas para la cara i
            c1 = fNbC(fNbCLoc * (i - 1) + 1);
            c2 = fNbC(fNbCLoc * (i - 1) + 2);
            
            % Añadir contribuciones de anb_internal a la diagonal de c1 y c2
            adiag(c1) = adiag(c1) - anb_internal(2 * (i - 1) + 1) ;
            adiag(c2) = adiag(c2) - anb_internal(2 * (i - 1) + 2) ;
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


    % Modificación para imprimir temperatura en celdas de borde y fantasma
    fprintf('Temperaturas en las celdas de borde y fantasma:\n');
    for j = 1:dom.nBf
        % Identificar las celdas física y fantasma para la cara de borde j
        cP = fNbC(fNbCLoc * (dom.nIf + j - 1) + 1); % Celda física
        cGC = fNbC(fNbCLoc * (dom.nIf + j - 1) + 2); % Celda fantasma
        
        % Imprimir la temperatura en la celda física y en la celda fantasma
        fprintf('Cara %d - Celda física %d: T = %.4f, Celda fantasma %d: T = %.4f\n', ...
                dom.nIf + j, cP, result.T.data(cP), cGC, result.T.data(cGC));
        fprintf('Celda física %d (Cara %d): bdata = %.4f\n', cP, dom.nIf + j, eqn.bdata(cP));
    end

end



    
    










    










 




