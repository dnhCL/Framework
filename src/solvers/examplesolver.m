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
 
    % Crear los objetos de campo
    T = Field(dom.allCells, 0);      % Temperatura [K] (escalar); campo vacío
    reset(T, 0);                     % Reiniciar con ceros
 
    % Crear un objeto de ecuación para contener la ecuación de conservación escalar
    eqn = ScalarFvEqn2(dom);
 
    % Definir propiedades del material
    k = casedef.material.k;  % Conductividad térmica [W/(m K)]
 
    % Configuración de la iteración
    iterate = true;
    niter = 0;
    while iterate
        niter = niter + 1;
 
        % Reiniciar todos los términos a cero
        reset(eqn);

        % Inicializar los coeficientes
        adiag = zeros(dom.nC, 1);
        anb_internal = zeros(2 * dom.nIf, 1);
        anb_boundary = zeros(2 * dom.nBf, 1);
        bdata = zeros(eqn.n, 1);

        % Calcular coeficientes para ecuaciones de celdas físicas y agregarlos al objeto eqn
        internal_idx = 1;
        boundary_idx = 1;
        for jF = 1:dom.nF
            % Obtener las celdas vecinas de la cara
            c1 = dom.fNbC((jF - 1) * dom.fNbCLoc + 1);
            c2 = dom.fNbC((jF - 1) * dom.fNbCLoc + 2);
 
            % Obtener el área de la cara y la distancia entre centroides
            Af = dom.fArea(jF);
            Lf = dom.fXiMag(jF);
 
            % Coeficiente de difusión
            D = k * Af / Lf;
 
            % Si c2 es una celda fantasma, entonces estamos en un borde
            if c2 > dom.nPc
                % Determinar el valor de la condición de borde usando casedef.BC
                bcval = 0; % Valor por defecto
                bcType = 'Dirichlet'; % Tipo por defecto
 
                % Buscar la condición de borde aplicable a la cara actual
                for jBC = 1:length(casedef.BC)
                    % Obtenemos el identificador de la zona de frontera
                    boundaryZone = getzone(dom, casedef.BC{jBC}.zoneID);
                    % Verificar si la cara pertenece a la zona de frontera usando el rango de la zona
                    if jF >= boundaryZone.range(1) && jF <= boundaryZone.range(2)
                        bcval = casedef.BC{jBC}.data.bcval;
                        if isfield(casedef.BC{jBC}.data, 'bcType')
                            bcType = casedef.BC{jBC}.data.bcType;
                        end
                        break;
                    end
                end

                if strcmp(bcType, 'Dirichlet')
                    % Aplicar la condición de borde Dirichlet a la celda física (c1)
                    adiag(c1) = adiag(c1) + D;
                    bdata(c1) = bdata(c1) + D * bcval;
                elseif strcmp(bcType, 'Neumann')
                    % Aplicar la condición de borde Neumann (flujo especificado)
                    bdata(c1) = bdata(c1) + bcval * Af;
                end
                
                anb_boundary(boundary_idx) = -D;
                anb_boundary(boundary_idx + 1) = -D;
                boundary_idx = boundary_idx + 2;
            else
                % Caso de cara interna: agregar los coeficientes de las dos celdas vecinas
                adiag(c1) = adiag(c1) + D;
                adiag(c2) = adiag(c2) + D;
                anb_internal(internal_idx) = -D;
                anb_internal(internal_idx + 1) = -D;
                internal_idx = internal_idx + 2;
            end
        end

        % Asegurar simetría y completar la diagonal de la matriz
        for i = 1:dom.nC
            if adiag(i) == 0
                adiag(i) = -sum(anb_internal(internal_idx - 2:internal_idx - 1)) - sum(anb_boundary(boundary_idx - 2:boundary_idx - 1));
            end
        end

        % Asignar los coeficientes a la ecuación
        eqn.adata = [adiag; anb_internal; anb_boundary];
        eqn.bdata = bdata;
 
        % Crear un sistema lineal esparso de matlab a partir del objeto eqn
        [A, b] = to_msparse(eqn);
        x = get(T);
        x = x';
        spy(A);
        title('Estructura de la matriz de coeficientes A');
        
        % Comprobar la tolerancia y el conteo de iteraciones
        TRes = b - A * x;
        TResnorm = norm(TRes);
        if TResnorm < casedef.iteration.TTol
            Tconverged = true;
            iterate = false;
        elseif niter > casedef.iteration.maxniter
            Tconverged = false;
            iterate = false;
        else
            % Resolver el sistema lineal
            x = A \ b;
            set(T, x'); % Colocar la solución algebraica en el campo
        end
    end
 
    % Preparar la estructura de resultados
    result.endtime = now; % Llamar a datestr(now) para mostrar esta hora
    result.Tconverged = Tconverged;
    result.niter = niter;
    result.TResnorm = TResnorm;
    result.TRes = Field(dom.allCells, 0);
    set(result.TRes, TRes');
    result.T = T;
 
    % Graficar el campo de temperatura resultante
    figure;
    hold on;
    axis off;
    axis equal;
    colormap(jet(50));
    fvmplotfield(result.T, 'lin', 0);
    title('Distribución de la Temperatura en el Dominio');

end



 




