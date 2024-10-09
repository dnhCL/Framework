function fXi = calc_fXi(nF, cCoord, fNbCLoc, fNbC)

    % Inicializamos la matriz para almacenar los vectores fXi
    fXi = zeros(2, nF);

    % Recorremos cada cara
    for iF = 1:nF
        
        % Obtenemos los índices de las celdas adyacentes
        c1 = fNbC(fNbCLoc * (iF - 1) + 1); % Primera celda
        c2 = fNbC(fNbCLoc * (iF - 1) + 2); % Segunda celda

        % Obtenemos las coordenadas de los centros de las celdas
        c1Coord = cCoord(:, c1);  % Coordenadas del centro de la celda c1
        c2Coord = cCoord(:, c2);  % Coordenadas del centro de la celda c2
        
        % Calculamos el vector Xi que va desde c1 hacia c2
        xiVector = c2Coord - c1Coord;
        
        % Asignamos el vector xi al resultado
        fXi(:, iF) = xiVector;
        
        % Aquí es donde se debe verificar que xi . n >= 0 y, de ser necesario,
        % ajustar el vector xi. (Esto se haría si tuviéramos el vector normal fNormal)
    end

end





