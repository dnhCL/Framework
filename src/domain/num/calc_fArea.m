function fArea = calc_fArea(nF, fNbVLoc, fNbV, vCoord)

    % Verificación del número de caras (nF)
    fprintf('Número de caras (nF) en calc_fArea: %d\n', nF);

    % Inicializamos el vector para almacenar las áreas (longitudes) de las caras
    fArea = zeros(1, nF);
    
    % Cada cara tiene 2 vértices
    numVerticesPerFace = fNbVLoc;  % Esto es 2 según la documentación
    
    % Verificación del tamaño de fNbV
    totalVertices = length(fNbV);
    fprintf('Tamaño de fNbV en calc_fArea: %d\n', totalVertices);
    
    if totalVertices == 0
        error('El vector fNbV está vacío en calc_fArea.');
    end
    
    % Recorremos cada cara
    for iF = 1:nF
        
        % El índice de los vértices para la cara iF está dado por:
        startIdx = (iF - 1) * numVerticesPerFace + 1;  % Índice inicial en fNbV
        endIdx = startIdx + numVerticesPerFace - 1;    % Índice final en fNbV
        
        % Verificación de que los índices están dentro de los límites
        if endIdx > totalVertices
            error('Índice fuera de los límites en calc_fArea: startIdx = %d, endIdx = %d, totalVertices = %d', startIdx, endIdx, totalVertices);
        end
        
        % Extraemos los índices de los vértices de la cara
        faceVertices = fNbV(startIdx:endIdx);
        
        % Extraemos las coordenadas de los dos vértices
        v1 = vCoord(:, faceVertices(1));
        v2 = vCoord(:, faceVertices(2));
        
        % Calculamos la longitud de la cara como la distancia entre los dos vértices
        fArea(iF) = sqrt((v2(1) - v1(1))^2 + (v2(2) - v1(2))^2);
    end
end







