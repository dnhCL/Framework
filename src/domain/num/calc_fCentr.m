function fCentr = calc_fCentr(nF, fNbVLoc, fNbV, vCoord)

    % Verificación del número de caras (nF)
    fprintf('Número de caras (nF) en calc_fCentr: %d\n', nF);

    % Inicializamos la matriz para almacenar los centroides de las caras
    fCentr = zeros(2, nF);
    
    % Cada cara tiene 2 vértices (según fNbVLoc)
    numVerticesPerFace = fNbVLoc;  % Esto es 2, de acuerdo con la documentación

    % Verificación del tamaño de fNbV
    totalVertices = length(fNbV);
    fprintf('Tamaño de fNbV: %d\n', totalVertices);
    
    if totalVertices == 0
        error('El vector fNbV está vacío en calc_fCentr.');
    end

    % Recorremos cada cara
    for iF = 1:nF
       
        % El índice de los vértices para la cara iF está dado por:
        startIdx = (iF - 1) * numVerticesPerFace + 1;  % Índice inicial en fNbV
        endIdx = startIdx + numVerticesPerFace - 1;    % Índice final en fNbV

        % Verificación de que los índices están dentro de los límites
        if endIdx > totalVertices
            error('Índice fuera de los límites en calc_fCentr: startIdx = %d, endIdx = %d, totalVertices = %d', startIdx, endIdx, totalVertices);
        end
        
        % Extraemos los índices de los vértices de la cara iF
        faceVertices = fNbV(startIdx:endIdx);
        
        % Obtenemos las coordenadas de estos vértices
        vCoordsForFace = vCoord(:, faceVertices);
        
        % Calculamos el centroide como el promedio de las coordenadas de los vértices
        fCentr(:, iF) = mean(vCoordsForFace, 2);
    end
end








