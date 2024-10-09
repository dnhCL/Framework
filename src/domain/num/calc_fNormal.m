function fNormal = calc_fNormal(fNbVLoc, fNbV, vCoord, fArea)

    % Número de caras nF dado por el tamaño de fNbV dividido por el número de vértices por cara (2)
    nF = length(fNbV) / fNbVLoc;
    fprintf('Número de caras (nF) en calc_fNormal: %d\n', nF);

    % Inicializamos la matriz para almacenar los vectores normales
    fNormal = zeros(2, nF);

    % Monitorear las primeras coordenadas de vCoord
    fprintf('Primeras 10 coordenadas de vCoord:\n');
    disp(vCoord(:, 1:min(10, size(vCoord, 2))));

    % Monitorear las primeras áreas de fArea
    fprintf('Primeras 10 áreas de fArea:\n');
    disp(fArea(1:min(10, length(fArea))));

    % Recorremos cada cara
    for iF = 1:nF
        fprintf('Procesando cara %d de %d\n', iF, nF);
        
        % Los índices de los vértices que definen la cara están en fNbV
        startIdx = (iF - 1) * fNbVLoc + 1;  % Índice inicial en fNbV
        endIdx = startIdx + fNbVLoc - 1;    % Índice final en fNbV

        % Verificación de los índices de los vértices
        fprintf('Índices de los vértices de la cara %d: [%d, %d]\n', iF, startIdx, endIdx);

        % Extraemos los dos vértices que forman la cara
        faceVertices = fNbV(startIdx:endIdx);
        v1 = vCoord(:, faceVertices(1));
        v2 = vCoord(:, faceVertices(2));

        % Monitorear las coordenadas de los vértices
        fprintf('Coordenadas del vértice 1 (v1) de la cara %d: [%f, %f]\n', iF, v1(1), v1(2));
        fprintf('Coordenadas del vértice 2 (v2) de la cara %d: [%f, %f]\n', iF, v2(1), v2(2));

        % Calculamos el vector tangente (de v1 a v2)
        tangent = v2 - v1;

        % El vector normal se obtiene rotando el vector tangente 90 grados
        normal = [-tangent(2); tangent(1)];

        % Monitorear el vector normal antes de la normalización
        fprintf('Vector normal no normalizado de la cara %d: [%f, %f]\n', iF, normal(1), normal(2));

        % Normalizamos el vector normal dividiendo por el área (longitud de la cara)
        fNormal(:, iF) = normal / fArea(iF);

        % Monitorear el vector normal normalizado
        fprintf('Vector normal normalizado de la cara %d: [%f, %f]\n', iF, fNormal(1, iF), fNormal(2, iF));
    end
end





