function fTangent = calc_fTangent(nF, fNormal)

    % Inicializamos la matriz para almacenar los vectores tangentes
    fTangent = zeros(2, nF);

    % Recorremos cada cara
    for iF = 1:nF
        % Obtenemos el vector normal de la cara iF
        normal = fNormal(:, iF);
        
        % Calculamos el vector tangente perpendicular al vector normal
        % Si estamos en 2D, simplemente rotamos la normal 90 grados
        tangent = [-normal(2); normal(1)];
        
        % Normalizamos el vector tangente
        fTangent(:, iF) = tangent / norm(tangent);
        
        % Monitorear el vector tangente normalizado
        %%fprintf('Vector tangente normalizado de la cara %d: [%f, %f]\n', iF, fTangent(1, iF), fTangent(2, iF));
    end
end





