function fTangent = calc_fTangent(nF, fNormal)
    % Function to calculate the tangent vectors of faces in a mesh
    % Inputs:
    % - nF: Number of faces in the mesh
    % - fNormal: Array of normal vectors (2 x nF), where each column represents the normal of a face
    %
    % Output:
    % - fTangent: Array of tangent vectors (2 x nF), where each column represents the tangent of a face

    % Initialize the array to store tangent vectors
    fTangent = zeros(2, nF);

    % Loop through each face to calculate its tangent vector
    for iF = 1:nF
        % Retrieve the normal vector for face iF
        normal = fNormal(:, iF);
        
        % Calculate the tangent vector as a 90-degree rotation of the normal
        % Since we are in 2D, we can rotate the normal vector by 90 degrees
        tangent = [-normal(2); normal(1)];
        
        % Normalize the tangent vector
        fTangent(:, iF) = tangent / norm(tangent);
        
        % Monitor the normalized tangent vector
        %%fprintf('Normalized tangent vector for face %d: [%f, %f]\n', iF, fTangent(1, iF), fTangent(2, iF));
    end
end





