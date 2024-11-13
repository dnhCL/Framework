function fNormal = calc_fNormal(fNbVLoc, fNbV, vCoord, fArea)
    % Function to calculate the unit normal vectors of faces in a mesh
    % Inputs:
    % - fNbVLoc: Number of vertices per face (typically 2 in 2D)
    % - fNbV: Array containing vertex indices for each face
    % - vCoord: Array of vertex coordinates
    % - fArea: Array of face lengths (calculated in calc_fArea)
    %
    % Output:
    % - fNormal: Array of unit normal vectors (2 x nF), where each column represents the normal of a face

    % Calculate the number of faces in the mesh
    nF = length(fNbV) / fNbVLoc;

    % Initialize the array to store normal vectors
    fNormal = zeros(2, nF);

    % Loop through each face to calculate its normal vector
    for iF = 1:nF
        % Determine the indices of the vertices for face iF
        startIdx = (iF - 1) * fNbVLoc + 1;
        endIdx = startIdx + fNbVLoc - 1;

        % Retrieve the vertex indices for the current face
        faceVertices = fNbV(startIdx:endIdx);
        v1 = vCoord(:, faceVertices(1));
        v2 = vCoord(:, faceVertices(2));

        % Calculate the tangent vector (from v1 to v2)
        tangent = v2 - v1;

        % Rotate the tangent vector by 90 degrees to get the normal vector
        normal = [-tangent(2); tangent(1)];

        % Normalize the normal vector by dividing by the face length (fArea)
        fNormal(:, iF) = normal / fArea(iF);
    end
end






