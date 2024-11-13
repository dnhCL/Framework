function fArea = calc_fArea(nF, fNbVLoc, fNbV, vCoord)
    % Function to calculate the areas (or lengths) of faces in a mesh
    % Inputs:
    % - nF: Number of faces in the mesh
    % - fNbVLoc: Number of vertices per face (typically 2 in 2D)
    % - fNbV: Array containing vertex indices for each face
    % - vCoord: Array of vertex coordinates
    %
    % Output:
    % - fArea: Array of face areas (1 x nF), where each entry represents the length of a face

    % Display the number of faces for debugging
    %fprintf('Number of faces (nF) in calc_fArea: %d\n', nF);

    % Initialize the area (or length) array
    fArea = zeros(1, nF);

    % Confirm that the number of vertices per face is set to 2
    numVerticesPerFace = fNbVLoc;

    % Check the size of fNbV to ensure it has data
    totalVertices = length(fNbV);
    %fprintf('Size of fNbV in calc_fArea: %d\n', totalVertices);

    if totalVertices == 0
        error('fNbV vector is empty in calc_fArea.');
    end

    % Loop through each face to calculate its length
    for iF = 1:nF
        % Calculate the indices of the vertices for face iF
        startIdx = (iF - 1) * numVerticesPerFace + 1;
        endIdx = startIdx + numVerticesPerFace - 1;

        % Ensure indices are within bounds
        if endIdx > totalVertices
            error('Index out of bounds in calc_fArea: startIdx = %d, endIdx = %d, totalVertices = %d', startIdx, endIdx, totalVertices);
        end

        % Retrieve the vertex indices for the current face
        faceVertices = fNbV(startIdx:endIdx);

        % Get the coordinates of the two vertices
        v1 = vCoord(:, faceVertices(1));
        v2 = vCoord(:, faceVertices(2));

        % Calculate the length of the face as the distance between the two vertices
        fArea(iF) = sqrt((v2(1) - v1(1))^2 + (v2(2) - v1(2))^2);
    end
end








