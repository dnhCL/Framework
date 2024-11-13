function fCentr = calc_fCentr(nF, fNbVLoc, fNbV, vCoord)
    % Function to calculate the centroids of faces in a mesh
    % Inputs:
    % - nF: Number of faces in the mesh
    % - fNbVLoc: Number of vertices per face (typically 2 in 2D)
    % - fNbV: Array containing vertex indices for each face
    % - vCoord: Array of vertex coordinates
    %
    % Output:
    % - fCentr: Array of face centroids (2 x nF), where each column represents the centroid of a face

    % Display the number of faces for debugging
    %%fprintf('Number of faces (nF) in calc_fCentr: %d\n', nF);

    % Initialize the centroid array
    fCentr = zeros(2, nF);

    % Confirm that the number of vertices per face is set to 2
    numVerticesPerFace = fNbVLoc;

    % Check the size of fNbV to ensure it has data
    totalVertices = length(fNbV);
    %%fprintf('Size of fNbV: %d\n', totalVertices);

    if totalVertices == 0
        error('fNbV vector is empty in calc_fCentr.');
    end

    % Loop through each face to calculate its centroid
    for iF = 1:nF
        % Calculate the indices of the vertices for face iF
        startIdx = (iF - 1) * numVerticesPerFace + 1;
        endIdx = startIdx + numVerticesPerFace - 1;

        % Ensure indices are within bounds
        if endIdx > totalVertices
            error('Index out of bounds in calc_fCentr: startIdx = %d, endIdx = %d, totalVertices = %d', startIdx, endIdx, totalVertices);
        end

        % Retrieve the vertex indices for the current face
        faceVertices = fNbV(startIdx:endIdx);

        % Get the coordinates of these vertices
        vCoordsForFace = vCoord(:, faceVertices);

        % Calculate the centroid as the average of vertex coordinates
        fCentr(:, iF) = mean(vCoordsForFace, 2);
    end
end









