function fXi = calc_fXi(nF, cCoord, fNbCLoc, fNbC)
    % Function to calculate the Xi vectors between cell centers across each face
    % Inputs:
    % - nF: Number of faces in the mesh
    % - cCoord: Array of cell centroids (2 x nC)
    % - fNbCLoc: Number of cells per face (typically 2 in 2D)
    % - fNbC: Array of neighboring cells for each face
    %
    % Output:
    % - fXi: Array of Xi vectors (2 x nF), where each column represents the vector from the owner cell to the neighbor cell for each face

    % Initialize the array to store Xi vectors
    fXi = zeros(2, nF);

    % Loop through each face to calculate the Xi vector
    for iF = 1:nF
        % Calculate the indices of the neighboring cells for face iF
        ownerIdx = fNbCLoc * (iF - 1) + 1;
        neighborIdx = ownerIdx + 1;

        % Retrieve the cell indices
        c1 = fNbC(ownerIdx); % Owner cell
        c2 = fNbC(neighborIdx); % Neighbor cell

        % Retrieve the coordinates of the cell centers
        c1Coord = cCoord(:, c1);  % Center coordinates of cell c1
        c2Coord = cCoord(:, c2);  % Center coordinates of cell c2
        
        % Calculate the Xi vector from c1 to c2
        xiVector = c2Coord - c1Coord;
        
        % Assign the Xi vector to the result
        fXi(:, iF) = xiVector;
    end
end







