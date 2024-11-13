function fXiMag = calc_fXiMag(nF, fXi)
   % Function to calculate the magnitudes of Xi vectors for each face
   % Inputs:
   % - nF: Number of faces in the mesh
   % - fXi: Array of Xi vectors (2 x nF), where each column represents the vector from the owner cell to the neighbor cell for each face
   %
   % Output:
   % - fXiMag: Array of magnitudes of Xi vectors (1 x nF)

   % Initialize the array to store Xi magnitudes
   fXiMag = zeros(1, nF);

   % Loop through each face to calculate the magnitude of Xi
   for jF = 1:nF
       % Calculate the magnitude of the Xi vector for face jF
       fXiMag(jF) = sqrt(fXi(1, jF)^2 + fXi(2, jF)^2);
   end
end




