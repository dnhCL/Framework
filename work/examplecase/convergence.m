clear; clc; close all;

% Definir los tamaños de malla
Nx_values = [50, 100, 200, 400];
Ny_values = [50, 100, 200, 400];

% Inicializar matriz de errores
errors = zeros(1, length(Nx_values));

fprintf("\nStarting Grid Convergence Study...\n");

% Ejecutar cada caso con diferentes tamaños de malla
for i = 1:length(Nx_values)
    Nx = Nx_values(i);
    Ny = Ny_values(i);

    fprintf("Running case with Nx = %d, Ny = %d...\n", Nx, Ny);
    result = ssolver1(Nx, Ny);

    % Almacenar error L2 para el análisis de convergencia
    errors(i) = result.L2Error;
    fprintf("Relative L2 Error: %10.4e\n", errors(i));
end

% Cálculo del Orden de Convergencia
r = 2; % Factor de refinamiento
p = log(errors(2)/errors(3)) / log(r);
fprintf("\nEstimated Order of Accuracy: p = %.2f\n", p);

% Grid Convergence Index (GCI)
Fs = 1.25;
GCI = Fs * abs(errors(3) / errors(2)) / (r^p - 1);
fprintf("GCI for finest grid: %.4f\n", GCI);

% Graficar convergencia
figure;
loglog(Nx_values, errors, '-o', 'LineWidth', 2);
xlabel('Number of Cells per Side');
ylabel('Relative Error (L2 Norm)');
title('Grid Convergence Study');
grid on;
