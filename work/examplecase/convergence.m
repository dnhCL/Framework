clear; clc; close all;

% Definir los tamaños de malla (agregamos 800x800 y 1600x1600)
Nx_values = [50, 100, 200, 400, 800, 1600];
Ny_values = [50, 100, 200, 400, 800, 1600];

% Inicializar matriz de errores
errors = zeros(1, length(Nx_values));

% Estructura para almacenar perfiles de temperatura
temperature_profiles_x = cell(1, length(Nx_values));  

fprintf("\nStarting Grid Convergence Study...\n");

for i = 1:length(Nx_values)
    Nx = Nx_values(i);
    Ny = Ny_values(i);

    fprintf("Running case with Nx = %d, Ny = %d...\n", Nx, Ny);
    result = ssolver1(Nx, Ny);

    % Almacenar error L2
    errors(i) = result.L2Error;
    fprintf("Relative L2 Error: %10.4e\n", errors(i));

    % Definir dimensiones del dominio
    Lx = 1.0;  
    Ly = 1.0;
    dx = Lx / Nx;
    dy = Ly / Ny;

    % Ubicación para extracción de perfiles
    x_target = (Lx / 2) - (dx / 2);
    y_target = (Ly / 2) - (dy / 2);
    tolerance = dx / 2;

    % Obtener coordenadas y temperaturas
    cell_centers = result.dom.cCoord;
    temperatures = result.T.data;

    % 🔹 Extraer Perfil Vertical en x = Lx/2
    cells_in_line_x = abs(cell_centers(1, :) - x_target) < tolerance;
    y_values_x = cell_centers(2, cells_in_line_x);
    temperature_profile_x = temperatures(cells_in_line_x);

    % Filtrar y ordenar datos
    valid_range_y = (y_values_x >= 0) & (y_values_x <= Ly);
    y_values_x = y_values_x(valid_range_y);
    temperature_profile_x = temperature_profile_x(valid_range_y);
    [y_values_x, sort_idx] = sort(y_values_x);
    temperature_profile_x = temperature_profile_x(sort_idx);
    
    % Guardar el perfil
    temperature_profiles_x{i} = struct('y', y_values_x, 'T', temperature_profile_x);

    % 🔹 Comparar con la malla anterior (Interpolación)
    if i > 1
        % Interpolamos la solución de la malla anterior a la malla actual
        T_interp = interp1(temperature_profiles_x{i-1}.y, temperature_profiles_x{i-1}.T, y_values_x, 'linear', 'extrap');

        % Calcular error de perfil
        error_profile = norm(temperature_profile_x - T_interp) / norm(T_interp);
        fprintf("Diferencia de temperatura entre %dx%d y %dx%d: %10.4e\n", Nx_values(i), Ny_values(i), Nx_values(i-1), Ny_values(i-1), error_profile);
    end
end

% 🔹 Graficar convergencia del error
figure;
loglog(Nx_values, errors, '-o', 'LineWidth', 2);
xlabel('Number of Cells per Side');
ylabel('Relative Error (L2 Norm)');
title('Grid Convergence Study');
grid on;

% 🔹 Graficar diferencias entre soluciones de mallas sucesivas
figure;
hold on;
for i = 2:length(Nx_values)
    if ~isempty(temperature_profiles_x{i}) && ~isempty(temperature_profiles_x{i-1})
        % Interpolación de la malla anterior a la malla actual
        T_interp = interp1(temperature_profiles_x{i-1}.y, temperature_profiles_x{i-1}.T, temperature_profiles_x{i}.y, 'linear', 'extrap');

        % Diferencia entre soluciones interpoladas
        diff_T = temperature_profiles_x{i}.T - T_interp;
        plot(temperature_profiles_x{i}.y, diff_T, '-o', 'DisplayName', sprintf('%dx%d vs %dx%d', Nx_values(i), Ny_values(i), Nx_values(i-1), Ny_values(i-1)));
    end
end
xlabel('y');
ylabel('Temperature Difference');
title('Diferencia de temperatura entre mallas sucesivas');
legend('show');
grid on;
hold off;
