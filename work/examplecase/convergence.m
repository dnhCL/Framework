clear; clc; close all;

% Definir los tamaños de malla (agregamos 800x800 y 1600x1600)
Nx_values = [50, 100, 200, 400, 800, 1600];
Ny_values = [50, 100, 200, 400, 800, 1600];

% Inicializar matriz de errores
errors_temp = zeros(1, length(Nx_values)); % Errores para temperatura
errors_pres = zeros(1, length(Nx_values)); % Errores para presión

% Estructuras para almacenar perfiles
temperature_profiles_x = cell(1, length(Nx_values));  
velocity_profiles_x = cell(1, length(Nx_values));  

fprintf("\nStarting Grid Convergence Study...\n");

for i = 1:length(Nx_values)
    Nx = Nx_values(i);
    Ny = Ny_values(i);

    fprintf("Running case with Nx = %d, Ny = %d...\n", Nx, Ny);

    % 🔹 Ejecutar Solver de Temperatura
    result_temp = ssolver1(Nx, Ny);
    errors_temp(i) = result_temp.L2Error;
    

    % 🔹 Ejecutar Solver de Presión
    result_pres = psolver(Nx, Ny);
    errors_pres(i) = result_pres.rel_error;
    

    % Definir dimensiones del dominio
    Lx = 1.0;  
    Ly = 1.0;
    dx = Lx / Nx;
    dy = Ly / Ny;

    % Ubicación para extracción de perfiles
    x_target = (Lx / 2) - (dx / 2);
    tolerance = dx / 2;

    % Obtener coordenadas y valores de temperatura
    cell_centers = result_temp.dom.cCoord;
    temperatures = result_temp.T.data;

    % 🔹 Extraer Perfil de Temperatura en x = Lx/2
    cells_in_line_x = abs(cell_centers(1, :) - x_target) < tolerance;
    y_values_x = cell_centers(2, cells_in_line_x);
    temperature_profile_x = temperatures(cells_in_line_x);

    % Filtrar y ordenar datos
    valid_range_y = (y_values_x >= 0) & (y_values_x <= Ly);
    y_values_x = y_values_x(valid_range_y);
    temperature_profile_x = temperature_profile_x(valid_range_y);
    [y_values_x, sort_idx] = sort(y_values_x);
    temperature_profile_x = temperature_profile_x(sort_idx);
    
    % Guardar perfil de temperatura
    temperature_profiles_x{i} = struct('y', y_values_x, 'T', temperature_profile_x);

    % 🔹 Extraer Perfil de Velocidad en x = Lx/2
    cell_centers_pres = result_pres.dom.cCoord;
    velocity_x = result_pres.U.data(1, :);

    cells_in_line_x = abs(cell_centers_pres(1, :) - x_target) < tolerance;
    y_values_vx = cell_centers_pres(2, cells_in_line_x);
    velocity_profile_x = velocity_x(cells_in_line_x);

    % Filtrar y ordenar datos
    valid_range_y = (y_values_vx >= 0) & (y_values_vx <= Ly);
    y_values_vx = y_values_vx(valid_range_y);
    velocity_profile_x = velocity_profile_x(valid_range_y);
    [y_values_vx, sort_idx] = sort(y_values_vx);
    velocity_profile_x = velocity_profile_x(sort_idx);

    % Guardar perfil de velocidad
    velocity_profiles_x{i} = struct('y', y_values_vx, 'U', velocity_profile_x);

    % 🔹 Comparar con la malla anterior (Interpolación)
    if i > 1
        % Temperatura
        T_interp = interp1(temperature_profiles_x{i-1}.y, temperature_profiles_x{i-1}.T, y_values_x, 'linear', 'extrap');
        error_temp_profile = norm(temperature_profile_x - T_interp) / norm(T_interp);
        fprintf("Diferencia de temperatura entre %dx%d y %dx%d: %10.4e\n", Nx_values(i), Ny_values(i), Nx_values(i-1), Ny_values(i-1), error_temp_profile);

        % Velocidad
        U_interp = interp1(velocity_profiles_x{i-1}.y, velocity_profiles_x{i-1}.U, y_values_vx, 'linear', 'extrap');
        error_vel_profile = norm(velocity_profile_x - U_interp) / norm(U_interp);
        fprintf("Diferencia de velocidad entre %dx%d y %dx%d: %10.4e\n", Nx_values(i), Ny_values(i), Nx_values(i-1), Ny_values(i-1), error_vel_profile);
    end
end

% 🔹 Graficar convergencia del error de temperatura
figure;
loglog(Nx_values, errors_temp, '-o', 'LineWidth', 2);
xlabel('Number of Cells per Side');
ylabel('Relative Error (L2 Norm)');
title('Grid Convergence Study - Temperature');
grid on;

% 🔹 Graficar convergencia del error de velocidad
figure;
loglog(Nx_values, errors_pres, '-o', 'LineWidth', 2);
xlabel('Number of Cells per Side');
ylabel('Relative Error (L2 Norm)');
title('Grid Convergence Study - Velocity');
grid on;

% 🔹 Graficar diferencias entre soluciones de mallas sucesivas (Temperatura)
figure;
hold on;
for i = 2:length(Nx_values)
    if ~isempty(temperature_profiles_x{i}) && ~isempty(temperature_profiles_x{i-1})
        T_interp = interp1(temperature_profiles_x{i-1}.y, temperature_profiles_x{i-1}.T, temperature_profiles_x{i}.y, 'linear', 'extrap');
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

% 🔹 Graficar diferencias entre soluciones de mallas sucesivas (Velocidad)
figure;
hold on;
for i = 2:length(Nx_values)
    if ~isempty(velocity_profiles_x{i}) && ~isempty(velocity_profiles_x{i-1})
        U_interp = interp1(velocity_profiles_x{i-1}.y, velocity_profiles_x{i-1}.U, velocity_profiles_x{i}.y, 'linear', 'extrap');
        diff_U = velocity_profiles_x{i}.U - U_interp;
        plot(velocity_profiles_x{i}.y, diff_U, '-o', 'DisplayName', sprintf('%dx%d vs %dx%d', Nx_values(i), Ny_values(i), Nx_values(i-1), Ny_values(i-1)));
    end
end
xlabel('y');
ylabel('Velocity Difference');
title('Diferencia de velocidad entre mallas sucesivas');
legend('show');
grid on;
hold off;
