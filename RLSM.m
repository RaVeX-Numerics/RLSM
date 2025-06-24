%% Random Level Set Method for Logistic Free-Boundary Problems
% This is the official implementation for the paper:
% "Tracking interfaces in a random logistic free-boundary diffusion problems: 
%  a Random Level Set Method"
%
% AUTHORS:
% Vera Egorova, M. C. Casabán, Rafael Company, Lucas Jódar
%
% REQUIREMENTS:
% MATLAB R2021a or newer.
% (Optional) Parallel Computing Toolbox for faster execution with Nsim > 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;
clc;
close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Simulation Setup ---
% Choose the initial habitat geometry:
% 1: Ellipse / Circle
% 2: Square
% 3: Equilateral Triangle
% 4: Rectangle
form = 1;

% Choose the numerical solver:
% 1: EFDM (Classic Explicit FDM) - Conditionally stable
% 2: IEFDM (Improved EFDM - Larkin) - Better stability
% 3: ADE (Alternating Direction Explicit - Saul'ev)
% 4: RALADE (Random Average Larkin ADE) - Recommended for stability and accuracy
method = 4;

% --- Monte Carlo Setup ---
% Number of Monte Carlo realizations.
% Set Nsim = 1 for a single, deterministic run.
Nsim = 100; % Options: 1, 25, 50, 100, 200, 400, 800

% --- Output Control ---
% Set to 1 to save figures to a 'results' subdirectory, 0 to disable.
save_flag = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Load or Define Random Parameters ---
if Nsim > 1
    % For stochastic runs, load pre-generated random variables.
    % NOTE: Ensure the .mat file corresponding to Nsim exists.
    % You may need to generate your own random data for other Nsim values.
    try
        load(['Dw_K', num2str(Nsim), '.mat']);
        load(['muw_K', num2str(Nsim), '.mat']);
        fprintf('Loaded random data for Nsim = %d.\n', Nsim);
    catch
        error('Random data file for Nsim = %d not found. Please generate it or choose a different Nsim.', Nsim);
    end
else
    % For deterministic run, define constant parameters.
    D = 1;  % Diffusion coefficient
    mu = 2; % Front expansion rate
    fprintf('Running a single deterministic simulation (Nsim = 1).\n');
end

% --- Simulation and Grid Parameters ---
T = 3;          % Total simulation time
Nx = 200; Ny = Nx; % Number of grid points
xmin = -10; xmax = 10;
ymin = -10; ymax = 10;

% Create grid
x = linspace(xmin, xmax, Nx + 1);
y = linspace(ymin, ymax, Ny + 1);
[X, Y] = meshgrid(x, y);
hx = (xmax - xmin) / Nx;
hy = (ymax - ymin) / Ny;
h = min(hx, hy);
c2 = (hx / hy)^2;

% --- Model Parameters (Logistic Growth) ---
% alpha: intrinsic growth rate
% beta: carrying capacity parameter
alpha = ones(size(X));
beta = ones(size(Y));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GEOMETRY AND INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch form
    case 1 % Ellipse / Circle
        l = 3; % Base radius for the circle case
        b = l; % Set b = l for a circle, or b != l for an ellipse
        a = l^2 / b;

        if b == l
            form_title = 'circle';
        else
            form_title = 'ellipse';
        end

        % Centered initial population
        U_initial = (X.^2 / a^2 + Y.^2 / b^2 < 1) .* (1 - (X.^2 / a^2 + Y.^2 / b^2));
        phi_initial = (sqrt(X.^2 / a^2 + Y.^2 / b^2) - 1);

    case 2 % Square
        a = 2.13; % Side length
        phi_initial = max(abs(X) - a, abs(Y) - a);
        U_initial = (X.^2 - a^2) .* (Y.^2 - a^2);
        U_initial(phi_initial > 0) = 0;
        form_title = 'square';

    case 3 % Equilateral Triangle
        a = 2.2 * sqrt(pi * 4 / sqrt(3)); % Side length to match area of circle with r=2.2
        ht = a * sqrt(3) / 2; % Triangle height
        v1 = [-a/2, -ht/3]; v2 = [a/2, -ht/3]; v3 = [0, 2*ht/3]; % Vertices

        % Barycentric coordinates
        lambda1 = ((Y-v2(2)).*(v3(1)-v2(1)) - (X-v2(1)).*(v3(2)-v2(2))) ./ ((v1(2)-v2(2)).*(v3(1)-v2(1)) - (v1(1)-v2(1)).*(v3(2)-v2(2)));
        lambda2 = ((Y-v3(2)).*(v1(1)-v3(1)) - (X-v3(1)).*(v1(2)-v3(2))) ./ ((v2(2)-v3(2)).*(v1(1)-v3(1)) - (v2(1)-v3(1)).*(v1(2)-v3(2)));
        lambda3 = 1 - lambda1 - lambda2;

        phi_initial = -min(cat(3, lambda1, lambda2, lambda3), [], 3);
        U_initial = lambda1 .* lambda2 .* lambda3;
        U_initial(phi_initial > 0) = 0;
        form_title = 'triangle';

    case 4 % Rectangle
        l = 2.2 * sqrt(pi); % Equivalent length to match area
        a = 1.3648; % Width
        b = l^2 / (4 * a); % Height
        
        phi_initial = max(abs(X) - a, abs(Y) - b);
        U_initial = (X.^2 - a^2) .* (Y.^2 - b^2);
        U_initial(phi_initial > 0) = 0;
        if a == b, form_title = 'square'; else, form_title = 'rectangle'; end
end

% Normalize initial population if desired (currently commented out)
% Population_initial = trapz(y,trapz(x,U_initial,2));
% U_initial = U_initial / Population_initial;
Population_initial = trapz(y, trapz(x, U_initial, 2));
habitat_area_initial = trapz(y, trapz(x, phi_initial < 0, 2));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  STABILITY ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate stable time step 'dt' based on the chosen method
D_max = max(D);
C0 = max(max(alpha ./ beta));
M0 = max(max(U_initial));
Umax = max(C0, M0);
kappa1 = min([min(alpha), min(min(beta))]);
kappa2 = max([max(alpha), max(max(beta))]);

switch method
    case 1 % EFDM
        Q1 = kappa1 * (2 * kappa2^2 / kappa1^2 - 1);
        Q2 = kappa2 * (2 * Umax - kappa1 / kappa2);
        cond = hx^2 / (2 * D_max * (1 + c2) + hx^2 * max(Q1, Q2));
    case {2, 3, 4} % IEFDM, ADE, RALADE
        cota1 = (kappa2 * Umax - kappa1);
        cota2 = (2 * kappa2 * Umax - kappa1);
        cond1 = hx^2 / (D_max * (1 + c2) + hx^2 * cota1);
        cond2 = hx^2 / (D_max * (1 + c2) + hx^2 * cota2);
        cond = min(cond1, cond2);
        if method > 2, cond = 2 * cond; end % ADE methods allow larger step
end

dt = min(0.5 * h^2, cond);
t_vector = 0:dt:T;
N_steps = length(t_vector) - 1;

fprintf('Stability analysis complete. Time step dt = %.4e\n', dt);
fprintf('Parabolic mesh ratio dt/h^2 = %.4f\n', dt / h^2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  MAIN SIMULATION LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Starting Monte Carlo simulation with %d realizations...\n', Nsim);

% --- Initialize storage for statistics ---
U_mean = zeros(size(U_initial));
U_squared_sum = zeros(size(U_initial));
phi_mean = zeros(size(phi_initial));
Population_mean = zeros(1, N_steps + 1);
Population_squared_sum = zeros(1, N_steps + 1);

tic; % Start timer

% Use parfor for parallel execution if Parallel Computing Toolbox is available
% for w = 1:Nsim
parfor w = 1:Nsim
    % --- Per-realization setup ---
    Population_w = zeros(1, N_steps + 1);
    Population_w(1) = Population_initial;
    DD = D(w);
    muu = mu(w);

    % Pre-calculate constants for efficiency
    if method > 2, W = 2 * hx^2 / dt; else, W = hx^2 / dt; end
    Cn = (W - DD * (1 + c2));
    Cn1 = (W + DD * (1 + c2));

    % Initialize solution variables for this realization
    U = U_initial;
    phi = phi_initial;
    V = U; Q = U; % Auxiliary variables for ADE methods

    % --- Time Integration Loop ---
    for n = 1:N_steps
        % 1. Compute gradient of the level set function (central differences)
        dphidy = (phi(3:end, 2:end-1) - phi(1:end-2, 2:end-1)) / (2 * hy);
        dphidx = (phi(2:end-1, 3:end) - phi(2:end-1, 1:end-2)) / (2 * hx);

        % 2. Compute gradient of population U to find advection velocity
        dudy = (U(3:end, 2:end-1) - U(1:end-2, 2:end-1)) / (2 * hy);
        dudx = (U(2:end-1, 3:end) - U(2:end-1, 1:end-2)) / (2 * hx);
        
        % 3. Update the level set function
        vdotn = muu .* sqrt(dudx.^2 + dudy.^2) .* sqrt(dphidx.^2 + dphidy.^2);
        phi(2:end-1, 2:end-1) = phi(2:end-1, 2:end-1) - dt * vdotn;

        % 4. Update the solution U using the chosen numerical method
        switch method
            case 1 % EFDM
                LapU = (U(3:end, 2:end-1) - 2*U(2:end-1, 2:end-1) + U(1:end-2, 2:end-1))/hy^2 ...
                     + (U(2:end-1, 3:end) - 2*U(2:end-1, 2:end-1) + U(2:end-1, 1:end-2))/hx^2;
                U(2:end-1, 2:end-1) = U(2:end-1, 2:end-1) + dt * DD * LapU ...
                    + dt * U(2:end-1, 2:end-1) .* (alpha(2:end-1, 2:end-1) - beta(2:end-1, 2:end-1) .* U(2:end-1, 2:end-1));
            
            case 2 % Improved EFDM (Larkin)
                U_new = U;
                for i = 2:Ny, for j = 2:Nx
                    U_new(i,j) = (Cn*U(i,j) + DD*(U(i,j+1) + U(i,j-1)) + DD*c2*(U(i+1,j) + U(i-1,j)) ...
                                + hx^2*U(i,j)*(alpha(i,j)-beta(i,j)*U(i,j))) / Cn1;
                end, end
                U = U_new;
                
            case 4 % RALADE (Larkin's ADE)
                % Forward sweep
                for i = 2:Ny, for j = 2:Nx
                    V(i,j) = (Cn*U(i,j) + DD*(U(i,j+1) + V(i,j-1)) + DD*c2*(U(i+1,j) + V(i-1,j)) ...
                            + hx^2*U(i,j)*(alpha(i,j)-beta(i,j)*U(i,j))) / Cn1;
                end, end
                % Backward sweep
                for i = Ny:-1:2, for j = Nx:-1:2
                    Q(i,j) = (Cn*U(i,j) + DD*(Q(i,j+1) + U(i,j-1)) + DD*c2*(Q(i+1,j) + U(i-1,j)) ...
                            + hx^2*U(i,j)*(alpha(i,j)-beta(i,j)*U(i,j))) / Cn1;
                end, end
                U = (V + Q) / 2; % Average the sweeps
        end

        % 5. Apply free boundary condition (population is zero outside the habitat)
        U(phi >= 0) = 0;

        % 6. Enforce Neumann boundary conditions on the level-set function
        phi(1, :) = phi(2, :); phi(end, :) = phi(end-1, :);
        phi(:, 1) = phi(:, 2); phi(:, end) = phi(:, end-1);

        % 7. Store total population for this realization
        Population_w(n + 1) = trapz(y, trapz(x, U, 2));
    end
    
        U_mean = U_mean + U;
        U_squared_sum = U_squared_sum + U.^2;
        phi_mean = phi_mean + phi;
        Population_mean = Population_mean + Population_w;
        Population_squared_sum = Population_squared_sum + Population_w.^2;
end
toc; % End timer
fprintf('Simulation finished. Total time: %.2f seconds.\n', toc);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  POST-PROCESSING AND STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calculating statistics and generating plots...\n');

% Calculate Mean and Standard Deviation
U_mean = U_mean / Nsim;
U_std = sqrt(U_squared_sum / Nsim - U_mean.^2);

phi_mean = phi_mean / Nsim;

Population_mean = Population_mean / Nsim;
Population_std = sqrt(Population_squared_sum / Nsim - Population_mean.^2);

% --- Display Final Statistics ---
Population_final = Population_mean(end);
habitat_area_final = trapz(y, trapz(x, phi_mean < 0, 2));

fprintf('\n--- Simulation Results ---\n');
fprintf('Initial Population: P(0) = %.4f\n', Population_initial);
fprintf('Final Mean Population: E[P(T)] = %.4f\n', Population_final);
fprintf('Initial Habitat Area: A(0) = %.4f\n', habitat_area_initial);
fprintf('Final Mean Habitat Area: E[A(T)] = %.4f\n', habitat_area_final);
fprintf('------------------------\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  VISUALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Define file saving parameters ---
if save_flag
    results_dir = 'results';
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    filename_base = ['Random_', form_title, '_T', num2str(T), '_K', num2str(Nsim), '_M', num2str(method)];
end

% --- Figure 1: Initial Population ---
figure('Name', 'Initial Conditions', 'NumberTitle', 'off');
surf(X, Y, U_initial, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
contour(X, Y, phi_initial, [0 0], 'b-', 'LineWidth', 2);
xlabel('x'); ylabel('y'); zlabel('Population Density');
title('Initial Population and Habitat');
grid on; shading interp; view(3); colorbar;

if save_flag, saveas(gcf, fullfile(results_dir, [filename_base, '_initial.png'])); end

% --- Figure 2: Front Evolution ---
figure('Name', 'Front Evolution', 'NumberTitle', 'off');
contour(X, Y, phi_initial, [0 0], 'b--', 'LineWidth', 2);
hold on;
contour(X, Y, phi_mean, [0 0], 'r-', 'LineWidth', 2);
xlabel('x'); ylabel('y');
title('Habitat Front Evolution');
legend('Initial Front (t=0)', sprintf('Mean Final Front (t=%.1f)', T));
axis equal; grid on;

if save_flag, saveas(gcf, fullfile(results_dir, [filename_base, '_front.png'])); end

% --- Figure 3: Mean Final Population ---
figure('Name', 'Mean Final Population', 'NumberTitle', 'off');
surf(X, Y, U_mean, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
hold on;
contour(X, Y, phi_mean, [0 0], 'k-', 'LineWidth', 2);
xlabel('x'); ylabel('y'); zlabel('Mean Population Density');
title(sprintf('Mean Population at t=%.1f (Nsim=%d)', T, Nsim));
shading interp;  colorbar; view(3);

if save_flag, saveas(gcf, fullfile(results_dir, [filename_base, '_final_mean.png'])); end

% --- Figure 4: Population Dynamics (P(t)) ---
figure('Name', 'Population Dynamics P(t)', 'NumberTitle', 'off');
plot(t_vector, Population_mean, 'r-', 'LineWidth', 2.5);
hold on;
plot(t_vector, Population_mean + Population_std, 'k--');
plot(t_vector, Population_mean - Population_std, 'k--');
grid on;
xlabel('Time (t)');
ylabel('Total Population P(t)');
title('Population Dynamics Over Time');
legend('Mean Population E[P(t)]', 'Mean \pm Std. Dev.');

if save_flag, saveas(gcf, fullfile(results_dir, [filename_base, '_population.png'])); end
