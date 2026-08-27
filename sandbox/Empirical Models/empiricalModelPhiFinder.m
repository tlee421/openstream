% --empiricalModelTest.m--%
%
%
% The purpose of this script is to evaluate the first iteration of the
% obstruction heat transfer model solver against wall temperature data for
% the results section of Tyler's ASME conference presentation.%
%
% As of mid summer 2026 the wall heat transfer model can be solved
% independently of OpenSTREAM by using the Wake width vs. Boiling number
% curve to compute the location of the dry wake. The local state pre
% obstruction is required to compute the dry wake location downstream of
% the obs. The HT model can then be solved knowing the location of the dry
% wake downstream of the obstruction. 
%
% This model is not meant to be used extensivley but is more of a quick
% workaround to quickly compare model results to data without needing to
% integrate the obstruction branch into openstream-database
%
% CHANGE LOG:
% - Added per-row optimization of phi (previously a fixed guess of 0.6).
%   phi is now fit for each test run by minimizing the SSE between
%   computed and measured thermocouple temperatures at the 5 target
%   downstream locations, using fminbnd. See fitPhiForRow / phiResidual
%   local functions at the bottom of this file.

import Solvers.Obstruction.WallHT.ftoglass_ann_optim

% Load in data regarding known values:
% - pre obs hydrodynamics
% - external wall temperatures
varNames = {'Param', 'T1', 'T2', 'T3', 'T4', 'T5', 'Tsat', 'G', 'x_in', ...
            'q_eff', 'Q_heater', 'delta', 'Re_f', 'SWR', 'mdot_f', 'tau_w', 'Re_v'};

T = readtable("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\07_Conferences\ASME SHTC\Data\WallHTModelComparison.xlsx", ...
    'VariableNamingRule', 'preserve');   % preserve original headers first
T.Properties.VariableNames = varNames;   % then overwrite with yours

%Compute the hydrodynamic wake width for each case:
T.W_adiabatic = (1.66 * T.Re_f.^0.085 .* T.SWR.^0.22)*2/1000;               %Closure relation for large obstruction wake width in [m]

%load in fluid properties
load("testfluidprops.mat")
k_film = fluid.KF;                                                         % [W/m-k] thermal conductivity of saturated liquid
k_vapor = fluid.KG;                                                        % [W/m-k] thermal conductivity of saturated vapor
rho_f = fluid.RHOF;
mu_f = fluid.MUF;
nu_f = mu_f/rho_f;
PR_f = fluid.PRANDTLF;
PR_v = fluid.PRANDTLG;
h_fg = fluid.HFG;

%Wall geometry:
W = 0.036;
H = 0.012;
per =2*W + 2*H;
A_c = W*H;
D_h = 4*A_c/per;
L = 0.525;
L_obs = 0.42; 
D_obs = 0.0025;

%Compute base film thickness using law of the wall
yplus = 5;                                                                 %Assumed yplus value base film occurs at
T.u_star = (T.tau_w/rho_f).^(0.5);
T.base_thick = yplus*(T.u_star).^(-1)*nu_f;                                %Base film thickness

%Compute wet wall HTC:
T.htc_wet = k_film * (T.base_thick).^(-1);

%Compute dry wall HTC using Dittus Boetler
T.htc_dry = k_vapor/D_h * 0.023 * (T.Re_v).^0.8 * PR_v^0.4;                %Dittus boetler

%Compute the dry wake width at N locations downstream
nodes = 22;                                                                % Number of axial nodes downstream of the obstruction. 
x_obs = L- L_obs;                                                          % Axial distance behind the obstruction
x = linspace(0,x_obs, nodes);                                              % vector of x locations behind obs

%Target downstream locations matching thermocouple positions
targets = [0.04 0.05 0.06 0.07 0.08];
tol = 1e-9;
meas_mask = any(abs(x - targets') < tol, 1);
meas_idx = find(meas_mask);

%Wall HT solver geometry/params (needed by phi fit as well as full solve)
M= 125;                                                                    %spanwise nodes
N = 125;                                                                   % wall normal nodes
W_HT = W/2;                                                                % half width of wall
H_wall = 0.003;                                                            % thickness of wall
k_wall = 1.2;                                                              % conductivity of wall 
htc_amb = 10;
T_inf = fluid.TSAT;
T_amb = 22.5 + 273.15;

%% Fit phi per test run
% Instead of assuming a fixed phi = 0.6 for every row, solve for the phi
% that minimizes the SSE between computed and measured thermocouple
% temperatures at the 5 target locations for each row independently.
phi_bounds = [0.01, 1];   % adjust if you have tighter physical bounds on phi
T.phi = nan(height(T), 1);

for i = 1:height(T)
    T.phi(i) = fitPhiForRow(i, T, targets, D_obs, per, h_fg, M, N, W, W_HT, ...
                            H_wall, k_wall, x_obs, htc_amb, T_inf, T_amb, phi_bounds);
end

%% Compute the dry wake width at N locations downstream using fitted phi
T.mdot_wake = T.mdot_f * (D_obs/per) .* T.phi;                                % 

scalarPart = (T.q_eff .* T.W_adiabatic) ./ (T.mdot_wake .* h_fg);  % N×1
Bo = scalarPart .* x;   % N×1 .* 1×M -> N×M via implicit expansion
T.Bo = num2cell(Bo, 2);  % store each row's array back into the table

wake_ratio_func = -5 + 0.0058*Bo + (6.3325 * Bo.^1.6617)./(0.1719^1.6617 + Bo.^1.6617); %Empirical curve fit for wake ratio as a function of Bo_obs
wake_ratio = max(0, wake_ratio_func); %limit wake ratio to a minimum value of 0
T.wake_ratio = num2cell(wake_ratio,2); 
dry_wake_W = T.W_adiabatic .* wake_ratio;
T.dry_wake_W = num2cell(dry_wake_W,2); 
w_free = W - dry_wake_W;                                                   %free stream width
T.w_free = num2cell(w_free,2);

% Solve wall HT solution at every downstream location:
N_rows = height(T);

% Pre-allocate cell arrays: N_rows x N_wfree (assuming each row has 100 w_free values)
T.x_full     = cell(N_rows, 1);
T.y          = cell(N_rows, 1);
T.T_full     = cell(N_rows, 1);
T.T_center   = cell(N_rows, 1);
T.q_bot_full = cell(N_rows, 1);
T.q_bot_total= cell(N_rows, 1);
T.q_top_full = cell(N_rows, 1);
T.q_top_total= cell(N_rows, 1);
T.x_out      = cell(N_rows, 1);   % renamed to avoid clashing with T.x
T.M2         = cell(N_rows, 1);

for i = 1:N_rows
    w_free_i = T.w_free(i,:);        % this row's 1x100 vector
    w_free_i = w_free_i{1};
    N_wfree = numel(w_free_i);

    % temp storage for this row's results across all w_free values
    x_full_row = cell(1, N_wfree);
    y_row = cell(1, N_wfree);
    T_full_row = cell(1, N_wfree);
    T_center_row = cell(1, N_wfree);
    q_bot_full_row = cell(1, N_wfree);
    q_bot_total_row = cell(1, N_wfree);
    q_top_full_row = cell(1, N_wfree);
    q_top_total_row = cell(1, N_wfree);
    x_out_row = cell(1, N_wfree);
    M2_row = cell(1, N_wfree);

    for j = 1:N_wfree
        [x_full_row{j}, y_row{j}, T_full_row{j}, T_center_row{j}, ...
         q_bot_full_row{j}, q_bot_total_row{j}, q_top_full_row{j}, q_top_total_row{j}, ...
         x_out_row{j}, M2_row{j}] = ftoglass_ann_optim( ...
            M, N, W_HT, H_wall, k_wall, x_obs, ...
            T.q_eff(i), w_free_i(j), T.htc_dry(i), T.htc_wet(i), ...
            htc_amb, T_inf, T_amb);
    end

    % store this row's full set of results into the table cell
    T.x_full{i}     = x_full_row;
    T.y{i}          = y_row;
    T.T_full{i}     = T_full_row;
    T.T_center{i}   = T_center_row;
    T.q_bot_full{i} = q_bot_full_row;
    T.q_bot_total{i}= q_bot_total_row;
    T.q_top_full{i} = q_top_full_row;
    T.q_top_total{i}= q_top_total_row;
    T.x_out{i}      = x_out_row;
    T.M2{i}         = M2_row;
end

%Fill computed values into T table
T.T_calc = cell(height(T), 1);  % preallocate new column
for i = 1:height(T)
    T.T_calc{i} = T.T_center{i}(meas_idx);  % pull matched values from that row's array
end

% Compare computed temperatures to measured temperatures 
Temps = table();
Temps.ParamName = T.Param;
Temps.T1_meas = T.T1;
Temps.T2_meas = T.T2;
Temps.T3_meas = T.T3;
Temps.T4_meas = T.T4;
Temps.T5_meas = T.T5;

Temps.T1_calc = cellfun(@(c) c(1), T.T_calc);
Temps.T2_calc = cellfun(@(c) c(2), T.T_calc);
Temps.T3_calc = cellfun(@(c) c(3), T.T_calc);
Temps.T4_calc = cellfun(@(c) c(4), T.T_calc);
Temps.T5_calc = cellfun(@(c) c(5), T.T_calc);

%%
%Plotting functions (AI)
% Loop through every column in the table
for colIdx = 1:width(Temps)
    colName = Temps.Properties.VariableNames{colIdx};
    
    % Keep the parameter names as text, but clean everything else
    if strcmp(colName, 'ParamName')
        if ~iscellstr(Temps.(colName)) && ~isstring(Temps.(colName))
            Temps.(colName) = string(Temps.(colName));
        end
    else
        % Force convert the temperature columns to standard numeric doubles
        rawCol = Temps.(colName);
        if iscell(rawCol)
            % If it's a cell array of strings/chars, convert to numeric
            if any(cellfun(@ischar, rawCol)) || any(cellfun(@isstring, rawCol))
                numericCol = str2double(rawCol);
            else
                numericCol = cell2mat(rawCol);
            end
        elseif isstring(rawCol) || ischar(rawCol)
            numericCol = str2double(rawCol);
        else
            numericCol = double(rawCol);
        end
        
        % Ensure it is a clean, flat column vector in the table
        Temps.(colName) = numericCol(:);
    end
end

% Verify all temp columns are now numeric double
disp('Data cleaning complete. Column types:');
summary(Temps)

function plotParity(Temps)
    % Extract and stack columns
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];

    % Create grouping vectors
    numRows = height(Temps);
    TC_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];

    % Calculate Performance Metrics
    mean_meas = mean(T_meas);
    SS_tot = sum((T_meas - mean_meas).^2);
    SS_res = sum((T_meas - T_calc).^2);
    R2 = 1 - (SS_res / SS_tot);
    RMSE = sqrt(mean((T_calc - T_meas).^2));
    MAE = mean(abs(T_calc - T_meas));

    % Plot setup
    figure('Name', 'Model Parity Analysis', 'Color', 'w', 'Position', [100, 100, 650, 550]);
    hold on;

    % --- 1. Draw the 20% Error Shaded Band (Background) ---
    lims = [min([T_meas; T_calc]) - 2, max([T_meas; T_calc]) + 2];
    x_band = linspace(lims(1), lims(2), 100);
    y_upper = 1.20 * x_band;
    y_lower = 0.80 * x_band;
   

    % --- 3. Draw the 1:1 Parity Line ---
    h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5);

    % --- 4. Plot the Actual Data Points (Foreground) ---
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(T_meas, T_calc, TC_labels, colors, 'o', 6, 'on');

    % Fill each group's markers with its own color
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end

    % --- 6. Formatting & Labels ---
    xlabel('Measured Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Calculated Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    % Update legend to include the new visual elements cleanly
    legendLabels = {'1:1 Agreement ($y=x$)', ...
                    'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend([h_parity, h_scatter'], legendLabels, ...
           'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 9);
    
    grid on; 
    axis equal; 
    xlim(lims); 
    ylim(lims);
end


function plotResiduals(Temps)
    % Stack and calculate residuals
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];
    residuals = T_calc - T_meas;

    figure('Name', 'Model Residuals', 'Color', 'w', 'Position', [150, 150, 650, 500]);
    histogram(residuals, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    hold on;

    % Fit and plot Normal Distribution
    pd = fitdist(residuals, 'Normal');
    x_val = linspace(min(residuals)-2, max(residuals)+2, 250);
    plot(x_val, pdf(pd, x_val), 'r-', 'LineWidth', 2);

    % LaTeX formatted stats block (uses double-backslash because of sprintf)
    statsStr = { ...
        sprintf('$\\mu_{\\mathrm{bias}} = %.3f\\,^{\\circ}\\mathrm{C}$', pd.mu), ...
        sprintf('$\\sigma = %.3f\\,^{\\circ}\\mathrm{C}$', pd.sigma) ...
    };
    annotation('textbox', [0.15, 0.70, 0.28, 0.13], 'String', statsStr, ...
               'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', [0.8 0.8 0.8], ...
               'Interpreter', 'latex', 'FontSize', 11);

    % Axis & Label formatting (uses single backslashes!)
    xlabel('Residual Error, $e = T_{\mathrm{calc}} - T_{\mathrm{meas}}$ ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Probability Density', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Distribution of Model Residuals}', 'Interpreter', 'latex', 'FontSize', 13);
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    grid on;
end

function plotSpatialError(Temps)
    % Calculate residuals per thermocouple location
    res_matrix = [Temps.T1_calc - Temps.T1_meas, ...
                  Temps.T2_calc - Temps.T2_meas, ...
                  Temps.T3_calc - Temps.T3_meas, ...
                  Temps.T4_calc - Temps.T4_meas, ...
                  Temps.T5_calc - Temps.T5_meas];

    figure('Name', 'Spatial Error Distribution', 'Color', 'w');
    
    % Create boxplot
    boxplot(res_matrix, 'Labels', {'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'}, ...
            'Widths', 0.5, 'Symbol', 'r+');
    hold on;

    % Benchmark line at 0 error
    yline(0, 'r--', 'LineWidth', 1.5);

    % Formatting
    ylabel('Temperature Error (T_{calc} - T_{meas}) [°C]');
    xlabel('Thermocouple Location');
    title('Model Residual Distribution across Spatial Coordinates', 'FontSize', 12);
    grid on;
end

function plotProfileComparison(Temps, selectedRuns)
    figure('Name', 'Profile Comparison', 'Color', 'w', 'Position', [200, 200, 780, 500]);
    hold on;
    x_tc = 1:5; 
    
    colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.929 0.694 0.125], ...
              [0.494 0.184 0.556], [0.466 0.674 0.188]};

    for i = 1:length(selectedRuns)
        runName = selectedRuns{i};
        
        rowIdx = find(strcmp(Temps.ParamName, runName));
        if isempty(rowIdx)
            warning('Parameter %s not found in Temps.', runName);
            continue;
        end
        
        % Extract profiles
        meas_vals = [Temps.T1_meas(rowIdx), Temps.T2_meas(rowIdx), ...
                     Temps.T3_meas(rowIdx), Temps.T4_meas(rowIdx), Temps.T5_meas(rowIdx)];
                 
        calc_vals = [Temps.T1_calc(rowIdx), Temps.T2_calc(rowIdx), ...
                     Temps.T3_calc(rowIdx), Temps.T4_calc(rowIdx), Temps.T5_calc(rowIdx)];
        
        c = colors{mod(i-1, length(colors)) + 1};
        
        % Escape underscores for LaTeX representation and set in monospaced font
        % (Use double backslash inside sprintf to output a single backslash to LaTeX)
        latexRunName = strrep(runName, '_', '\_');
        displayNameMeas = sprintf('\\texttt{%s} ($T_{\\mathrm{meas}}$)', latexRunName);
        displayNameCalc = sprintf('\\texttt{%s} ($T_{\\mathrm{calc}}$)', latexRunName);
        
        % Plot lines and symbols
        plot(x_tc, meas_vals, 'o', 'Color', c, 'MarkerFaceColor', c, ...
             'MarkerSize', 8, 'DisplayName', displayNameMeas);
        plot(x_tc, calc_vals, '-', 'Color', c, 'LineWidth', 2, ...
             'DisplayName', displayNameCalc);
    end

    % Formatting with LaTeX (uses single backslashes!)
    xlabel('Thermocouple Location', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Comparison of Spatial Thermal Profiles}', 'Interpreter', 'latex', 'FontSize', 13);
    
    % Set ticks and style them cleanly using LaTeX Roman font
    xticks(1:5);
    xticklabels({'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'});
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);

    legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 10);
    grid on;
end

%% --- Phi-fitting helper functions ---

function phi_fit = fitPhiForRow(i, T, x_targets, D_obs, per, h_fg, M, N, ...
                                 W, W_HT, H_wall, k_wall, x_obs, htc_amb, ...
                                 T_inf, T_amb, phi_bounds)
    % Finds the phi value for row i of T that minimizes the SSE between
    % computed and measured thermocouple temperatures at x_targets.

    q_eff_i   = T.q_eff(i);
    Wad_i     = T.W_adiabatic(i);
    mdot_f_i  = T.mdot_f(i);
    htc_dry_i = T.htc_dry(i);
    htc_wet_i = T.htc_wet(i);
    Tmeas_i   = [T.T1(i) T.T2(i) T.T3(i) T.T4(i) T.T5(i)];

    objFun = @(phi) phiResidual(phi, q_eff_i, Wad_i, mdot_f_i, htc_dry_i, htc_wet_i, ...
                    x_targets, D_obs, per, h_fg, M, N, W, W_HT, H_wall, k_wall, ...
                    x_obs, htc_amb, T_inf, T_amb, Tmeas_i);

    opts = optimset('TolX', 1e-4, 'Display', 'off');
    phi_fit = fminbnd(objFun, phi_bounds(1), phi_bounds(2), opts);
end

function sse = phiResidual(phi, q_eff_i, Wad_i, mdot_f_i, htc_dry_i, htc_wet_i, ...
                            x_targets, D_obs, per, h_fg, M, N, W, W_HT, H_wall, ...
                            k_wall, x_obs, htc_amb, T_inf, T_amb, Tmeas_i)
    % Forward model evaluated only at the 5 target downstream locations
    % (rather than all 22 nodes) for speed, since only these locations
    % contribute to the fit objective.

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    mdot_wake_i  = mdot_f_i * (D_obs/per) * phi;
    scalarPart_i = (q_eff_i * Wad_i) / (mdot_wake_i * h_fg);
    Bo_i         = scalarPart_i * x_targets;

    wake_ratio_i = max(0, -5 + 0.0058*Bo_i + ...
                   (6.3325*Bo_i.^1.6617)./(0.1719^1.6617 + Bo_i.^1.6617));
    dry_wake_W_i = Wad_i .* wake_ratio_i;
    w_free_i     = W - dry_wake_W_i;

    Tcalc_i = zeros(1, numel(x_targets));
    for k = 1:numel(x_targets)
        [~,~,~, Tc, ~,~,~,~,~,~] = ftoglass_ann_optim(M, N, W_HT, H_wall, k_wall, ...
            x_obs, q_eff_i, w_free_i(k), htc_dry_i, htc_wet_i, htc_amb, T_inf, T_amb);
        Tcalc_i(k) = Tc;
    end

    sse = sum((Tcalc_i - Tmeas_i).^2);
end

%%
% execute plots:
% Assumes your table is loaded in your workspace as "Temps"
plotParity(Temps);
plotResiduals(Temps);
plotSpatialError(Temps);

% Plot specific high/low-temperature scenarios to compare profiles
plotProfileComparison(Temps, {'PARAM01_04', 'PARAM04_03', 'PARAM05_04'});


%% --- Hydrodynamic Correlation for Phi (Heated Cases Only) ---

% 1. Filter out NaN, non-positive phi, AND zero-heating runs
% We use Q_heater > 5 to safely exclude adiabatic/zero-heat runs.
valid_idx = ~isnan(T.phi) & (T.phi > 0) & (T.Q_heater > 5);

phi_data = T.phi(valid_idx);
Re_data  = T.Re_f(valid_idx) .* (T.delta(valid_idx).^2);
SWR_data = T.SWR(valid_idx);

if isempty(phi_data)
    error('No valid heated runs with optimized phi values found.');
end

% 2. Direct Non-Linear Fit in Physical Space: phi = C * Re_f^a * SWR^b
% Initial guess for [C, a, b] -- adjust based on rough expectations
x0 = [1, 1, 1];
% Define the objective function to minimize the sum of squared errors in physical space
obj_fun = @(p) sum((phi_data - (p(1) * (Re_data.^p(2)) .* (SWR_data.^p(3)))).^2);

% Run fminsearch to get the optimal physical parameters
opts = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'Display', 'off');
fitted_params = fminsearch(obj_fun,x0, opts);

% Extract optimized parameters
C_opt = fitted_params(1);
a_opt = fitted_params(2);
b_opt = fitted_params(3);

% Calculate optimized predicted phi and the true physical R-squared
phi_pred_opt = C_opt * (Re_data.^a_opt) .* (SWR_data.^b_opt);
SS_tot = sum((phi_data - mean(phi_data)).^2);
SS_res = sum((phi_data - phi_pred_opt).^2);
R2_phi = 1 - (SS_res / SS_tot);

fprintf('\n==================================================\n');
fprintf('    PHI CORRELATION   \n');
fprintf('==================================================\n');
fprintf('Model Form: phi = C * Re_f^a * SWR^b\n\n');
fprintf('  C = %.6f\n', C_opt);
fprintf('  a = %.4f  [Re_f exponent]\n', a_opt);
fprintf('  b = %.4f  [SWR exponent]\n', b_opt);
fprintf('  R^2 Accuracy = %.4f (on %d active data points)\n', R2_phi, length(phi_data));
fprintf('==================================================\n\n');
%% --- Plotting the Dependency & Surface ---
figure('Name', 'Phi Hydrodynamic Correlation (Heated Runs Only)', 'Color', 'w', 'Position', [100, 100, 1100, 450]);

% --- Subplot 1: Phi vs Re_f (Colored by SWR) ---
subplot(1, 3, 1);
scatter(Re_data, phi_data, 50, SWR_data, 'filled', 'MarkerEdgeColor', 'k');
colormap(gca, 'parula');
cb1 = colorbar;
ylabel(cb1, 'SWR (Shear-to-Weight Ratio)', 'Interpreter', 'latex');
grid on;
xlabel('Film Reynolds Number, $Re_f$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Optimized $\phi$', 'Interpreter', 'latex', 'FontSize', 11);
title('Influence of $Re_f$ on $\phi$ (Heated Only)', 'Interpreter', 'latex', 'FontSize', 12);

% --- Subplot 2: Phi vs SWR (Colored by Re_f) ---
subplot(1, 3, 2);
scatter(SWR_data, phi_data, 50, Re_data, 'filled', 'MarkerEdgeColor', 'k');
colormap(gca, 'jet');
cb2 = colorbar;
ylabel(cb2, 'Re_f', 'Interpreter', 'latex');
grid on;
xlabel('Shear-to-Weight Ratio, $SWR$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Optimized $\phi$', 'Interpreter', 'latex', 'FontSize', 11);
title('Influence of $SWR$ on $\phi$ (Heated Only)', 'Interpreter', 'latex', 'FontSize', 12);

% --- Subplot 3: 3D Visualization of Correlation Surface ---
subplot(1, 3, 3);
hold on;

% Define grids for plotting the fitted surface
Re_grid = linspace(min(Re_data)*0.9, max(Re_data)*1.1, 40);
SWR_grid = linspace(min(SWR_data)*0.9, max(SWR_data)*1.1, 40);
[RE, SWR_m] = meshgrid(Re_grid, SWR_grid);
PHI_mesh = C_opt * (RE.^a_opt) .* (SWR_m.^b_opt);

% Plot fitted power-law surface
surf(RE, SWR_m, PHI_mesh, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
colormap(gca, 'viridis');

% Plot original optimized data points (excluding the 0.9999 artifacts!)
stem3(Re_data, SWR_data, phi_data, 'filled', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'Color', [0.3 0.3 0.3]);

grid on; view(3);
xlabel('$Re_f$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$SWR$', 'Interpreter', 'latex', 'FontSize', 11);
zlabel('$\phi$', 'Interpreter', 'latex', 'FontSize', 11);
title(sprintf('Correlation Fit ($R^2 = %.3f$)', R2_phi), 'Interpreter', 'latex', 'FontSize', 12);

% Format the surface equation onto the 3D plot
fitStr = sprintf('$\\phi = %.4f \\cdot Re_f^{%.3f} \\cdot SWR^{%.3f}$', C_opt, a_opt, b_opt);
annotation('textbox', [0.71, 0.02, 0.25, 0.08], 'String', fitStr, ...
           'Interpreter', 'latex', 'FontSize', 11, 'FitBoxToText', 'on', ...
           'BackgroundColor', 'w', 'EdgeColor', [0.8 0.8 0.8]);


%% --- Simplified Side-by-Side Diagnostic Plots ---
figure('Name', 'Phi Hydrodynamic Trends', 'Color', 'w', 'Position', [150, 150, 950, 420]);

% Style settings for clean, professional slide graphics
markerColor = [0.00, 0.45, 0.74]; % Solid professional blue
lineColor   = [0.85, 0.33, 0.10]; % Clean contrast red/orange

% --- Subplot 1: Phi vs. Re_f ---
subplot(1, 2, 1);
hold on;

% Plot raw data points as a single solid color
plot(Re_data, phi_data, 'o', 'MarkerFaceColor', markerColor, ...
     'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerSize', 7, 'DisplayName', 'Data Points');

% Simple single-variable power law fit: phi = C1 * Re_f^a1
fit_Re = fminsearch(@(p) sum((phi_data - p(1)*Re_data.^p(2)).^2), [0.1, 0.2]);
Re_smooth = linspace(min(Re_data), max(Re_data), 100);
phi_trend_Re = fit_Re(1) * Re_smooth.^fit_Re(2);

plot(Re_smooth, phi_trend_Re, '-', 'Color', lineColor, 'LineWidth', 2.5, ...
     'DisplayName', 'Power-Law Trend');

grid on;
xlabel('Film Reynolds Number, Re_f', 'FontSize', 12);
ylabel('Bypass Fraction, \phi', 'FontSize', 12);
title('\phi vs. Film Reynolds Number', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 10);
ylim([0 1.05]);
xlim([0 max(Re_data)*1.05]);

% --- Subplot 2: Phi vs. SWR ---
subplot(1, 2, 2);
hold on;

% Plot raw data points as a single solid color
plot(SWR_data, phi_data, 'o', 'MarkerFaceColor', markerColor, ...
     'MarkerEdgeColor', [0.2 0.2 0.2], 'MarkerSize', 7, 'DisplayName', 'Data Points');

% Simple single-variable power law fit: phi = C2 * SWR^b2
fit_SWR = fminsearch(@(p) sum((phi_data - p(1)*SWR_data.^p(2)).^2), [0.1, 0.4]);
SWR_smooth = linspace(min(SWR_data), max(SWR_data), 100);
phi_trend_SWR = fit_SWR(1) * SWR_smooth.^fit_SWR(2);

plot(SWR_smooth, phi_trend_SWR, '-', 'Color', lineColor, 'LineWidth', 2.5, ...
     'DisplayName', 'Power-Law Trend');

grid on;
xlabel('Shear-to-Weight Ratio, SWR', 'FontSize', 12);
ylabel('Bypass Fraction, \phi', 'FontSize', 12);
title('\phi vs. Shear-to-Weight Ratio', 'FontSize', 13);
legend('Location', 'southeast', 'FontSize', 10);
ylim([0 1.05]);
xlim([min(SWR_data)*0.9 max(SWR_data)*1.05]);
%% --- Plot: Phi vs. Film Thickness (delta) ---
% Extract valid film thickness data (matching your heated-only indices)
delta_data = T.delta(valid_idx) * 1e6; % Convert from meters to micrometers (um)

figure('Name', 'Phi vs Film Thickness', 'Color', 'w', 'Position', [200, 200, 550, 420]);
hold on;

% 1. Scatter plot of optimized phi vs. film thickness, colored by SWR
scatter(delta_data, phi_data, 60, SWR_data, 'filled', 'MarkerEdgeColor', 'k');
colormap(gca, 'parula');
cb = colorbar;
ylabel(cb, 'SWR (Shear-to-Weight)', 'Interpreter', 'latex');

% 2. Fit a quick trend line directly to delta to show the general trajectory
p_fit = polyfit(delta_data, phi_data, 1);
delta_smooth = linspace(min(delta_data), max(delta_data), 100);
phi_trend = polyval(p_fit, delta_smooth);

plot(delta_smooth, phi_trend, 'r--', 'LineWidth', 2, ...
     'DisplayName', 'Linear Trend');

grid on;
xlabel('Film Thickness, $\delta$ ($\mu\mathrm{m}$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Bypass Fraction, $\phi$', 'Interpreter', 'latex', 'FontSize', 12);
title('$\phi$ Sensitivity to Liquid Film Thickness', 'Interpreter', 'latex', 'FontSize', 13);
legend('Location', 'northwest', 'Interpreter', 'latex');