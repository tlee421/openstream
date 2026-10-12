% =========================================================================
% -- HTModelEvalFull.m --
% =========================================================================
% Purpose:
%   Evaluates a 2D obstruction wall heat transfer model (ftoglass_ann_optim)
%   against experimental temperature data and dry wake width measurements.
%
% Capabilities:
%   1. Parses measured dry wake widths from an experimental data structure.
%   2. Calculates empirical wake width predictions using closure relation.
%   3. Solves the 2D wall heat transfer solver across downstream locations.
%   4. Generates parity, residual, spatial, and cross-variable diagnostic plots.
%   5. Performs sensitivity analysis (±10% variations in h_wet, h_dry, W_dry).
%   6. Runs an inverse heat transfer optimization (fminbnd) to extract 
%      optimal h_dry values and correlates them against hydrodynamic parameters.
% =========================================================================

import Solvers.Obstruction.WallHT.ftoglass_ann_optim

% -------------------------------------------------------------------------
% 1. LOAD EXPERIMENTAL DATA & FLUID PROPERTIES
% -------------------------------------------------------------------------
% Map raw Excel spreadsheet headers to standardized script variable names
varNames = {'Param', 'T1', 'T2', 'T3', 'T4', 'T5', 'Tsat', 'G', 'x_in', ...
            'q_eff', 'Q_heater', 'delta', 'Re_f', 'SWR', 'mdot_f', 'tau_w', 'Re_v'};

% Read Excel file with externally computed information of flow state at obstruction
T = readtable("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\07_Conferences\ASME SHTC\Data\WallHTModelComparison.xlsx", ...
    'VariableNamingRule', 'preserve');   
T.Properties.VariableNames = varNames; % Overwrite table headers with standard names

% Load fluid thermophysical properties and experimental wake data structure
load("testfluidprops.mat") % Contains fluid property struct (KF, KG, RHOF, etc.)
load("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\02_MFVAL\Data Processing\Jason_FT_Processing_Program\film-thickness-master\Scripts\wake_width_datastructure.mat", "Master")

% Extract key fluid properties for saturated liquid and vapor states
k_film  = fluid.KF;         % [W/m-K] Thermal conductivity of saturated liquid
k_vapor = fluid.KG;         % [W/m-K] Thermal conductivity of saturated vapor
rho_f   = fluid.RHOF;       % [kg/m^3] Density of saturated liquid
mu_f    = fluid.MUF;        % [Pa-s] Dynamic viscosity of saturated liquid
nu_f    = mu_f/rho_f;       % [m^2/s] Kinematic viscosity of saturated liquid
PR_f    = fluid.PRANDTLF;   % [-] Liquid Prandtl number
PR_v    = fluid.PRANDTLG;   % [-] Vapor Prandtl number
h_fg    = fluid.HFG;        % [J/kg] Latent heat of vaporization

% Test Section & Obstruction Geometry
W     = 0.036;              % [m] Full channel width
H     = 0.012;              % [m] Full channel height
per   = 2*W + 2*H;          % [m] Wetted perimeter
A_c   = W*H;                % [m^2] Cross-sectional area
D_h   = 4*A_c/per;          % [m] Hydraulic diameter
L     = 0.525;              % [m] Total channel length
L_obs = 0.42;               % [m] Obstruction location from channel inlet
D_obs = 0.0025;             % [m] Obstruction diameter
x_obs = L - L_obs;          % [m] Axial distance from obstruction to end of test section (0.105 m)

% Compute Base Liquid Film Thickness using Law of the Wall (y+ = 5)
yplus        = 5;                                                          % Assumed viscous sublayer threshold
T.u_star     = (T.tau_w / rho_f).^(0.5);                                   % [m/s] Friction velocity based on wall shear stress
T.base_thick = yplus * (T.u_star).^(-1) * nu_f;                            % [m] Base film thickness

% Baseline Heat Transfer Coefficients (HTC)
T.htc_wet    = k_film * (T.base_thick).^(-1);                               % [W/m^2-K] Convective HTC in wet wall region
T.htc_dry    = k_vapor / D_h * 0.023 * (T.Re_v).^0.8 * PR_v^0.4;            % [W/m^2-K] Convective HTC in dry wake (Dittus-Boelter)

% -------------------------------------------------------------------------
% 2. PARSE MEASURED DRY WAKE WIDTHS FROM MASTER DATA STRUCTURE
% -------------------------------------------------------------------------
% Preallocate cell arrays to store 1x5 vectors corresponding to TC 1 through TC 5
T.dry_wake_W = cell(height(T), 1);
T.w_free     = cell(height(T), 1);

for i = 1:height(T)
    param_str   = string(T.Param(i));
    found_match = false;
    
    % Search nested data structure for image files matching current test parameter name
    for e = 1:length(Master.datasets(2).entries)
        wakes = Master.datasets(2).entries(e).wakes;
        for w = 1:length(wakes)
            if contains(string(wakes(w).imgFileName), param_str)
                % Sum left and right wake half-widths for TC1-TC5 (indices 2 through 6)
                measured_widths = wakes(w).widths(2:6, 1) + wakes(w).widths(2:6, 2);
                
                % Store measured dry wake width [m] and compute wetted free-stream width [m]
                T.dry_wake_W{i} = measured_widths(:)' * 1e-3; % Convert mm to meters
                T.w_free{i}     = W - T.dry_wake_W{i};       % W_free = Channel_Width - Dry_Wake_Width
                found_match     = true;
                break;
            end
        end
        if found_match; break; end
    end
    
    % Fallback handling if parameter name is missing from Master dataset
    if ~found_match
        warning('No wake data found in Master structure for parameter: %s', param_str);
        T.dry_wake_W{i} = zeros(1, 5);      % Default to zero dry wake
        T.w_free{i}     = W * ones(1, 5);  % Default to fully wetted wall
    end
end

% -------------------------------------------------------------------------
% 3. CALCULATE EMPIRICAL PREDICTED WAKE WIDTHS (CLOSURE CORRELATION)
% -------------------------------------------------------------------------
x_tc = [0.04, 0.05, 0.06, 0.07, 0.08]; % [m] Axial distances downstream of obstruction for TC1-TC5
T.dry_wake_W_calc = cell(height(T), 1);

% Compute adiabatic wake width using empirical correlation [m]
T.W_adiabatic     = (1.66 * T.Re_f.^0.085 .* T.SWR.^0.22) * 2 / 1000;   %Empirical adiabatic wake width closure relation
T.phi             = 0.53 * ones(height(T), 1);                          % Empirical wake partition factor

% Mass flow in wake region and Boiling number scaling
T.mdot_wake = T.mdot_f * (D_obs / per) .* T.phi;
scalarPart  = (T.q_eff .* T.W_adiabatic) ./ (T.mdot_wake .* h_fg);

% Compute predicted dry wake width profile for each test condition across TC locations
for i = 1:height(T)
    Bo_tc           = scalarPart(i) .* x_tc;                            % Boiling number at each TC x-location
    wake_ratio_func = -5 + 0.0058*Bo_tc + (6.3325 * Bo_tc.^1.6617) ./ (0.1719^1.6617 + Bo_tc.^1.6617); %dry wake width ratio closure relation
    wake_ratio      = max(0, wake_ratio_func);                          % Enforce non-negative wake ratio
    T.dry_wake_W_calc{i} = T.W_adiabatic(i) .* wake_ratio;             % [m] Calculated dry wake width
end

% -------------------------------------------------------------------------
% 4. SOLVE BASELINE HEAT TRANSFER MODEL AT DOWNSTREAM TC LOCATIONS
% -------------------------------------------------------------------------
% Grid discretization & wall parameters for 2D solver
M       = 125;              % Number of spanwise nodes
N       = 125;              % Number of wall-normal nodes
W_HT    = W/2;              % [m] Half-width of wall (symmetry boundary condition)
H_wall  = 0.003;            % [m] Glass wall thickness
k_wall  = 1.2;              % [W/m-K] Thermal conductivity of glass wall
htc_amb = 10;               % [W/m^2-K] External ambient natural convection HTC
T_inf   = fluid.TSAT;       % [K] Fluid saturation temperature
T_amb   = 22.5 + 273.15;    % [K] Ambient environment temperature
N_rows  = height(T);

T.T_center = cell(N_rows, 1);

% Execute 2D conduction solver for each test condition using measured wake widths
for i = 1:N_rows
    w_free_i     = T.w_free{i};          
    N_wfree      = numel(w_free_i);
    T_center_row = zeros(1, N_wfree);
    
    % Solve for outer wall centerline temperature at each thermocouple x-location
    for j = 1:N_wfree
        [~, ~, ~, T_center_row(j), ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
            M, N, W_HT, H_wall, k_wall, x_obs, ...
            T.q_eff(i), w_free_i(j), T.htc_dry(i), T.htc_wet(i), ...
            htc_amb, T_inf, T_amb);
    end
    T.T_center{i} = T_center_row;
end

% Construct structured evaluation table for comparative plotting
Temps = table();
Temps.ParamName = string(T.Param);
Temps.T1_meas   = double(T.T1);
Temps.T2_meas   = double(T.T2);
Temps.T3_meas   = double(T.T3);
Temps.T4_meas   = double(T.T4);
Temps.T5_meas   = double(T.T5);

% Extract calculated temperatures into matrix form and assign to table
T_calc_mat      = vertcat(T.T_center{:});
Temps.T1_calc   = T_calc_mat(:, 1);
Temps.T2_calc   = T_calc_mat(:, 2);
Temps.T3_calc   = T_calc_mat(:, 3);
Temps.T4_calc   = T_calc_mat(:, 4);
Temps.T5_calc   = T_calc_mat(:, 5);

disp('Data pre-processing and baseline solution complete.');

% -------------------------------------------------------------------------
% 5. EXECUTE ALL PLOTS, SENSITIVITY, AND OPTIMIZATION
% -------------------------------------------------------------------------
% Baseline Temperature Parity & Error Distribution Plots
plotParity(Temps);
%plotResiduals(Temps);
plotSpatialError(Temps);

% Spatial Thermal Profile Comparison Plots
%plotProfileComparison(Temps, {'PARAM08_02', 'PARAM08_03', 'PARAM08_04'});
plotAllProfilesSummary(Temps);

% Wake Width Parity & Spatial Profile Plots
plotWakeWidthParity(T);
%plotWakeWidthProfileComparison(T, {'PARAM01_04', 'PARAM04_03', 'PARAM05_03'});
plotAllWakeWidthProfilesSummary(T);

% Cross-Variable Interaction Plots (Temperature vs. Wake Width)
%plotWakeWidthVsTempMeas(T, Temps);
%plotWakeWidthVsTempCalc(T, Temps);
%plotWakeWidthCalcVsTempCalc(T, Temps);

% Perturbation Sensitivity Analysis (+/- 10% on h_wet, h_dry, W_meas)
runSensitivityAnalysis(T, Temps, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb);

% Inverse Optimization (Extract optimal h_dry and plot parametric correlations)
T = optimizeDryHTC(T, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb);
plotHTCParametricCorrelations(T, fluid, D_h);


%% =========================================================================
% PLOTTING & ANALYSIS SUBFUNCTIONS
% =========================================================================

% -------------------------------------------------------------------------
% PLOT: Temperature Parity (Measured vs Calculated)
% -------------------------------------------------------------------------
function plotParity(Temps)
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];
    numRows = height(Temps);
    
    % Grouping labels for TC 1 through TC 5
    TC_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Model Parity Analysis', 'Color', 'w', 'Position', [100, 100, 650, 550]);
    hold on;
    
    lims = [min([T_meas; T_calc]) - 2, max([T_meas; T_calc]) + 2];
    h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5); % 1:1 Agreement line
    
    % Distinct color palette for TC 1 to TC 5
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(T_meas, T_calc, TC_labels, colors, 'o', 6, 'on');
    
    % Fill scatter marker centers
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end
    
    xlabel('Measured Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Calculated Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legendLabels = {'1:1 Agreement ($y=x$)', 'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend([h_parity, h_scatter'], legendLabels, 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 9);
    grid on; axis equal; xlim(lims); ylim(lims);
end

% -------------------------------------------------------------------------
% PLOT: Model Residual Distribution Histogram & Gaussian Fit
% -------------------------------------------------------------------------
function plotResiduals(Temps)
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];
    residuals = T_calc - T_meas;
    
    figure('Name', 'Model Residuals', 'Color', 'w', 'Position', [150, 150, 650, 500]);
    histogram(residuals, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    hold on;
    
    % Fit normal probability density function to residual vector
    pd = fitdist(residuals, 'Normal');
    x_val = linspace(min(residuals)-2, max(residuals)+2, 250);
    plot(x_val, pdf(pd, x_val), 'r-', 'LineWidth', 2);
    
    % Annotation box displaying mean bias error (mu) and standard deviation (sigma)
    statsStr = { sprintf('$\\mu_{\\mathrm{bias}} = %.3f\\,^{\\circ}\\mathrm{C}$', pd.mu), ...
                 sprintf('$\\sigma = %.3f\\,^{\\circ}\\mathrm{C}$', pd.sigma) };
    annotation('textbox', [0.15, 0.70, 0.28, 0.13], 'String', statsStr, ...
               'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', [0.8 0.8 0.8], ...
               'Interpreter', 'latex', 'FontSize', 11);
               
    xlabel('Residual Error, $e = T_{\mathrm{calc}} - T_{\mathrm{meas}}$ ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Probability Density', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Distribution of Model Residuals}', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Spatial Residual Error Distribution (Boxplots per TC)
% -------------------------------------------------------------------------
function plotSpatialError(Temps)
    res_matrix = [Temps.T1_calc - Temps.T1_meas, Temps.T2_calc - Temps.T2_meas, ...
                  Temps.T3_calc - Temps.T3_meas, Temps.T4_calc - Temps.T4_meas, ...
                  Temps.T5_calc - Temps.T5_meas];
              
    figure('Name', 'Spatial Error Distribution', 'Color', 'w');
    boxplot(res_matrix, 'Labels', {'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'}, 'Widths', 0.5, 'Symbol', 'r+');
    hold on;
    
    yline(0, 'r--', 'Zero Error Baseline', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right', 'FontSize', 10);
    ylabel('Temperature Error (T_{calc} - T_{meas}) [°C]', 'FontWeight', 'bold');
    xlabel('Thermocouple Location', 'FontWeight', 'bold');
    title('Model Residual Distribution across Spatial Coordinates', 'FontSize', 12);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Spatial Profile Comparisons for Specific Test Runs
% -------------------------------------------------------------------------
function plotProfileComparison(Temps, selectedRuns)
    figure('Name', 'Profile Comparison', 'Color', 'w', 'Position', [200, 200, 780, 500]);
    hold on;
    x_tc = 1:5; 
    colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.929 0.694 0.125], [0.494 0.184 0.556], [0.466 0.674 0.188]};
    
    for i = 1:length(selectedRuns)
        runName = selectedRuns{i};
        rowIdx = find(strcmp(Temps.ParamName, runName));
        if isempty(rowIdx)
            warning('Parameter %s not found in Temps.', runName);
            continue;
        end
        
        meas_vals = [Temps.T1_meas(rowIdx), Temps.T2_meas(rowIdx), Temps.T3_meas(rowIdx), Temps.T4_meas(rowIdx), Temps.T5_meas(rowIdx)];
        calc_vals = [Temps.T1_calc(rowIdx), Temps.T2_calc(rowIdx), Temps.T3_calc(rowIdx), Temps.T4_calc(rowIdx), Temps.T5_calc(rowIdx)];
        c = colors{mod(i-1, length(colors)) + 1};
        
        latexRunName = strrep(runName, '_', '\_');
        displayNameMeas = sprintf('\\texttt{%s} ($T_{\\mathrm{meas}}$)', latexRunName);
        displayNameCalc = sprintf('\\texttt{%s} ($T_{\\mathrm{calc}}$)', latexRunName);
        
        % Plot discrete symbols for measurement, solid lines for model
        plot(x_tc, meas_vals, 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 8, 'DisplayName', displayNameMeas);
        plot(x_tc, calc_vals, '-', 'Color', c, 'LineWidth', 2, 'DisplayName', displayNameCalc);
    end
    
    xlabel('Thermocouple Location', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Comparison of Spatial Thermal Profiles}', 'Interpreter', 'latex', 'FontSize', 13);
    xticks(1:5);
    xticklabels({'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'});
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 10);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Aggregate Ensemble Profile Comparison (All Runs + Ensemble Means)
% -------------------------------------------------------------------------
function plotAllProfilesSummary(Temps)
    T_meas_all = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas];
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc];
    
    figure('Name', 'Aggregate Profile Comparison', 'Color', 'w', 'Position', [250, 250, 800, 500]);
    hold on;
    
    % Plot low-opacity individual run traces
    for i = 1:height(Temps)
        plot(1:5, T_meas_all(i,:), '-', 'Color', [0 0.447 0.741 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(1:5, T_calc_all(i,:), '-', 'Color', [0.850 0.325 0.098 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    % Overlay thick ensemble mean trends
    h_meas = plot(1:5, mean(T_meas_all, 1, 'omitnan'), '-o', 'Color', [0 0.447 0.741], ...
        'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', [0 0.447 0.741]);
    h_calc = plot(1:5, mean(T_calc_all, 1, 'omitnan'), '-^', 'Color', [0.850 0.325 0.098], ...
        'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', [0.850 0.325 0.098]);
    
    xlabel('Thermocouple Location', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Aggregate Spatial Thermal Profiles (All Data)}', 'Interpreter', 'latex', 'FontSize', 13);
    
    xticks(1:5);
    xticklabels({'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'});
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legend([h_meas, h_calc], {'Mean Measured', 'Mean Calculated'}, ...
        'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 11);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Parity Analysis for Dry Wake Widths (Measured vs Empirical Model)
% -------------------------------------------------------------------------
function plotWakeWidthParity(T)
    if ~ismember('dry_wake_W_calc', T.Properties.VariableNames) || ~ismember('dry_wake_W', T.Properties.VariableNames)
        error('No wake width data found. Ensure T.dry_wake_W and T.dry_wake_W_calc are populated in table T.');
    end
    meas_all  = [];
    calc_all  = [];
    tc_labels = [];
    
    % Flatten cell array entries into single vectors
    for i = 1:height(T)
        meas_i = T.dry_wake_W{i};
        calc_i = T.dry_wake_W_calc{i};
        
        if ~isempty(meas_i) && ~isempty(calc_i)
            meas_all  = [meas_all, meas_i];
            calc_all  = [calc_all, calc_i];
            tc_labels = [tc_labels, 1:5];
        end
    end
    
    if isempty(meas_all)
        error('No non-empty wake width entries were found to plot.');
    end
    
    figure('Name', 'Wake Width Parity Analysis', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    lims = [min([meas_all, calc_all])*0.9, max([meas_all, calc_all])*1.1];
    if min(lims) == 0; lims(1) = -0.001; end 
    
    h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5);
    
    colors    = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(meas_all, calc_all, tc_labels, colors, 'o', 6, 'on');
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end
    
    xlabel('Measured Dry Wake Width, $W_{\mathrm{meas}}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Empirical Predicted Width, $W_{\mathrm{calc}}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Parity Plot: Dry Wake Width}', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legendLabels = {'1:1 Agreement ($y=x$)', 'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend([h_parity, h_scatter'], legendLabels, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    
    grid on; axis equal; xlim(lims); ylim(lims);
end

% -------------------------------------------------------------------------
% PLOT: Dry Wake Width Spatial Profiles for Selected Runs
% -------------------------------------------------------------------------
function plotWakeWidthProfileComparison(T, selectedRuns)
    figure('Name', 'Wake Width Profile Comparison', 'Color', 'w', 'Position', [200, 200, 780, 500]);
    hold on;
    x_tc   = 1:5; 
    colors = {[0 0.447 0.741], [0.85 0.325 0.098], [0.929 0.694 0.125], ...
              [0.494 0.184 0.556], [0.466 0.674 0.188]};
          
    for i = 1:length(selectedRuns)
        runName = selectedRuns{i};
        rowIdx = find(strcmp(string(T.Param), runName));
        if isempty(rowIdx)
            warning('Parameter %s not found in T.', runName);
            continue;
        end
        
        meas_vals = T.dry_wake_W{rowIdx};
        calc_vals = T.dry_wake_W_calc{rowIdx};
        c = colors{mod(i-1, length(colors)) + 1};
        
        latexRunName    = strrep(runName, '_', '\_');
        displayNameMeas = sprintf('\\texttt{%s} ($W_{\\mathrm{meas}}$)', latexRunName);
        displayNameCalc = sprintf('\\texttt{%s} ($W_{\\mathrm{calc}}$)', latexRunName);
        
        plot(x_tc, meas_vals, 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 8, 'DisplayName', displayNameMeas);
        plot(x_tc, calc_vals, '-', 'Color', c, 'LineWidth', 2, 'DisplayName', displayNameCalc);
    end
    
    xlabel('Thermocouple Location', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Dry Wake Width, $W$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Comparison of Spatial Dry Wake Width Profiles}', 'Interpreter', 'latex', 'FontSize', 13);
    xticks(1:5);
    xticklabels({'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'});
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 10);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Ensemble Aggregate Dry Wake Width Profiles
% -------------------------------------------------------------------------
function plotAllWakeWidthProfilesSummary(T)
    W_meas_all = vertcat(T.dry_wake_W{:});
    W_calc_all = vertcat(T.dry_wake_W_calc{:});
    
    figure('Name', 'Aggregate Wake Width Profile Comparison', 'Color', 'w', 'Position', [250, 250, 800, 500]);
    hold on;
    
    % Individual dataset traces
    for i = 1:height(T)
        if ~isempty(T.dry_wake_W{i}) && ~isempty(T.dry_wake_W_calc{i})
            plot(1:5, T.dry_wake_W{i}, '-', 'Color', [0 0.447 0.741 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(1:5, T.dry_wake_W_calc{i}, '-', 'Color', [0.850 0.325 0.098 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    
    % Ensemble averages across all runs
    h_meas = plot(1:5, mean(W_meas_all, 1, 'omitnan'), '-o', 'Color', [0 0.447 0.741], ...
        'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', [0 0.447 0.741]);
    h_calc = plot(1:5, mean(W_calc_all, 1, 'omitnan'), '-^', 'Color', [0.850 0.325 0.098], ...
        'LineWidth', 3, 'MarkerSize', 8, 'MarkerFaceColor', [0.850 0.325 0.098]);
    
    xlabel('Thermocouple Location', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Dry Wake Width, $W$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Aggregate Spatial Dry Wake Width Profiles (All Data)}', 'Interpreter', 'latex', 'FontSize', 13);
    
    xticks(1:5);
    xticklabels({'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'});
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legend([h_meas, h_calc], {'Mean Measured $W_{\mathrm{meas}}$', 'Mean Calculated $W_{\mathrm{calc}}$'}, ...
        'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 11);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Measured Wake Width vs. Measured Temperature
% -------------------------------------------------------------------------
function plotWakeWidthVsTempMeas(T, Temps)
    W_meas_all = vertcat(T.dry_wake_W{:}); 
    T_meas_all = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas]; 
    
    W_flat = W_meas_all(:);
    T_flat = T_meas_all(:);
    
    numRows   = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Wake Width vs Measured Temp', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    colors    = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(W_flat, T_flat, tc_labels, colors, 'o', 6, 'on');
    
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end
    
    xlabel('Measured Dry Wake Width, $W_{\mathrm{meas}}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Measured Wall Temperature, $T_{\mathrm{meas}}$ ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Measured Dry Wake Width vs. Wall Temperature}', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legendLabels = {'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend(h_scatter, legendLabels, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Measured Wake Width vs. Calculated Temperature
% -------------------------------------------------------------------------
function plotWakeWidthVsTempCalc(T, Temps)
    W_meas_all = vertcat(T.dry_wake_W{:}); 
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc]; 
    
    W_flat = W_meas_all(:);
    T_flat = T_calc_all(:);
    
    numRows   = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Wake Width vs Calculated Temp', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    colors    = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(W_flat, T_flat, tc_labels, colors, 'o', 6, 'on');
    
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end
    
    xlabel('Measured Dry Wake Width, $W_{\mathrm{meas}}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Calculated Wall Temperature, $T_{\mathrm{calc}}$ ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Measured Dry Wake Width vs. Calculated Wall Temperature}', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legendLabels = {'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend(h_scatter, legendLabels, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
end

% -------------------------------------------------------------------------
% PLOT: Calculated Wake Width vs. Calculated Temperature
% -------------------------------------------------------------------------
function plotWakeWidthCalcVsTempCalc(T, Temps)
    W_calc_all = vertcat(T.dry_wake_W_calc{:}); 
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc]; 
    
    W_flat = W_calc_all(:);
    T_flat = T_calc_all(:);
    
    numRows   = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Calc Wake Width vs Calc Temp', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    colors    = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(W_flat, T_flat, tc_labels, colors, 'o', 6, 'on');
    
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end
    
    xlabel('Calculated Dry Wake Width, $W_{\mathrm{calc}}$ (m)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Calculated Wall Temperature, $T_{\mathrm{calc}}$ ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Calculated Dry Wake Width vs. Calculated Wall Temperature}', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legendLabels = {'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend(h_scatter, legendLabels, 'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
end

% -------------------------------------------------------------------------
% SENSITIVITY ANALYSIS: Perturb Input Variables by +/-10% & Evaluate Metrics
% -------------------------------------------------------------------------
function runSensitivityAnalysis(T, Temps, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb)

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    N_rows     = height(T);
    T_meas_mat = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas];
    
    % Perturbation test matrix: {Scenario Name, h_wet factor, h_dry factor, W_meas factor}
    scenarios = {
        'Baseline',       1.00, 1.00, 1.00;
        'h_wet +10%',     1.10, 1.00, 1.00;
        'h_wet -10%',     0.90, 1.00, 1.00;
        'h_dry +10%',     1.00, 1.10, 1.00;
        'h_dry -10%',     1.00, 0.90, 1.00;
        'W_wake +10%',    1.00, 1.00, 1.10;
        'W_wake -10%',    1.00, 1.00, 0.90;
    };
    
    numScenarios = size(scenarios, 1);
    results      = struct();
    
    fprintf('\n=========================================================\n');
    fprintf('  HEAT TRANSFER MODEL SENSITIVITY ANALYSIS (PM 10%%)\n');
    fprintf('=========================================================\n');
    
    % Sweep across all scenarios
    for s = 1:numScenarios
        sName  = scenarios{s, 1};
        f_wet  = scenarios{s, 2};
        f_dry  = scenarios{s, 3};
        f_wake = scenarios{s, 4};
        
        T_calc_mat = zeros(N_rows, 5);
        
        for i = 1:N_rows
            W_wake_i  = T.dry_wake_W{i} * f_wake;
            w_free_i  = W - W_wake_i;
            htc_wet_i = T.htc_wet(i) * f_wet;
            htc_dry_i = T.htc_dry(i) * f_dry;
            
            % Solve thermal conduction model for perturbed inputs
            for j = 1:5
                [~, ~, ~, T_center, ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
                    M, N, W_HT, H_wall, k_wall, x_obs, ...
                    T.q_eff(i), w_free_i(j), htc_dry_i, htc_wet_i, ...
                    htc_amb, T_inf, T_amb);
                T_calc_mat(i, j) = T_center;
            end
        end
        
        % Compute summary statistical error metrics
        err      = T_calc_mat - T_meas_mat;
        rmse_val = sqrt(mean(err(:).^2, 'omitnan'));                     % Root Mean Square Error [°C]
        mae_val  = mean(abs(err(:)), 'omitnan');                         % Mean Absolute Error [°C]
        mbe_val  = mean(err(:), 'omitnan');                              % Mean Bias Error [°C]
        
        results(s).name   = sName;
        results(s).T_calc = T_calc_mat;
        results(s).RMSE   = rmse_val;
        results(s).MAE    = mae_val;
        results(s).MBE    = mbe_val;
        
        fprintf('%-14s | RMSE: %6.2f °C | MAE: %6.2f °C | MBE: %6.2f °C\n', ...
            sName, rmse_val, mae_val, mbe_val);
    end
    fprintf('=========================================================\n\n');
    
    % --- PLOT 1: INDIVIDUAL PARITY SENSITIVITY PLOTS (1x3 TILED LAYOUT) ---
    figure('Name', 'Sensitivity Parity Bounds - Individual Parameters', ...
           'Color', 'w', 'Position', [100, 150, 1350, 430]);
    t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Subplot configurations for each variable
    param_configs = {
        '$h_{\mathrm{wet}} \pm 10\%$ Sensitivity',     3, 2, [0.85 0.33 0.10];
        '$h_{\mathrm{dry}} \pm 10\%$ Sensitivity',     5, 4, [0.93 0.69 0.13];
        '$W_{\mathrm{meas}} \pm 10\%$ Sensitivity',    7, 6, [0.49 0.18 0.56];
    };
    
    T_meas_flat = T_meas_mat(:);
    lims        = [min(T_meas_flat) - 2, max(T_meas_flat) + 2];
    
    for p = 1:3
        nexttile;
        hold on;
        
        pTitle   = param_configs{p, 1};
        idx_low  = param_configs{p, 2};
        idx_high = param_configs{p, 3};
        pColor   = param_configs{p, 4};
        
        h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Agreement');
        
        % Draw vertical sensitivity range lines per data point
        plot([T_meas_flat, T_meas_flat]', [results(idx_low).T_calc(:), results(idx_high).T_calc(:)]', ...
            '-', 'Color', [pColor, 0.40], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        h_base  = scatter(T_meas_flat, results(1).T_calc(:), 25, [0 0.45 0.74], 'filled', 'DisplayName', 'Baseline');
        h_range = plot(nan, nan, '-', 'Color', pColor, 'LineWidth', 2, 'DisplayName', '$\pm 10\%$ Range');
        
        xlabel('Measured Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 11);
        ylabel('Calculated Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 11);
        title(['\textbf{', pTitle, '}'], 'Interpreter', 'latex', 'FontSize', 12);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
        legend([h_parity, h_base, h_range], 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 8);
        grid on; axis equal; xlim(lims); ylim(lims);
    end
    
    title(t, '\textbf{Individual Parameter Sensitivity Analysis on Temperature Parity}', ...
        'Interpreter', 'latex', 'FontSize', 14);
    
    % --- PLOT 2: RMSE BAR CHART ---
    figure('Name', 'Sensitivity RMSE Impact', 'Color', 'w', 'Position', [200, 200, 750, 450]);
    names = {results.name};
    rmses = [results.RMSE];
    
    b = bar(rmses, 'FaceColor', 'flat');
    b.CData(1, :)   = [0 0.45 0.74];          % Baseline bar color
    b.CData(2:3, :) = repmat([0.85 0.33 0.10], 2, 1); % h_wet bars
    b.CData(4:5, :) = repmat([0.93 0.69 0.13], 2, 1); % h_dry bars
    b.CData(6:7, :) = repmat([0.49 0.18 0.56], 2, 1); % W_wake bars
    
    % Overlay numerical RMSE values above each bar
    for k = 1:numScenarios
        text(k, rmses(k) + max(rmses)*0.02, sprintf('%.2f°C', rmses(k)), ...
            'HorizontalAlignment', 'center', 'Interpreter', 'latex', 'FontSize', 9);
    end
    
    set(gca, 'XTick', 1:numScenarios, 'XTickLabel', names, 'TickLabelInterpreter', 'none', 'FontSize', 10);
    ylabel('Root Mean Square Error, $\mathrm{RMSE}$ ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title('\textbf{Model Sensitivity Ranking (RMSE Impact)}', 'Interpreter', 'latex', 'FontSize', 13);
    ylim([0, max(rmses)*1.15]);
    grid on;
end

% -------------------------------------------------------------------------
% OPTIMIZATION: Row-by-Row Inverse Solver for Optimal Dry HTC (h_dry)
% -------------------------------------------------------------------------
function T = optimizeDryHTC(T, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb)
    
    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    N_rows          = height(T);
    T.htc_dry_opt   = zeros(N_rows, 1);
    T.RMSE_htc_opt  = zeros(N_rows, 1);
    
    h_min   = 5;       % [W/m^2-K] Lower bound for optimization search
    h_max   = 5000;    % [W/m^2-K] Upper bound for optimization search
    options = optimset('Display', 'off', 'TolX', 1e-1); % Bounded scalar optimization settings
    
    fprintf('\n=========================================================\n');
    fprintf('  RUNNING INVERSE HTC OPTIMIZATION (ROW-BY-ROW)\n');
    fprintf('=========================================================\n');
    
    for i = 1:N_rows
        T_meas_i = [T.T1(i), T.T2(i), T.T3(i), T.T4(i), T.T5(i)];
        w_free_i = T.w_free{i};
        
        % Objective handle: minimize temperature SSE across TC 1-5 for condition i
        objFun = @(h_dry) calcRowSSE(h_dry, T_meas_i, w_free_i, T.q_eff(i), T.htc_wet(i), ...
            M, N, W_HT, H_wall, k_wall, x_obs, htc_amb, T_inf, T_amb);
        
        % Run Golden Section / parabolic interpolation bounded solver
        [h_opt, sse_val] = fminbnd(objFun, h_min, h_max, options);
        
        T.htc_dry_opt(i)  = h_opt;
        T.RMSE_htc_opt(i) = sqrt(sse_val / 5); % Remaining RMSE per run after optimization
        
        if mod(i, 10) == 0 || i == N_rows
            fprintf('  Processed %d/%d test conditions...\n', i, N_rows);
        end
    end
    fprintf('Inverse HTC Optimization Complete!\n');
    fprintf('=========================================================\n\n');
end

% -------------------------------------------------------------------------
% OBJECTIVE FUNCTION: Sum of Squared Errors (SSE) across TC 1-5
% -------------------------------------------------------------------------
function sse = calcRowSSE(h_dry, T_meas, w_free, q_eff, htc_wet, ...
    M, N, W_HT, H_wall, k_wall, x_obs, htc_amb, T_inf, T_amb)

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    T_calc = zeros(1, 5);
    for j = 1:5
        [~, ~, ~, T_center, ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
            M, N, W_HT, H_wall, k_wall, x_obs, ...
            q_eff, w_free(j), h_dry, htc_wet, ...
            htc_amb, T_inf, T_amb);
        T_calc(j) = T_center;
    end
    sse = sum((T_calc - T_meas).^2); % Sum of squared temperature differences [°C^2]
end

% -------------------------------------------------------------------------
% PLOT: Multi-Panel Parametric Correlations for Backed-Out Dry HTC
% -------------------------------------------------------------------------
function plotHTCParametricCorrelations(T, fluid, D_h)
    k_vapor = fluid.KG;
    PR_v    = fluid.PRANDTLG;
    
    % Calculate Nusselt Numbers for optimized and baseline correlations
    Nu_opt = (T.htc_dry_opt .* D_h) ./ k_vapor;
    Nu_DB  = 0.023 .* (T.Re_v).^0.8 .* (PR_v)^0.4;
    
    figure('Name', 'Optimized HTC Parametric Correlations', 'Color', 'w', 'Position', [100, 100, 1200, 750]);
    t = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % 1. h_dry vs Re_v (Log-Log Forced Convection Power Fit)
    nexttile;
    loglog(T.Re_v, T.htc_dry_opt, 'o', 'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74]);
    hold on;
    loglog(T.Re_v, T.htc_dry, 's', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
    
    % Power law fit: log(h) = m*log(Re_v) + b
    valid_idx = T.htc_dry_opt > 10;
    p = polyfit(log(T.Re_v(valid_idx)), log(T.htc_dry_opt(valid_idx)), 1);
    Re_fit = linspace(min(T.Re_v), max(T.Re_v), 100);
    loglog(Re_fit, exp(p(2)) * Re_fit.^p(1), 'k--', 'LineWidth', 1.5);
    
    xlabel('Vapor Reynolds Number, $\mathrm{Re}_v$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry, opt}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title(sprintf('\\textbf{1. vs. $\\mathrm{Re}_v$ (Power Exponent: %.2f)}', p(1)), 'Interpreter', 'latex', 'FontSize', 12);
    legend({'Optimized', 'Dittus-Boelter', sprintf('Fit: $\\propto \\mathrm{Re}_v^{%.2f}$', p(1))}, ...
        'Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    % 2. h_dry vs Shear to Weight Ratio (SWR)
    nexttile;
    plot(T.SWR, T.htc_dry_opt, 'd', 'Color', [0.49 0.18 0.56], 'MarkerFaceColor', [0.49 0.18 0.56]);
    xlabel('Sheat to Weight Ratio, $\mathrm{SWR}$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry, opt}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{2. vs. Shear to Weight Ratio (SWR)}', 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    
    % 3. h_dry vs Effective Heat Flux (q_eff)
    nexttile;
    plot(T.q_eff / 1000, T.htc_dry_opt, '^', 'Color', [0.85 0.33 0.10], 'MarkerFaceColor', [0.85 0.33 0.10]);
    xlabel('Effective Heat Flux, $q''_{\mathrm{eff}}$ ($\mathrm{kW/m}^2$)', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry, opt}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{3. vs. Heat Flux $q''_{\mathrm{eff}}$}', 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    
    % 4. h_dry vs Inlet Vapor Quality (x_in)
    nexttile;
    plot(T.x_in, T.htc_dry_opt, 's', 'Color', [0.47 0.67 0.19], 'MarkerFaceColor', [0.47 0.67 0.19]);
    xlabel('Inlet Quality, $x_{\mathrm{in}}$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry, opt}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{4. vs. Inlet Quality $x_{\mathrm{in}}$}', 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    
    % 5. h_dry vs Liquid Film Reynolds Number (Re_f)
    nexttile;
    plot(T.Re_f, T.htc_dry_opt, 'v', 'Color', [0.93 0.69 0.13], 'MarkerFaceColor', [0.93 0.69 0.13]);
    xlabel('Liquid Film Reynolds Number, $\mathrm{Re}_f$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry, opt}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{5. vs. Liquid Film $\mathrm{Re}_f$}', 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    
    % 6. Nusselt Enhancement Ratio (Nu_opt / Nu_DB) vs Re_v
    nexttile;
    plot(T.Re_v, Nu_opt ./ Nu_DB, 'p', 'Color', [0.30 0.75 0.93], 'MarkerFaceColor', [0.30 0.75 0.93], 'MarkerSize', 7);
    yline(1.0, 'r--', 'Dittus-Boelter Baseline', 'LineWidth', 1.5);
    xlabel('Vapor Reynolds Number, $\mathrm{Re}_v$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('Enhancement Factor, $\mathrm{Nu}_{\mathrm{opt}} / \mathrm{Nu}_{\mathrm{DB}}$', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{6. Nusselt Enhancement Ratio}', 'Interpreter', 'latex', 'FontSize', 12);
    grid on;
    
    title(t, '\textbf{Parametric Scaling Analysis of Backed-Out Dry Heat Transfer Coefficient}', ...
        'Interpreter', 'latex', 'FontSize', 14);
end