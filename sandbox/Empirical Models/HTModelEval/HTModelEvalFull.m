% =========================================================================
% -- HTModelEvalFull.m --
% =========================================================================
% Purpose:
%   Evaluates a 2D obstruction wall heat transfer model (ftoglass_ann_optim)
%   against experimental temperature data and dry wake width measurements.
%
% Capabilities:
%   1. Parses measured dry wake widths from an experimental data structure and filters outliers.
%   2. Calculates empirical wake width predictions using closure relation.
%   3. Solves the 2D wall heat transfer solver across downstream locations.
%   4. Generates parity, residual, spatial, and cross-variable diagnostic plots.
%   5. Performs sensitivity analysis (±10% variations in h_wet, h_dry, W_dry).
%   6. Runs an inverse heat transfer optimization (fminbnd) to extract 
%      optimal h_dry values and correlates them against hydrodynamic parameters.
% =========================================================================

import Solvers.Obstruction.WallHT.ftoglass_ann_optim

% -------------------------------------------------------------------------
% CONTROL FLAGS & FILTER CONFIGURATION (EDIT FILTERING HERE)
% -------------------------------------------------------------------------
filter_outliers = true; % Set to 'true' to enable wake width & TC outlier filtering

% 1. Parameters to exclude ENTIRELY (All 5 thermocouples discarded and set to NaN)
full_exclude_params = { ...
    'PARAM_07_03', ...
    'PARAM_07_04', ...
    'PARAM_08_02', ...
    'PARAM_10_02', ...
    'PARAM_15_02', ...
    'PARAM_13_02', ...
    'PARAM_11_04', ...
    'PARAM_12_03', ... 
    'PARAM_06_04', ...
};

% 2. Specific Thermocouple Exclusions per Parameter
% Format: {'PARAM_NAME', [TC_indices_to_exclude]} (TC indices range from 1 to 5)
specific_tc_exclude = { ...
    'PARAM_05_04', 5;          % Exclude TC 5
    'PARAM_03_03', 2;          % Exclude TC 2
    'PARAM_09_02', [1, 2];     % Exclude TC 1 and TC 2
    'PARAM_12_04', [1, 2, 3];    % Exclude TC 1,2,3
};


% -------------------------------------------------------------------------
% 1. LOAD EXPERIMENTAL DATA & FLUID PROPERTIES
% -------------------------------------------------------------------------
% Standardized table header mapping for input spreadsheet
varNames = {'Param', 'T1', 'T2', 'T3', 'T4', 'T5', 'Tsat', 'G', 'x_in', ...
            'q_eff', 'Q_heater', 'delta', 'Re_f', 'SWR', 'mdot_f', 'tau_w', 'Re_v','We_c', 'pi_e'};

% Read main experimental Excel table
T = readtable("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\07_Conferences\ASME SHTC\Data\WallHTModelComparison.xlsx", ...
    'VariableNamingRule', 'preserve');   
T.Properties.VariableNames = varNames;

% Load fluid property structure and measured wake width Master structure
load("testfluidprops.mat")
load("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\02_MFVAL\Data Processing\Jason_FT_Processing_Program\film-thickness-master\Scripts\wake_width_datastructure.mat", "Master")

% Extract key fluid thermophysical properties
k_film  = fluid.KF;         % [W/m-K] Thermal conductivity of saturated liquid
k_vapor = fluid.KG;         % [W/m-K] Thermal conductivity of saturated vapor
rho_f   = fluid.RHOF;       % [kg/m^3] Density of saturated liquid
mu_f    = fluid.MUF;        % [Pa-s] Dynamic viscosity of saturated liquid
nu_f    = mu_f/rho_f;       % [m^2/s] Kinematic viscosity of saturated liquid
PR_f    = fluid.PRANDTLF;   % [-] Liquid Prandtl number
PR_v    = fluid.PRANDTLG;   % [-] Vapor Prandtl number
h_fg    = fluid.HFG;        % [J/kg] Latent heat of vaporization

% Test Section & Obstruction Geometry Dimensions
W     = 0.036;              % [m] Full channel width
H     = 0.012;              % [m] Full channel height
per   = 2*W + 2*H;          % [m] Wetted channel perimeter
A_c   = W*H;                % [m^2] Cross-sectional flow area
D_h   = 4*A_c/per;          % [m] Channel hydraulic diameter
L     = 0.525;              % [m] Total test channel length
L_obs = 0.42;               % [m] Distance from channel INLET to obstruction pin
D_obs = 0.0025;             % [m] Obstruction pin diameter
x_obs = L - L_obs;          % [m] Distance from obstruction pin to channel EXIT (0.105 m)

% Base Liquid Film Thickness via Law of the Wall (y+ = 5 viscous sublayer assumption)
yplus        = 5;                                                           
T.u_star     = (T.tau_w / rho_f).^(0.5);                                    % [m/s] Friction velocity from wall shear stress
T.base_thick = yplus * (T.u_star).^(-1) * nu_f;                            % [m] Calculated base film thickness

% Baseline Convective Heat Transfer Coefficients (HTC)
T.htc_wet    = k_film * (T.base_thick).^(-1);                               % [W/m^2-K] Wet wall liquid film HTC
T.htc_dry    = k_vapor / D_h * 0.023 * (T.Re_v).^0.8 * PR_v^0.4;            % [W/m^2-K] Dry wake vapor core HTC (Dittus-Boelter)

% -------------------------------------------------------------------------
% 2. PARSE MEASURED DRY WAKE WIDTHS & APPLY OUTLIER FILTERS
% -------------------------------------------------------------------------
T.dry_wake_W = cell(height(T), 1);
T.w_free     = cell(height(T), 1);

for i = 1:height(T)
    param_str  = string(T.Param(i));
    param_norm = upper(replace(param_str, "_", "")); % Normalize string (e.g., 'PARAM_05_04' -> 'PARAM0504')
    
    found_match = false;
    measured_widths = zeros(1, 5);
    
    % Search nested Master struct for image filenames matching current parameter name
    for e = 1:length(Master.datasets(2).entries)
        wakes = Master.datasets(2).entries(e).wakes;
        for w = 1:length(wakes)
            if contains(string(wakes(w).imgFileName), param_str)
                % Extract left + right wake half-widths for TC1-TC5 (indices 2 through 6)
                raw_widths = wakes(w).widths(2:6, 1) + wakes(w).widths(2:6, 2);
                measured_widths = raw_widths(:)' * 1e-3; % Convert mm to meters
                found_match = true;
                break;
            end
        end
        if found_match; break; end
    end
    
    if ~found_match
        warning('No wake data found in Master structure for parameter: %s', param_str);
    end
    
    % --- APPLY OUTLIER FILTERING RULES ---
    if filter_outliers
        % Rule 1: Mask zero or negative wake width measurements (no dry wake present)
        measured_widths(measured_widths <= 0) = NaN;
        
        % Rule 2: Mask entire parameter run if listed in full_exclude_params
        norm_full_exclude = upper(replace(full_exclude_params, "_", ""));
        if any(strcmp(param_norm, norm_full_exclude))
            measured_widths(:) = NaN;
        end
        
        % Rule 3: Mask specific TC positions defined in specific_tc_exclude array
        for k = 1:size(specific_tc_exclude, 1)
            target_param = upper(replace(string(specific_tc_exclude{k, 1}), "_", ""));
            if strcmp(param_norm, target_param)
                bad_tcs = specific_tc_exclude{k, 2};
                measured_widths(bad_tcs) = NaN;
            end
        end
    end
    
    T.dry_wake_W{i} = measured_widths;
    T.w_free{i}     = W - measured_widths; % W_free automatically becomes NaN where wake width is NaN
end

% -------------------------------------------------------------------------
% 3. CALCULATE EMPIRICAL PREDICTED WAKE WIDTHS (CLOSURE RELATION)
% -------------------------------------------------------------------------
x_tc = [0.04, 0.05, 0.06, 0.07, 0.08]; % [m] Downstream axial distances from obstruction for TC1-TC5
T.dry_wake_W_calc = cell(height(T), 1);

% Compute adiabatic wake width using empirical correlation [m]
T.W_adiabatic = (1.66 * T.Re_f.^0.085 .* T.SWR.^0.22) * 2 / 1000;
T.phi         = 0.53 * ones(height(T), 1);                          % Partition factor

% Wake region mass flow and Boiling number scaling
T.mdot_wake = T.mdot_f * (D_obs / per) .* T.phi;
scalarPart  = (T.q_eff .* T.W_adiabatic) ./ (T.mdot_wake .* h_fg);

for i = 1:height(T)
    Bo_tc           = scalarPart(i) .* x_tc;                            
    wake_ratio_func = -5 + 0.0058*Bo_tc + (6.3325 * Bo_tc.^1.6617) ./ (0.1719^1.6617 + Bo_tc.^1.6617);
    wake_ratio      = max(0, wake_ratio_func);                          
    
    calc_wakes = T.W_adiabatic(i) .* wake_ratio;
    
    % Mask calculated wakes at positions where measured wakes are NaN for fair comparison
    if filter_outliers
        calc_wakes(isnan(T.dry_wake_W{i})) = NaN;
    end
    T.dry_wake_W_calc{i} = calc_wakes;
end

% -------------------------------------------------------------------------
% 4. SOLVE BASELINE HEAT TRANSFER MODEL AT DOWNSTREAM TC LOCATIONS
% -------------------------------------------------------------------------
% Domain discretization & boundary parameters
M       = 125;              % Number of spanwise wall nodes
N       = 125;              % Number of wall-normal thickness nodes
W_HT    = W/2;              % [m] Half-width of thermal conduction domain (symmetry line)
H_wall  = 0.003;            % [m] Substrate glass wall thickness
k_wall  = 1.2;              % [W/m-K] Thermal conductivity of glass substrate
htc_amb = 10;               % [W/m^2-K] Ambient natural convection HTC
T_inf   = fluid.TSAT;       % [K] Fluid saturation temperature
T_amb   = 22.5 + 273.15;    % [K] Ambient environment temperature
N_rows  = height(T);

T.T_center = cell(N_rows, 1);

% Solve 2D conduction solver for each test condition
for i = 1:N_rows
    w_free_i     = T.w_free{i};          
    N_wfree      = numel(w_free_i);
    T_center_row = zeros(1, N_wfree);
    
    for j = 1:N_wfree
        if isnan(w_free_i(j))
            T_center_row(j) = NaN; % Bypass computationally heavy solver if point is filtered
        else
            [~, ~, ~, T_center_row(j), ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
                M, N, W_HT, H_wall, k_wall, x_obs, ...
                T.q_eff(i), w_free_i(j), T.htc_dry(i), T.htc_wet(i), ...
                htc_amb, T_inf, T_amb);
        end
    end
    T.T_center{i} = T_center_row;
end

% Construct structured evaluation table for comparative plotting
Temps = table();
Temps.ParamName = string(T.Param);

T_meas_mat = [double(T.T1), double(T.T2), double(T.T3), double(T.T4), double(T.T5)];
W_meas_mat = vertcat(T.dry_wake_W{:});

% Mask measured wall temperatures for filtered wake points
T_meas_mat(isnan(W_meas_mat)) = NaN;

Temps.T1_meas = T_meas_mat(:, 1);
Temps.T2_meas = T_meas_mat(:, 2);
Temps.T3_meas = T_meas_mat(:, 3);
Temps.T4_meas = T_meas_mat(:, 4);
Temps.T5_meas = T_meas_mat(:, 5);

T_calc_mat    = vertcat(T.T_center{:});
Temps.T1_calc = T_calc_mat(:, 1);
Temps.T2_calc = T_calc_mat(:, 2);
Temps.T3_calc = T_calc_mat(:, 3);
Temps.T4_calc = T_calc_mat(:, 4);
Temps.T5_calc = T_calc_mat(:, 5);

disp('Data pre-processing and baseline solution complete.');

% -------------------------------------------------------------------------
% 5. EXECUTE ALL PLOTS, SENSITIVITY, AND OPTIMIZATION
% -------------------------------------------------------------------------
% Temperature Parity, Residuals, & Spatial Error
plotParity(Temps);
%plotResiduals(Temps);
%plotSpatialError(Temps);

% Spatial Thermal Profile Comparisons
%plotProfileComparison(Temps, {'PARAM01_04', 'PARAM04_03', 'PARAM05_03'});
plotAllProfilesSummary(Temps);

% Wake Width Parity & Spatial Profiles
plotWakeWidthParity(T);
%plotWakeWidthProfileComparison(T, {'PARAM01_04', 'PARAM04_03', 'PARAM05_03'});
%plotAllWakeWidthProfilesSummary(T);

% Cross-Variable Scatter Plots (Temperature vs Wake Width)
%plotWakeWidthVsTempMeas(T, Temps);
%plotWakeWidthVsTempCalc(T, Temps);
%plotWakeWidthCalcVsTempCalc(T, Temps);

% Perturbation Sensitivity Analysis (+/- 10%)
runSensitivityAnalysis(T, Temps, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb);

% Inverse HTC Optimization & Multi-Panel Correlation Analysis
T = optimizeDryHTC(T, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb);

plotHTCParametricCorrelations(T, fluid, D_h);

% Droplet Entrainment Modeling & Correlation Analysis
T = computeEntrainedFraction(T, fluid, A_c);
plotEntrainmentCorrelations(T);


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
    
    TC_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    % Omit NaNs resulting from outlier filtering
    valid = ~isnan(T_meas) & ~isnan(T_calc);
    T_meas    = T_meas(valid);
    T_calc    = T_calc(valid);
    TC_labels = TC_labels(valid);
    
    figure('Name', 'Model Parity Analysis', 'Color', 'w', 'Position', [100, 100, 650, 550]);
    hold on;
    
    lims = [min([T_meas; T_calc]) - 2, max([T_meas; T_calc]) + 2];
    h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5); 
    
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(T_meas, T_calc, TC_labels, colors, 'o', 6, 'on');
    
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
% PLOT: Residual Error Histogram & Gaussian PDF Fit
% -------------------------------------------------------------------------
function plotResiduals(Temps)
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];
    residuals = T_calc - T_meas;
    residuals = residuals(~isnan(residuals)); % Strip NaNs
    
    figure('Name', 'Model Residuals', 'Color', 'w', 'Position', [150, 150, 650, 500]);
    histogram(residuals, 'Normalization', 'pdf', 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'w', 'FaceAlpha', 0.8);
    hold on;
    
    pd = fitdist(residuals, 'Normal');
    x_val = linspace(min(residuals)-2, max(residuals)+2, 250);
    plot(x_val, pdf(pd, x_val), 'r-', 'LineWidth', 2);
    
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
% PLOT: Spatial Boxplot Distribution of Residual Errors
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
% PLOT: Spatial Thermal Profiles for Selected Runs
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
% PLOT: Aggregate Thermal Profiles (Ensemble Mean + Individual Traces)
% -------------------------------------------------------------------------
function plotAllProfilesSummary(Temps)
    T_meas_all = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas];
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc];
    
    figure('Name', 'Aggregate Profile Comparison', 'Color', 'w', 'Position', [250, 250, 800, 500]);
    hold on;
    
    for i = 1:height(Temps)
        plot(1:5, T_meas_all(i,:), '-', 'Color', [0 0.447 0.741 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
        plot(1:5, T_calc_all(i,:), '-', 'Color', [0.850 0.325 0.098 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
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
% PLOT: Dry Wake Width Parity (Measured vs Calculated)
% -------------------------------------------------------------------------
function plotWakeWidthParity(T)
    meas_all  = [];
    calc_all  = [];
    tc_labels = [];
    
    for i = 1:height(T)
        meas_i = T.dry_wake_W{i};
        calc_i = T.dry_wake_W_calc{i};
        
        if ~isempty(meas_i) && ~isempty(calc_i)
            meas_all  = [meas_all, meas_i];
            calc_all  = [calc_all, calc_i];
            tc_labels = [tc_labels, 1:5];
        end
    end
    
    valid     = ~isnan(meas_all) & ~isnan(calc_all);
    meas_all  = meas_all(valid);
    calc_all  = calc_all(valid);
    tc_labels = tc_labels(valid);
    
    figure('Name', 'Wake Width Parity Analysis', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    lims = [min([meas_all, calc_all])*0.9, max([meas_all, calc_all])*1.1];
    if isempty(lims) || min(lims) == 0; lims = [0, 0.02]; end 
    
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
% PLOT: Aggregate Wake Width Spatial Profiles
% -------------------------------------------------------------------------
function plotAllWakeWidthProfilesSummary(T)
    W_meas_all = vertcat(T.dry_wake_W{:});
    W_calc_all = vertcat(T.dry_wake_W_calc{:});
    
    figure('Name', 'Aggregate Wake Width Profile Comparison', 'Color', 'w', 'Position', [250, 250, 800, 500]);
    hold on;
    
    for i = 1:height(T)
        if ~isempty(T.dry_wake_W{i}) && ~isempty(T.dry_wake_W_calc{i})
            plot(1:5, T.dry_wake_W{i}, '-', 'Color', [0 0.447 0.741 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
            plot(1:5, T.dry_wake_W_calc{i}, '-', 'Color', [0.850 0.325 0.098 0.15], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    
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
% PLOT: Measured Wake Width vs Measured Temperature
% -------------------------------------------------------------------------
function plotWakeWidthVsTempMeas(T, Temps)
    W_meas_all = vertcat(T.dry_wake_W{:}); 
    T_meas_all = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas]; 
    
    W_flat = W_meas_all(:);
    T_flat = T_meas_all(:);
    
    numRows   = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
             
    valid  = ~isnan(W_flat) & ~isnan(T_flat);
    W_flat = W_flat(valid);
    T_flat = T_flat(valid);
    tc_labels = tc_labels(valid);
    
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
% PLOT: Measured Wake Width vs Calculated Temperature
% -------------------------------------------------------------------------
function plotWakeWidthVsTempCalc(T, Temps)
    W_meas_all = vertcat(T.dry_wake_W{:}); 
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc]; 
    
    W_flat = W_meas_all(:);
    T_flat = T_calc_all(:);
    
    numRows   = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
             
    valid  = ~isnan(W_flat) & ~isnan(T_flat);
    W_flat = W_flat(valid);
    T_flat = T_flat(valid);
    tc_labels = tc_labels(valid);
    
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
% PLOT: Calculated Wake Width vs Calculated Temperature
% -------------------------------------------------------------------------
function plotWakeWidthCalcVsTempCalc(T, Temps)
    W_calc_all = vertcat(T.dry_wake_W_calc{:}); 
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc]; 
    
    W_flat = W_calc_all(:);
    T_flat = T_calc_all(:);
    
    numRows   = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
             
    valid  = ~isnan(W_flat) & ~isnan(T_flat);
    W_flat = W_flat(valid);
    T_flat = T_flat(valid);
    tc_labels = tc_labels(valid);
    
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
% SENSITIVITY ANALYSIS: Input Perturbations (+/-10% on h_wet, h_dry, W_meas)
% -------------------------------------------------------------------------
function runSensitivityAnalysis(T, Temps, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb)

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    N_rows     = height(T);
    T_meas_mat = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas];
    
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
            
            for j = 1:5
                if isnan(w_free_i(j))
                    T_calc_mat(i, j) = NaN;
                else
                    [~, ~, ~, T_center, ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
                        M, N, W_HT, H_wall, k_wall, x_obs, ...
                        T.q_eff(i), w_free_i(j), htc_dry_i, htc_wet_i, ...
                        htc_amb, T_inf, T_amb);
                    T_calc_mat(i, j) = T_center;
                end
            end
        end
        
        err      = T_calc_mat - T_meas_mat;
        rmse_val = sqrt(mean(err(:).^2, 'omitnan'));                     
        mae_val  = mean(abs(err(:)), 'omitnan');                         
        mbe_val  = mean(err(:), 'omitnan');                              
        
        results(s).name   = sName;
        results(s).T_calc = T_calc_mat;
        results(s).RMSE   = rmse_val;
        results(s).MAE    = mae_val;
        results(s).MBE    = mbe_val;
        
        fprintf('%-14s | RMSE: %6.2f °C | MAE: %6.2f °C | MBE: %6.2f °C\n', ...
            sName, rmse_val, mae_val, mbe_val);
    end
    fprintf('=========================================================\n\n');
    
    % --- PLOT 1: Individual Parameter Parity Sensitivity Envelopes ---
    figure('Name', 'Sensitivity Parity Bounds - Individual Parameters', ...
           'Color', 'w', 'Position', [100, 150, 1350, 430]);
    t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
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
    
    % --- PLOT 2: RMSE Impact Bar Chart ---
    figure('Name', 'Sensitivity RMSE Impact', 'Color', 'w', 'Position', [200, 200, 750, 450]);
    names = {results.name};
    rmses = [results.RMSE];
    
    b = bar(rmses, 'FaceColor', 'flat');
    b.CData(1, :)   = [0 0.45 0.74];          
    b.CData(2:3, :) = repmat([0.85 0.33 0.10], 2, 1); 
    b.CData(4:5, :) = repmat([0.93 0.69 0.13], 2, 1); 
    b.CData(6:7, :) = repmat([0.49 0.18 0.56], 2, 1); 
    
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
% OPTIMIZATION: Row-by-Row Inverse Solver (fminbnd) for Dry HTC (h_dry)
% -------------------------------------------------------------------------
function T = optimizeDryHTC(T, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb)
    N_rows          = height(T);
    T.htc_dry_opt   = zeros(N_rows, 1);
    T.RMSE_htc_opt  = zeros(N_rows, 1);
    
    h_min   = 5;       % [W/m^2-K] Lower search bound
    h_max   = 10000;    % [W/m^2-K] Upper search bound
    options = optimset('Display', 'off', 'TolX', 1e-1); 
    
    fprintf('\n=========================================================\n');
    fprintf('  RUNNING INVERSE HTC OPTIMIZATION (ROW-BY-ROW)\n');
    fprintf('=========================================================\n');
    
    for i = 1:N_rows
        T_meas_i = [T.T1(i), T.T2(i), T.T3(i), T.T4(i), T.T5(i)];
        w_free_i = T.w_free{i};
        
        valid_tc = ~isnan(T_meas_i) & ~isnan(w_free_i);
        if ~any(valid_tc)
            T.htc_dry_opt(i)  = NaN;
            T.RMSE_htc_opt(i) = NaN;
            continue;
        end
        
        objFun = @(h_dry) calcRowSSE(h_dry, T_meas_i, w_free_i, T.q_eff(i), T.htc_wet(i), ...
            M, N, W_HT, H_wall, k_wall, x_obs, htc_amb, T_inf, T_amb);
        
        [h_opt, sse_val] = fminbnd(objFun, h_min, h_max, options);
        
        T.htc_dry_opt(i)  = h_opt;
        T.RMSE_htc_opt(i) = sqrt(sse_val / sum(valid_tc)); 
        
        if mod(i, 10) == 0 || i == N_rows
            fprintf('  Processed %d/%d test conditions...\n', i, N_rows);
        end
    end
    fprintf('Inverse HTC Optimization Complete!\n');
    fprintf('=========================================================\n\n');
end

% -------------------------------------------------------------------------
% OBJECTIVE FUNCTION: Minimize Sum of Squared Errors across valid TCs
% -------------------------------------------------------------------------
function sse = calcRowSSE(h_dry, T_meas, w_free, q_eff, htc_wet, ...
    M, N, W_HT, H_wall, k_wall, x_obs, htc_amb, T_inf, T_amb)

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    valid_tc = ~isnan(T_meas) & ~isnan(w_free);
    T_calc   = zeros(1, 5);
    
    for j = find(valid_tc)
        [~, ~, ~, T_center, ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
            M, N, W_HT, H_wall, k_wall, x_obs, ...
            q_eff, w_free(j), h_dry, htc_wet, ...
            htc_amb, T_inf, T_amb);
        T_calc(j) = T_center;
    end
    sse = sum((T_calc(valid_tc) - T_meas(valid_tc)).^2); 
end

% -------------------------------------------------------------------------
% PLOT: Multi-Panel Parametric Correlations for Backed-Out Dry HTC
% -------------------------------------------------------------------------
function plotHTCParametricCorrelations(T, fluid, D_h)
    figure('Name', 'Optimized HTC Parametric Correlations', 'Color', 'w', 'Position', [100, 100, 1200, 750]);
    t = tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Per-panel distinct color scheme for optimized h_dry
    c1 = [0.00 0.45 0.74];  % Blue (Re_v)
    c2 = [0.49 0.18 0.56];  % Purple (SWR)
    c3 = [0.85 0.33 0.10];  % Orange (q_eff)
    c4 = [0.47 0.67 0.19];  % Green (x_in)
    c5 = [0.93 0.69 0.13];  % Gold/Yellow (Re_f)
    c6 = [0.30 0.75 0.93];  % Cyan/Teal (delta)
    
    % Consistent baseline color across all subplots
    c_db = [0.40 0.40 0.40]; % Dark Gray for Dittus-Boelter Baseline
    
    % Mask out NaN rows resulting from outlier filtering
    valid = ~isnan(T.htc_dry_opt);
    
    % 1. h_dry vs Re_v
    nexttile;
    semilogx(T.Re_v(valid), T.htc_dry_opt(valid), 'o', 'Color', c1, 'MarkerFaceColor', c1, 'DisplayName', 'Optimized $h_{\mathrm{dry}}$');
    hold on;
    semilogx(T.Re_v(valid), T.htc_dry(valid), 's', 'Color', c_db, 'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Dittus-Boelter Baseline');
    xlabel('Vapor Reynolds Number, $\mathrm{Re}_v$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{1. vs. Vapor Reynolds Number ($\mathrm{Re}_v$)}', 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    % 2. h_dry vs Core Weber Number
    nexttile;
    plot(T.We_c(valid), T.htc_dry_opt(valid), 'o', 'Color', c2, 'MarkerFaceColor', c2, 'DisplayName', 'Optimized $h_{\mathrm{dry}}$');
    hold on;
    plot(T.We_c(valid), T.htc_dry(valid), 's', 'Color', c_db, 'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Dittus-Boelter Baseline');
    xlabel('Core Weber Number, $\mathrm{We_c}$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{2. vs. Core Weber Number}', 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    % 3. h_dry vs Effective Heat Flux (q_eff)
    nexttile;
    plot(T.q_eff(valid) / 1000, T.htc_dry_opt(valid), 'o', 'Color', c3, 'MarkerFaceColor', c3, 'DisplayName', 'Optimized $h_{\mathrm{dry}}$');
    hold on;
    plot(T.q_eff(valid) / 1000, T.htc_dry(valid), 's', 'Color', c_db, 'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Dittus-Boelter Baseline');
    xlabel('Effective Heat Flux, $q''_{\mathrm{eff}}$ ($\mathrm{kW/m}^2$)', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{3. vs. Heat Flux $q''_{\mathrm{eff}}$}', 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    % 4. h_dry vs Inlet Vapor Quality (x_in)
    nexttile;
    plot(T.x_in(valid), T.htc_dry_opt(valid), 'o', 'Color', c4, 'MarkerFaceColor', c4, 'DisplayName', 'Optimized $h_{\mathrm{dry}}$');
    hold on;
    plot(T.x_in(valid), T.htc_dry(valid), 's', 'Color', c_db, 'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Dittus-Boelter Baseline');
    xlabel('Inlet Quality, $x_{\mathrm{in}}$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{4. vs. Inlet Quality $x_{\mathrm{in}}$}', 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    % 5. h_dry vs Non Dimensional Entrainment Number (Okawa 2003 publication)
    nexttile;
    plot(T.pi_e(valid), T.htc_dry_opt(valid), 'o', 'Color', c5, 'MarkerFaceColor', c5, 'DisplayName', 'Optimized $h_{\mathrm{dry}}$');
    hold on;
    plot(T.pi_e(valid), T.htc_dry(valid), 's', 'Color', c_db, 'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Dittus-Boelter Baseline');
    xlabel('Okawa Entrainment Number, $\mathrm{pi}_e$', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{5. vs. Okawa Entrainment Number $\mathrm{pi}_e$}', 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    % 6. h_dry vs Pre-Obstruction Measured Film Thickness (delta)
    nexttile;
    plot(T.delta(valid), T.htc_dry_opt(valid), 'o', 'Color', c6, 'MarkerFaceColor', c6, 'DisplayName', 'Optimized $h_{\mathrm{dry}}$');
    hold on;
    plot(T.delta(valid), T.htc_dry(valid), 's', 'Color', c_db, 'LineWidth', 1.2, 'MarkerFaceColor', 'none', 'DisplayName', 'Dittus-Boelter Baseline');
    xlabel('Pre-Obstruction Measured Film Thickness, $\delta$ (m)', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('$h_{\mathrm{dry}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{6. vs. Pre-Obstruction Film Thickness $\delta$}', 'Interpreter', 'latex', 'FontSize', 12);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 8);
    grid on;
    
    title(t, '\textbf{Parametric Scaling Analysis of Backed-Out Dry Heat Transfer Coefficient}', ...
        'Interpreter', 'latex', 'FontSize', 14);
end

% -------------------------------------------------------------------------
% ENTRAINMENT: Energy Balance Calculation for Droplet Fraction E
% -------------------------------------------------------------------------
function T = computeEntrainedFraction(T, fluid, A_c)
    h_fg       = fluid.HFG;
    mdot_total = T.G * A_c; % [kg/s] Total mass flow rate
    
    % Local thermodynamic vapor quality at obstruction using measured Q_heater [W]
    x_local     = T.x_in + T.Q_heater ./ (mdot_total .* h_fg);
    mdot_liquid = mdot_total .* (1 - x_local);                  % [kg/s] Liquid flow rate
    mdot_drop   = max(0, mdot_liquid - T.mdot_f);               % [kg/s] Entrained droplet flow rate
    
    T.x_local     = x_local;
    T.mdot_liquid = mdot_liquid;
    T.mdot_drop   = mdot_drop;
    T.E           = mdot_drop ./ mdot_liquid;                   % [-] Entrained fraction
end

% -------------------------------------------------------------------------
% PLOT: Droplet Entrainment Correlations against Dry HTC & Ratio
% -------------------------------------------------------------------------
function plotEntrainmentCorrelations(T)
    figure('Name', 'HTC Droplet Entrainment Correlations', 'Color', 'w', 'Position', [150, 150, 1000, 450]);
    t = tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    valid = ~isnan(T.pi_e) & ~isinf(T.pi_e) & (T.pi_e >= 0) & ~isnan(T.htc_dry_opt);
    
    % Panel 1: Optimized h_dry vs Entrained Fraction E
    nexttile;
    plot(T.pi_e(valid), T.htc_dry_opt(valid), 'o', 'Color', [0 0.45 0.74], ...
        'MarkerFaceColor', [0 0.45 0.74], 'MarkerSize', 7);
    hold on;
    
    if any(valid)
        p1 = polyfit(T.pi_e(valid), T.htc_dry_opt(valid), 1);
        pi_e_fit = linspace(min(T.pi_e(valid)), max(T.pi_e(valid)), 100);
        plot(pi_e_fit, polyval(p1, pi_e_fit), 'k--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Linear Fit (Slope: %.1f)', p1(1)));
        legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
    end
    
    xlabel('Okawa Entrainment Number $\mathrm{pi}_e$}', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('Optimized Dry HTC, $h_{\mathrm{dry, opt}}$ ($\mathrm{W/m}^2\mathrm{K}$)', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{1. $h_{\mathrm{dry, opt}}$ vs. Okawa Entrainment Number $\mathrm{pi}_e$', 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    grid on;
    
    % Panel 2: HTC Ratio (h_dry_opt / h_wet) vs Entrained Fraction E
    nexttile;
    ratio_h = T.htc_dry_opt ./ T.htc_dry;
    plot(T.pi_e(valid), ratio_h(valid), 's', 'Color', [0.85 0.33 0.10], ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 7);
    hold on;
    
    if any(valid)
        p2 = polyfit(T.pi_e(valid), ratio_h(valid), 1);
        plot(pi_e_fit, polyval(p2, pi_e_fit), 'k--', 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Linear Fit (Slope: %.2f)', p2(1)));
        legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 9);
    end
    
    xlabel('Okawa Entrainment Number $\mathrm{pi}_e$}', 'Interpreter', 'latex', 'FontSize', 11);
    ylabel('HTC Ratio, $h_{\mathrm{dry, opt}} / h_{\mathrm{dry}}$', 'Interpreter', 'latex', 'FontSize', 11);
    title('\textbf{2. HTC Ratio ($h_{\mathrm{dry, opt}} / h_{\mathrm{wet}}$) vs. Okawa Entrainment Number $\mathrm{pi}_e$', 'Interpreter', 'latex', 'FontSize', 12);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 10);
    grid on;
    
    title(t, '\textbf{Impact of Liquid Droplet Entrainment on Dry Wake Heat Transfer}', ...
        'Interpreter', 'latex', 'FontSize', 14);
end