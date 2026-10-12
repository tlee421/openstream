% --HTModelSensitivity.m--%
import Solvers.Obstruction.WallHT.ftoglass_ann_optim

% Load in data regarding known values:
varNames = {'Param', 'T1', 'T2', 'T3', 'T4', 'T5', 'Tsat', 'G', 'x_in', ...
            'q_eff', 'Q_heater', 'delta', 'Re_f', 'SWR', 'mdot_f', 'tau_w', 'Re_v'};
T = readtable("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\07_Conferences\ASME SHTC\Data\WallHTModelComparison.xlsx", ...
    'VariableNamingRule', 'preserve');   
T.Properties.VariableNames = varNames;   

% Load fluid properties and measured wake data structure
load("testfluidprops.mat")
load("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\02_MFVAL\Data Processing\Jason_FT_Processing_Program\film-thickness-master\Scripts\wake_width_datastructure.mat", "Master")

k_film  = fluid.KF;                                                         
k_vapor = fluid.KG;                                                        
rho_f   = fluid.RHOF;
mu_f    = fluid.MUF;
nu_f    = mu_f/rho_f;
PR_f    = fluid.PRANDTLF;
PR_v    = fluid.PRANDTLG;
h_fg    = fluid.HFG;

% Wall geometry
W     = 0.036;
H     = 0.012;
per   = 2*W + 2*H;
A_c   = W*H;
D_h   = 4*A_c/per;
L     = 0.525;
L_obs = 0.42; 
D_obs = 0.0025;
x_obs = L - L_obs;                                                          

% Base film thickness and HTCs
yplus        = 5;                                                                 
T.u_star     = (T.tau_w/rho_f).^(0.5);
T.base_thick = yplus*(T.u_star).^(-1)*nu_f;                                
T.htc_wet    = k_film * (T.base_thick).^(-1);
T.htc_dry    = k_vapor/D_h * 0.023 * (T.Re_v).^0.8 * PR_v^0.4;                

% -------------------------------------------------------------------------
% PARSE MEASURED DRY WAKE WIDTHS
% -------------------------------------------------------------------------
T.dry_wake_W = cell(height(T), 1);
T.w_free     = cell(height(T), 1);

for i = 1:height(T)
    param_str   = string(T.Param(i));
    found_match = false;
    
    for e = 1:length(Master.datasets(2).entries)
        wakes = Master.datasets(2).entries(e).wakes;
        for w = 1:length(wakes)
            if contains(string(wakes(w).imgFileName), param_str)
                % Extract widths for Thermocouples 1-5 (convert mm to m)
                measured_widths = wakes(w).widths(2:6, 1) + wakes(w).widths(2:6, 2);
                T.dry_wake_W{i} = measured_widths(:)' * 1e-3;                
                T.w_free{i}     = W - T.dry_wake_W{i};
                found_match     = true;
                break;
            end
        end
        if found_match; break; end
    end
    
    if ~found_match
        warning('No wake data found in Master structure for parameter: %s', param_str);
        T.dry_wake_W{i} = zeros(1, 5);
        T.w_free{i}     = W * ones(1, 5);
    end
end

% -------------------------------------------------------------------------
% SOLVE BASELINE HEAT TRANSFER MODEL
% -------------------------------------------------------------------------
M       = 125;                                                                    
N       = 125;                                                                   
W_HT    = W/2;                                                                
H_wall  = 0.003;                                                                 
k_wall  = 1.2;                                                              
htc_amb = 10;
T_inf   = fluid.TSAT;
T_amb   = 22.5 + 273.15;
N_rows  = height(T);

T.T_center = cell(N_rows, 1);

for i = 1:N_rows
    w_free_i  = T.w_free{i};          
    N_wfree   = numel(w_free_i);
    T_center_row = zeros(1, N_wfree);
    
    for j = 1:N_wfree
        [~, ~, ~, T_center_row(j), ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
            M, N, W_HT, H_wall, k_wall, x_obs, ...
            T.q_eff(i), w_free_i(j), T.htc_dry(i), T.htc_wet(i), ...
            htc_amb, T_inf, T_amb);
    end
    T.T_center{i} = T_center_row;
end

% Build Temps evaluation table
Temps = table();
Temps.ParamName = string(T.Param);
Temps.T1_meas   = double(T.T1);
Temps.T2_meas   = double(T.T2);
Temps.T3_meas   = double(T.T3);
Temps.T4_meas   = double(T.T4);
Temps.T5_meas   = double(T.T5);

T_calc_mat      = vertcat(T.T_center{:});
Temps.T1_calc   = T_calc_mat(:, 1);
Temps.T2_calc   = T_calc_mat(:, 2);
Temps.T3_calc   = T_calc_mat(:, 3);
Temps.T4_calc   = T_calc_mat(:, 4);
Temps.T5_calc   = T_calc_mat(:, 5);

%% Execute Baseline Plots & Sensitivity Analysis
plotParity(Temps);
runSensitivityAnalysis(T, Temps, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb);

%% =========================================================================
% PLOTTING & SENSITIVITY FUNCTIONS
% =========================================================================

function plotParity(Temps)
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];
    numRows = height(Temps);
    TC_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
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

function runSensitivityAnalysis(T, Temps, M, N, W_HT, H_wall, k_wall, x_obs, W, htc_amb, T_inf, T_amb)

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    N_rows     = height(T);
    T_meas_mat = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas];
    
    % Scenarios: {Scenario Name, h_wet factor, h_dry factor, W_meas factor}
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
                [~, ~, ~, T_center, ~, ~, ~, ~, ~, ~] = ftoglass_ann_optim( ...
                    M, N, W_HT, H_wall, k_wall, x_obs, ...
                    T.q_eff(i), w_free_i(j), htc_dry_i, htc_wet_i, ...
                    htc_amb, T_inf, T_amb);
                T_calc_mat(i, j) = T_center;
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
    
    % -------------------------------------------------------------------------
    % PLOT 1: INDIVIDUAL PARITY SENSITIVITY PLOTS (1x3 TILED LAYOUT)
    % -------------------------------------------------------------------------
    figure('Name', 'Sensitivity Parity Bounds - Individual Parameters', ...
           'Color', 'w', 'Position', [100, 150, 1350, 430]);
    t = tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    % Configuration matrix: {Subplot Title, Low Index, High Index, Line Color}
    param_configs = {
        '$h_{\mathrm{wet}} \pm 10\%$ Sensitivity',     3, 2, [0.85 0.33 0.10];
        '$h_{\mathrm{dry}} \pm 10\%$ Sensitivity',     5, 4, [0.93 0.69 0.13];
        '$W_{\mathrm{meas}} \pm 10\%$ Sensitivity',    7, 6, [0.49 0.18 0.56];
    };
    
    T_meas_flat = T_meas_mat(:);
    lims        = [min(T_meas_flat) - 2, max(T_meas_flat) + 20];
    
    for p = 1:3
        nexttile;
        hold on;
        
        pTitle   = param_configs{p, 1};
        idx_low  = param_configs{p, 2};
        idx_high = param_configs{p, 3};
        pColor   = param_configs{p, 4};
        
        % 1:1 Parity Line
        h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Agreement');
        
        % Sensitivity range lines for this specific parameter
        plot([T_meas_flat, T_meas_flat]', [results(idx_low).T_calc(:), results(idx_high).T_calc(:)]', ...
            '-', 'Color', [pColor, 0.40], 'LineWidth', 1.5, 'HandleVisibility', 'off');
        
        % Overlay Baseline points
        h_base = scatter(T_meas_flat, results(1).T_calc(:), 25, [0 0.45 0.74], 'filled', 'DisplayName', 'Baseline');
        
        % Dummy line for legend entry
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
    
    % -------------------------------------------------------------------------
    % PLOT 2: RMSE BAR CHART
    % -------------------------------------------------------------------------
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