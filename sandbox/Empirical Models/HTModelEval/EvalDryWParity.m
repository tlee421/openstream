% --EvalDryWParity.m--%
import Solvers.Obstruction.WallHT.ftoglass_ann_optim

% Load in data regarding known values:
varNames = {'Param', 'T1', 'T2', 'T3', 'T4', 'T5', 'Tsat', 'G', 'x_in', ...
            'q_eff', 'Q_heater', 'delta', 'Re_f', 'SWR', 'mdot_f', 'tau_w', 'Re_v'};
T = readtable("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\07_Conferences\ASME SHTC\Data\WallHTModelComparison.xlsx", ...
    'VariableNamingRule', 'preserve');   
T.Properties.VariableNames = varNames;   

% load in fluid properties and structural data
load("testfluidprops.mat")
load("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\02_MFVAL\Data Processing\Jason_FT_Processing_Program\film-thickness-master\Scripts\wake_width_datastructure.mat", "Master") % Ensure this file is on path

k_film = fluid.KF;                                                         
k_vapor = fluid.KG;                                                        
rho_f = fluid.RHOF;
mu_f = fluid.MUF;
nu_f = mu_f/rho_f;
PR_f = fluid.PRANDTLF;
PR_v = fluid.PRANDTLG;
h_fg = fluid.HFG;

% Wall geometry:
W = 0.036;
H = 0.012;
per = 2*W + 2*H;
A_c = W*H;
D_h = 4*A_c/per;
L = 0.525;
L_obs = 0.42; 
D_obs = 0.0025;
x_obs = L - L_obs;                                                          

% Compute base film thickness using law of the wall
yplus = 5;                                                                 
T.u_star = (T.tau_w/rho_f).^(0.5);
T.base_thick = yplus*(T.u_star).^(-1)*nu_f;                                

% Compute wet wall HTC:
T.htc_wet = k_film * (T.base_thick).^(-1);

% Compute dry wall HTC using Dittus Boelter
T.htc_dry = k_vapor/D_h * 0.023 * (T.Re_v).^0.8 * PR_v^0.4;                

% -------------------------------------------------------------------------
% PARSE MEASURED DRY WAKE WIDTHS FROM DATA STRUCTURE
% -------------------------------------------------------------------------
T.dry_wake_W = cell(height(T), 1);
T.w_free = cell(height(T), 1);

for i = 1:height(T)
    param_str = string(T.Param(i));
    found_match = false;
    
    for e = 1:length(Master.datasets(2).entries)
        wakes = Master.datasets(2).entries(e).wakes;
        for w = 1:length(wakes)
            if contains(string(wakes(w).imgFileName), param_str)
                % Extract widths for Thermocouples 1 through 5
                measured_widths = wakes(w).widths(2:6, 1) + wakes(w).widths(2:6, 2);
                T.dry_wake_W{i} = measured_widths(:)' *1e-3;                %convert to [m] from [mm]
                T.w_free{i} = W - T.dry_wake_W{i};
                
                found_match = true;
                break;
            end
        end
        if found_match; break; end
    end
    
    if ~found_match
        warning('No wake data found in Master structure for parameter: %s', param_str);
        T.dry_wake_W{i} = zeros(1, 5);
        T.w_free{i} = W * ones(1, 5);
    end
end

% -------------------------------------------------------------------------
% CALCULATE PREDICTED WAKE WIDTHS (ORIGINAL EMPIRICAL CORRELATION)
% -------------------------------------------------------------------------
x_tc = [0.04, 0.05, 0.06, 0.07, 0.08]; 
T.dry_wake_W_calc = cell(height(T), 1);

T.W_adiabatic = (1.66 * T.Re_f.^0.085 .* T.SWR.^0.22)*2/1000;
T.phi = 0.53 * ones(height(T),1);
%Use phi values from optimizing phi for fitting with measured temps from empiricalModelPhiFinder.m
% T.phi = [
% 0.9999
% 0.9999
% 0.5842
% 0.5516
% 0.9999
% 0.6219
% 0.5658
% 0.6110
% 0.9999
% 0.5235
% 0.5690
% 0.5826
% 0.9999
% 0.6664
% 0.6261
% 0.6917
% 0.9999
% 0.9999
% 1.0000
% 0.0799
% 0.9999
% 0.9999
% 0.6219
% 0.6196
% 0.9999
% 0.6219
% 0.6555
% 0.6656
% 0.9999
% 0.6505
% 0.6876
% 0.7209
% 0.9999
% 0.6526
% 0.7574
% 0.8894
% 0.9999
% 0.4277
% 0.9999
% 0.9999
% 0.9999
% 0.9999
% 0.7167
% 0.7581
% 0.9999
% 0.6219
% 0.7687
% 0.7633
% 0.9999
% 0.7573
% 0.7768
% 0.7729
% 0.9999
% 0.7055
% 0.7878
% 0.8674
% 0.9999
% 0.6724
% 0.9999
% 0.9999
% ];

T.mdot_wake = T.mdot_f * (D_obs/per) .* T.phi;
scalarPart = (T.q_eff .* T.W_adiabatic) ./ (T.mdot_wake .* h_fg);

for i = 1:height(T)
    Bo_tc = scalarPart(i) .* x_tc; 
    wake_ratio_func = -5 + 0.0058*Bo_tc + (6.3325 * Bo_tc.^1.6617)./(0.1719^1.6617 + Bo_tc.^1.6617);
    wake_ratio = max(0, wake_ratio_func);
    T.dry_wake_W_calc{i} = T.W_adiabatic(i) .* wake_ratio;
end

% -------------------------------------------------------------------------
% SOLVE HEAT TRANSFER AT DOWNSTREAM LOCATIONS
% -------------------------------------------------------------------------
M = 125;                                                                    
N = 125;                                                                   
W_HT = W/2;                                                                
H = 0.003;                                                                 
k_wall = 1.2;                                                              
htc_amb = 10;
T_inf = fluid.TSAT;
T_amb = 22.5 + 273.15;
N_rows = height(T);

T.x_full     = cell(N_rows, 1);
T.y          = cell(N_rows, 1);
T.T_full     = cell(N_rows, 1);
T.T_center   = cell(N_rows, 1);
T.q_bot_full = cell(N_rows, 1);
T.q_bot_total= cell(N_rows, 1);
T.q_top_full = cell(N_rows, 1);
T.q_top_total= cell(N_rows, 1);
T.x_out      = cell(N_rows, 1);   
T.M2         = cell(N_rows, 1);

for i = 1:N_rows
    w_free_i = T.w_free(i,:);        
    w_free_i = w_free_i{1};          
    N_wfree = numel(w_free_i);
    
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
            M, N, W_HT, H, k_wall, x_obs, ...
            T.q_eff(i), w_free_i(j), T.htc_dry(i), T.htc_wet(i), ...
            htc_amb, T_inf, T_amb);
    end
    
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

%% Format results for plotting
T.T_calc = cell(height(T), 1);  
for i = 1:height(T)
    T.T_calc{i} = T.T_center{i};  
end

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

% Clean data types
for colIdx = 1:width(Temps)
    colName = Temps.Properties.VariableNames{colIdx};
    if strcmp(colName, 'ParamName')
        if ~iscellstr(Temps.(colName)) && ~isstring(Temps.(colName))
            Temps.(colName) = string(Temps.(colName));
        end
    else
        rawCol = Temps.(colName);
        if iscell(rawCol)
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
        Temps.(colName) = numericCol(:);
    end
end

disp('Data cleaning complete. Column types:');
summary(Temps)

%% Execute Plots
%Temp Parity
plotParity(Temps);
%plotResiduals(Temps);
plotSpatialError(Temps);
plotProfileComparison(Temps, {'PARAM08_02', 'PARAM08_03', 'PARAM08_04'});
plotAllProfilesSummary(Temps);

%Dry wake width parity
plotWakeWidthParity(T);
% Execute spatial wake width profile plots
plotWakeWidthProfileComparison(T, {'PARAM01_04', 'PARAM04_03', 'PARAM05_03'});
plotAllWakeWidthProfilesSummary(T);

%%
%Plot dry wake width vs temperature
plotWakeWidthVsTempMeas(T, Temps)
plotWakeWidthVsTempCalc(T, Temps)
plotWakeWidthCalcVsTempCalc(T, Temps)
%%

%% Plotting functions
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

function plotResiduals(Temps)
    T_meas = [Temps.T1_meas; Temps.T2_meas; Temps.T3_meas; Temps.T4_meas; Temps.T5_meas];
    T_calc = [Temps.T1_calc; Temps.T2_calc; Temps.T3_calc; Temps.T4_calc; Temps.T5_calc];
    residuals = T_calc - T_meas;
    
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

function plotWakeWidthParity(T)
    if ~ismember('dry_wake_W_calc', T.Properties.VariableNames) || ~ismember('dry_wake_W', T.Properties.VariableNames)
        error('No wake width data found. Ensure T.dry_wake_W and T.dry_wake_W_calc are populated in table T.');
    end

    meas_all = [];
    calc_all = [];
    tc_labels = [];
    
    for i = 1:height(T)
        meas_i = T.dry_wake_W{i};
        calc_i = T.dry_wake_W_calc{i};
        
        if ~isempty(meas_i) && ~isempty(calc_i)
            meas_all = [meas_all, meas_i];
            calc_all = [calc_all, calc_i];
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
    
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
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

function plotWakeWidthProfileComparison(T, selectedRuns)
    figure('Name', 'Wake Width Profile Comparison', 'Color', 'w', 'Position', [200, 200, 780, 500]);
    hold on;
    x_tc = 1:5; 
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
        
        latexRunName = strrep(runName, '_', '\_');
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

function plotWakeWidthVsTempMeas(T, Temps)
    % Extract measured wake widths (N x 5) and measured temperatures (N x 5)
    W_meas_all = vertcat(T.dry_wake_W{:}); 
    T_meas_all = [Temps.T1_meas, Temps.T2_meas, Temps.T3_meas, Temps.T4_meas, Temps.T5_meas]; 
    
    % Unroll column-wise to maintain TC 1 to TC 5 grouping
    W_flat = W_meas_all(:);
    T_flat = T_meas_all(:);
    
    numRows = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Wake Width vs Measured Temp', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
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

function plotWakeWidthVsTempCalc(T, Temps)
    % Extract measured wake widths (N x 5) and calculated temperatures (N x 5)
    W_meas_all = vertcat(T.dry_wake_W{:}); 
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc]; 
    
    % Unroll column-wise to maintain TC 1 to TC 5 grouping
    W_flat = W_meas_all(:);
    T_flat = T_calc_all(:);
    
    numRows = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Wake Width vs Calculated Temp', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
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

function plotWakeWidthCalcVsTempCalc(T, Temps)
    % Extract calculated wake widths (N x 5) and calculated temperatures (N x 5)
    W_calc_all = vertcat(T.dry_wake_W_calc{:}); 
    T_calc_all = [Temps.T1_calc, Temps.T2_calc, Temps.T3_calc, Temps.T4_calc, Temps.T5_calc]; 
    
    % Unroll column-wise to maintain TC 1 to TC 5 grouping
    W_flat = W_calc_all(:);
    T_flat = T_calc_all(:);
    
    numRows = height(Temps);
    tc_labels = [ones(numRows, 1); 2*ones(numRows, 1); 3*ones(numRows, 1); ...
                 4*ones(numRows, 1); 5*ones(numRows, 1)];
    
    figure('Name', 'Calc Wake Width vs Calc Temp', 'Color', 'w', 'Position', [300, 300, 650, 550]);
    hold on;
    
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
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