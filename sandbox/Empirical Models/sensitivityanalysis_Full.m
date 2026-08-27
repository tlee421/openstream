% -- sensitivityAnalysis_Full.m --
%
% Complete script including:
% 1) Baseline Parity Plot (T_calc vs T_meas for TC1–TC5).
% 2) OAT Sensitivity Sweep generating dual Tornado Plots (RMSE and MBE).
% 3) Implicit iterative entrainment solver in computeAnnularFlow.

import Solvers.Obstruction.WallHT.ftoglass_ann_optim

%% ---------------- User-configurable settings -----------------------------
dp              = 0.1;       % fractional perturbation size (10%)
USE_COARSE_MESH = false;     % coarser wall-conduction mesh -> faster
COARSE_MN       = 40;        % M=N if USE_COARSE_MESH
ROW_SUBSET      = [];        % e.g. 1:10; [] = all rows
OUTDIR          = pwd;       

%% ---------------- Load data ----------------------------------------------
varNames = {'Param', 'T1', 'T2', 'T3', 'T4', 'T5', 'Tsat', 'G', 'x_in', ...
            'q_eff', 'Q_heater', 'delta', 'Re_f', 'SWR', 'mdot_f', 'tau_w', 'Re_v'};

T = readtable("C:\Users\tlee\OneDrive - UW-Madison\Documents\01_Research\07_Conferences\ASME SHTC\Data\WallHTModelComparison.xlsx", ...
    'VariableNamingRule', 'preserve');
T.Properties.VariableNames = varNames;

if ~isempty(ROW_SUBSET)
    T = T(ROW_SUBSET, :);
end

load("testfluidprops.mat")
geom.k_film  = fluid.KF;
geom.k_vapor = fluid.KG;
geom.rho_f   = fluid.RHOF;
mu_f         = fluid.MUF;
geom.mu_f    = mu_f;
geom.nu_f    = mu_f/geom.rho_f;
geom.PR_v    = fluid.PRANDTLG;
geom.h_fg    = fluid.HFG;
geom.T_inf   = fluid.TSAT;
geom.T_amb   = 22.5 + 273.15;

if isfield(fluid, 'RHOG'),  geom.rho_v   = fluid.RHOG;  else, geom.rho_v   = 14.8; end
if isfield(fluid, 'MUG'),   geom.mu_v    = fluid.MUG;   else, geom.mu_v    = 1.05e-5; end
if isfield(fluid, 'SIGMA'), geom.sigma_l = fluid.SIGMA; else, geom.sigma_l = 0.014; end

% Wall / channel geometry
geom.W      = 0.036;
Hchan       = 0.012;
per_        = 2*geom.W + 2*Hchan;
A_c         = geom.W*Hchan;
geom.A_c    = A_c;
geom.D_h    = 4*A_c/per_;
L           = 0.525;
L_obs       = 0.42;
geom.L      = L;
geom.L_obs  = L_obs;
geom.D_obs  = 0.0025;
geom.per    = per_;

geom.x_obs  = L - L_obs;
geom.x_targets = [0.04 0.05 0.06 0.07 0.08];   

if USE_COARSE_MESH
    geom.M = COARSE_MN; geom.N = COARSE_MN;
else
    geom.M = 125; geom.N = 125;
end
geom.W_HT   = geom.W/2;
geom.H_wall = 0.003;
geom.k_wall = 1.2;
geom.htc_amb= 10;

Tmeas = [T.T1 T.T2 T.T3 T.T4 T.T5];   

%% ---------------- Baseline parameter set ---------------------------------
par0.mdot_f   = T.mdot_f;
par0.Re_f     = T.Re_f;
par0.SWR      = T.SWR;
par0.tau_w    = T.tau_w;
par0.Re_v     = T.Re_v;

par0.G        = T.G;
par0.x_in     = T.x_in;
par0.Q_heater = T.Q_heater;

par0.q_eff    = T.q_eff;
par0.phi      = 0.6 * ones(height(T),1);   

par0.Wad_C    = 1.475;
par0.Wad_a    = 0.085;
par0.Wad_b    = 0.22;
par0.yplus    = 5;
par0.DB_C     = 0.023;
par0.DB_m     = 0.8;
par0.DB_n     = 0.4;
par0.wr_p1    = -5;
par0.wr_p2    = 0.0058;
par0.wr_p3    = 6.3325;
par0.wr_p4    = 1.6617;
par0.wr_p5    = 0.1719;

par0.scale_Wad     = 1;
par0.scale_dryW    = 1;
par0.scale_htc_wet = 1;
par0.scale_htc_dry = 1;

%% ---------------- Baseline model evaluation & Parity Plot -----------------
fprintf('Running baseline model (%d rows)...\n', height(T));
[Tcalc0, RMSE0, resid0, RMSE_row0, MBE0] = evalModel(par0, geom, Tmeas);
fprintf('Baseline RMSE = %.4f °C\n', RMSE0);
fprintf('Baseline MBE  = %+.4f °C\n\n', MBE0);

% Build table structure required for plotParity
TempsTable = table(Tmeas(:,1), Tmeas(:,2), Tmeas(:,3), Tmeas(:,4), Tmeas(:,5), ...
                   Tcalc0(:,1), Tcalc0(:,2), Tcalc0(:,3), Tcalc0(:,4), Tcalc0(:,5), ...
                   'VariableNames', {'T1_meas','T2_meas','T3_meas','T4_meas','T5_meas', ...
                                     'T1_calc','T2_calc','T3_calc','T4_calc','T5_calc'});

% Generate baseline Parity Plot
plotParity(TempsTable);

%% ---------------- Parameter perturbation list -----------------------------
paramList = {
    'x_in',         'Input/Experimental', 'x_in (inlet quality)',          'vector'
    'G',            'Input/Experimental', 'G (mass flux)',                  'vector'
    'q_eff',        'Input/Experimental', 'q_eff (effective heat flux)',    'vector'
    'phi',          'Closure Relation',   '\phi (mass division fraction)', 'vector'
    'scale_Wad',    'Closure Relation',   'Adiabatic wake width',          'scalar'
    'scale_dryW',   'Closure Relation',   'Dry wake width',                'scalar'
    'scale_htc_wet','Closure Relation',   'Free-stream (wet) HTC',         'scalar'
    'scale_htc_dry','Closure Relation',   'Wake (dry) HTC',                'scalar'
};
n_par = size(paramList, 1);

%% ---------------- OAT finite-difference sensitivity sweep -----------------
dRMSE_plus   = nan(n_par, 1);
dRMSE_minus  = nan(n_par, 1);
dMBE_plus    = nan(n_par, 1);
dMBE_minus   = nan(n_par, 1);

fprintf('Running OAT sensitivity sweep on %d parameters...\n', n_par);
for p = 1:n_par
    field = paramList{p,1};

    par_plus  = par0;
    par_minus = par0;
    par_plus.(field)  = par0.(field) * (1 + dp);
    par_minus.(field) = par0.(field) * (1 - dp);

    if strcmp(field, 'G') || strcmp(field, 'x_in')
        par_plus.is_hydro_perturbed  = true;
        par_minus.is_hydro_perturbed = true;
    end

    [~, RMSE_plus,  ~, ~, MBE_plus]  = evalModel(par_plus,  geom, Tmeas);
    [~, RMSE_minus, ~, ~, MBE_minus] = evalModel(par_minus, geom, Tmeas);

    % Absolute change in RMSE [°C]
    dRMSE_plus(p)  = RMSE_plus  - RMSE0;
    dRMSE_minus(p) = RMSE_minus - RMSE0;

    % Absolute change in MBE [°C]
    dMBE_plus(p)  = MBE_plus  - MBE0;
    dMBE_minus(p) = MBE_minus - MBE0;

    fprintf('  [%2d/%2d] %-14s | dRMSE (+10%%: %+.3f °C | -10%%: %+.3f °C) | dMBE (+10%%: %+.3f °C | -10%%: %+.3f °C)\n', ...
        p, n_par, field, dRMSE_plus(p), dRMSE_minus(p), dMBE_plus(p), dMBE_minus(p));
end

%% ---------------- Consolidated Results Table ------------------------------
Results = table(paramList(:,1), paramList(:,2), paramList(:,3), ...
    dRMSE_plus, dRMSE_minus, dMBE_plus, dMBE_minus, ...
    'VariableNames', {'Parameter','Category','Label','dRMSE_plus','dRMSE_minus','dMBE_plus','dMBE_minus'});

writetable(Results, fullfile(OUTDIR, 'sensitivityResults_Full.csv'));

catColors = containers.Map({'Input/Experimental','Closure Relation'}, ...
                            {[0.20 0.45 0.75], [0.85 0.35 0.10]});

%% =========================================================================
%% FIGURE 1: RMSE TORNADO PLOT
%% =========================================================================
[~, ordRMSE] = sort(max(abs(Results.dRMSE_plus), abs(Results.dRMSE_minus)), 'descend');
ResRMSE = Results(ordRMSE, :);

figure('Name','Sensitivity Tornado Plot (RMSE)','Color','w','Position',[100 100 800 500]);
hold on;

nR = height(ResRMSE);
yPos = nR:-1:1;
barWidth = 0.6;

for i = 1:nR
    c = catColors(ResRMSE.Category{i});
    valPlus  = ResRMSE.dRMSE_plus(i);
    valMinus = ResRMSE.dRMSE_minus(i);
    
    if abs(valPlus) >= abs(valMinus)
        barh(yPos(i), valPlus,  barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
        barh(yPos(i), valMinus, barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.35);
    else
        barh(yPos(i), valMinus, barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.35);
        barh(yPos(i), valPlus,  barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
    end
end

xline(0, 'k-', 'LineWidth', 1.2);
yticks(1:nR);
yticklabels(flip(ResRMSE.Label));
xlabel(sprintf('Change in Root Mean Square Error \\DeltaRMSE [°C] for a ±%.0f%% perturbation', dp*100));
title(sprintf('Sensitivity of Wall-Temperature RMSE (Baseline RMSE = %.2f °C)', RMSE0));
grid on; box on;

h1 = patch(NaN, NaN, catColors('Input/Experimental'), 'EdgeColor', 'k', 'FaceAlpha', 0.85);
h2 = patch(NaN, NaN, catColors('Closure Relation'),   'EdgeColor', 'k', 'FaceAlpha', 0.85);
lgd1 = legend([h1 h2], {'Input / Experimental', 'Closure Relation'}, 'Location', 'southeast');
title(lgd1, sprintf('Dark: +%.0f%% | Light: -%.0f%%', dp*100, dp*100), 'FontWeight', 'normal', 'FontSize', 9);
set(gca, 'FontSize', 10, 'YDir', 'normal');

%% =========================================================================
%% FIGURE 2: MBE TORNADO PLOT
%% =========================================================================
[~, ordMBE] = sort(max(abs(Results.dMBE_plus), abs(Results.dMBE_minus)), 'descend');
ResMBE = Results(ordMBE, :);

figure('Name','Sensitivity Tornado Plot (MBE)','Color','w','Position',[150 150 800 500]);
hold on;

for i = 1:nR
    c = catColors(ResMBE.Category{i});
    valPlus  = ResMBE.dMBE_plus(i);
    valMinus = ResMBE.dMBE_minus(i);
    
    if abs(valPlus) >= abs(valMinus)
        barh(yPos(i), valPlus,  barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
        barh(yPos(i), valMinus, barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.35);
    else
        barh(yPos(i), valMinus, barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.35);
        barh(yPos(i), valPlus,  barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
    end
end

xline(0, 'k-', 'LineWidth', 1.2);
yticks(1:nR);
yticklabels(flip(ResMBE.Label));
xlabel(sprintf('Change in Mean Bias Error \\DeltaMBE [°C] for a ±%.0f%% perturbation', dp*100));
title(sprintf('Sensitivity of Wall-Temperature MBE (Baseline MBE = %+.2f °C)', MBE0));
grid on; box on;

h1 = patch(NaN, NaN, catColors('Input/Experimental'), 'EdgeColor', 'k', 'FaceAlpha', 0.85);
h2 = patch(NaN, NaN, catColors('Closure Relation'),   'EdgeColor', 'k', 'FaceAlpha', 0.85);
lgd2 = legend([h1 h2], {'Input / Experimental', 'Closure Relation'}, 'Location', 'southeast');
title(lgd2, sprintf('Dark: +%.0f%% | Light: -%.0f%%', dp*100, dp*100), 'FontWeight', 'normal', 'FontSize', 9);
set(gca, 'FontSize', 10, 'YDir', 'normal');


%% ---------------- Local Functions -----------------------------------------

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
    
    % Axis limits
    lims = [min([T_meas; T_calc]) - 2, max([T_meas; T_calc]) + 2];
    x_band = linspace(lims(1), lims(2), 100);
    y_upper = 1.20 * x_band;
    y_lower = 0.80 * x_band;
   
    % --- 1. Draw the 1:1 Parity Line ---
    h_parity = plot(lims, lims, 'k--', 'LineWidth', 1.5);
    
    % --- 2. Plot the Actual Data Points ---
    colors = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.93 0.69 0.13; 0.49 0.18 0.56; 0.47 0.67 0.19];
    h_scatter = gscatter(T_meas, T_calc, TC_labels, colors, 'o', 6, 'on');
    
    % Fill each group's markers with its own color
    for i = 1:numel(h_scatter)
        h_scatter(i).MarkerFaceColor = h_scatter(i).Color;
    end
    
    % --- 3. Formatting & Labels ---
    xlabel('Measured Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('Calculated Temperature ($^{\circ}\mathrm{C}$)', 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Baseline Parity (RMSE = %.2f \\circC, R^2 = %.3f)', RMSE, R2), 'Interpreter', 'tex');
 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 11);
    
    legendLabels = {'1:1 Agreement ($y=x$)', ...
                    'TC 1', 'TC 2', 'TC 3', 'TC 4', 'TC 5'};
    legend([h_parity, h_scatter'], legendLabels, ...
           'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 9);
    
    grid on; 
    axis equal; 
    xlim(lims); 
    ylim(lims);
end

function [T_calc, RMSE, resid, RMSE_row, MBE] = evalModel(par, geom, Tmeas)
    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    if isfield(par, 'is_hydro_perturbed') && par.is_hydro_perturbed
        [mdot_f, Re_f, SWR, tau_w, Re_v, ~] = computeAnnularFlow(par.G, par.x_in, par.Q_heater, geom);
    else
        mdot_f = par.mdot_f;
        Re_f   = par.Re_f;
        SWR    = par.SWR;
        tau_w  = par.tau_w;
        Re_v   = par.Re_v;
    end

    n = numel(mdot_f);
    n_x = numel(geom.x_targets);

    % --- Closure Calculations ---
    W_adiabatic = par.Wad_C .* Re_f.^par.Wad_a .* SWR.^par.Wad_b * 2/1000;
    W_adiabatic = W_adiabatic .* par.scale_Wad;

    u_star     = sqrt(tau_w / geom.rho_f);
    base_thick = par.yplus ./ u_star * geom.nu_f;
    htc_wet    = geom.k_film ./ base_thick;
    htc_wet    = htc_wet .* par.scale_htc_wet;

    htc_dry = geom.k_vapor/geom.D_h * par.DB_C .* Re_v.^par.DB_m * geom.PR_v.^par.DB_n;
    htc_dry = htc_dry .* par.scale_htc_dry;

    mdot_wake  = mdot_f * (geom.D_obs/geom.per) .* par.phi;
    scalarPart = (par.q_eff .* W_adiabatic) ./ (mdot_wake .* geom.h_fg);
    Bo = scalarPart .* geom.x_targets;

    wake_ratio = max(0, par.wr_p1 + par.wr_p2 .* Bo + ...
                 (par.wr_p3 .* Bo.^par.wr_p4) ./ (par.wr_p5.^par.wr_p4 + Bo.^par.wr_p4));
    dry_wake_W = W_adiabatic .* wake_ratio;
    dry_wake_W = dry_wake_W .* par.scale_dryW;
    w_free     = geom.W - dry_wake_W;

    % --- Solve 2-D Wall Conduction ---
    T_calc = zeros(n, n_x);
    for i = 1:n
        for k = 1:n_x
            [~,~,~, Tc, ~,~,~,~,~,~] = ftoglass_ann_optim( ...
                geom.M, geom.N, geom.W_HT, geom.H_wall, geom.k_wall, geom.x_obs, ...
                par.q_eff(i), w_free(i,k), htc_dry(i), htc_wet(i), ...
                geom.htc_amb, geom.T_inf, geom.T_amb);
            T_calc(i,k) = Tc;
        end
    end

    resid    = T_calc - Tmeas;
    RMSE_row = sqrt(mean(resid.^2, 2));
    RMSE     = sqrt(mean(resid(:).^2));
    MBE      = mean(resid(:));
end

function [mdot_f, Re_f, SWR, tau_w, Re_v, delta] = computeAnnularFlow(G, x_in, Q_heater, geom)
    n = numel(G);
    mdot_f = zeros(n,1); Re_f = zeros(n,1); SWR = zeros(n,1);
    tau_w = zeros(n,1);  Re_v = zeros(n,1); delta = zeros(n,1);
    g = 9.81;

    for i = 1:n
        % 1. Total and phase mass flow rates
        m_dot_total = G(i) * geom.A_c;
        q_in_obs    = Q_heater(i); % Matching EES model definition
        x_local     = x_in(i) + q_in_obs / (m_dot_total * geom.h_fg);
        x_local     = max(min(x_local, 0.99), 0.01);

        m_dot_l = m_dot_total * (1 - x_local);
        m_dot_v = m_dot_total * x_local;

        u_v     = m_dot_v / (geom.A_c * geom.rho_v);
        Re_v(i) = geom.rho_v * u_v * geom.D_h / geom.mu_v;

        % 2. Implicit Iterative Entrainment Solver
        E_guess   = 0.5;
        tol       = 1e-6;
        maxIter   = 100;
        converged = false;

        for iter = 1:maxIter
            rho_c = (x_local + E_guess * (1 - x_local)) / ...
                    ((x_local / geom.rho_v) + (E_guess * (1 - x_local) / geom.rho_f));
            
            We_c = rho_c * (u_v^2) * geom.D_h / geom.sigma_l;
            E_calc = (1 + 279.6 * We_c^(-0.8395))^(-2.209);
            
            if abs(E_calc - E_guess) < tol
                converged = true;
                break;
            end
            E_guess = 0.5 * E_guess + 0.5 * E_calc;
        end

        if ~converged
            E_calc = E_guess;
        end

        m_dot_d  = m_dot_l * E_calc;
        mdot_f(i)= m_dot_l - m_dot_d;
        Re_f(i)  = mdot_f(i) / (geom.mu_f * geom.per);

        % 3. Film thickness and wall shear stress
        f_w = max(16/Re_f(i), 0.005);
        objFun = @(d) (0.5 * f_w * geom.rho_f * (mdot_f(i)/(geom.rho_f * geom.per * d))^2) - ...
                      (0.5 * geom.rho_v * (u_v - (mdot_f(i)/(geom.rho_f * geom.per * d)))^2 * ...
                       (0.005 * (1 + 300 * (d / geom.D_h))));
        
        try
            d_sol = fzero(objFun, [1e-6, geom.D_h/2]);
        catch
            d_grid = linspace(1e-6, geom.D_h/2, 500);
            [~, idx] = min(abs(arrayfun(objFun, d_grid)));
            d_sol = d_grid(idx);
        end

        delta(i) = d_sol;
        u_f      = mdot_f(i) / (geom.rho_f * geom.per * delta(i));
        f_i      = 0.005 * (1 + 300 * (delta(i) / geom.D_h));
        tau_int  = 0.5 * geom.rho_v * (u_v - u_f)^2 * f_i;
        
        tau_w(i) = tau_int;
        SWR(i)   = tau_int / (geom.rho_f * delta(i) * g);
    end
end