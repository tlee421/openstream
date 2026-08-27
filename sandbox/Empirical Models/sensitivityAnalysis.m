% --sensitivityAnalysis.m--%
%
% Local (one-at-a-time, finite-difference) sensitivity analysis of the
% obstructed annular flow wall HT model against the WallHTModelComparison
% dataset.
%
% This is the "physical mechanism" variant of sensitivityAnalysis.m. The
% Input/Experimental category is unchanged (mdot_f, Re_f, SWR, tau_w,
% Re_v, q_eff, perturbed as a uniform +/-dp bias). The Closure Relation
% category is different: instead of perturbing the internal fit
% coefficients/exponents of each correlation, this script perturbs the
% physical QUANTITY that correlation produces, directly, as a
% multiplicative scale factor applied after the baseline closure is
% evaluated. E.g. instead of asking "what if the Dittus-Boelter leading
% coefficient is off by 5%", it asks "what if the wake HTC itself is off
% by 5%" -- which is agnostic to which part of the correlation's
% functional form is wrong and maps directly onto "how much can I trust
% this physical mechanism", which is the more useful question when you
% don't have a specific reason to suspect one exponent over another.
%
% WHAT IS PERTURBED
% ------------------
% Category "Input/Experimental"  (per-row columns, perturbed as a
% uniform +/-dp fractional bias across all rows simultaneously, to mimic
% a systematic measurement/upstream-closure error):
%   mdot_f, Re_f, SWR, tau_w, Re_v, q_eff
%
%   NOTE: G, x_in, and delta are NOT perturbed directly here. In this
%   script they only enter through Re_f, SWR, mdot_f, tau_w (already
%   computed upstream by your entrainment/hydrodynamic closure and
%   stored in the data table). Perturbing those proxies is the
%   traceable stand-in for G/x_in/delta uncertainty. If you want the
%   literal dG, dx_in, ddelta -> dRMSE mapping, that requires calling
%   your entrainment closure function inside evalModel() so the
%   perturbation actually propagates through Re_f/SWR/mdot_f/tau_w -
%   flag it and I'll wire that function in.
%
% Category "Closure Relation" (physical quantities, perturbed +/-dp
% fractionally as a direct multiplicative scale on the computed value,
% applied to every row):
%   W_adiabatic     - adiabatic (hydrodynamic) wake width
%   dry_wake_W      - dry wake width (post Bo/wake-ratio calc, scaled
%                      independently of W_adiabatic so you can isolate
%                      "the wake grows/shrinks vs Bo faster/slower than
%                      modeled" from "the base adiabatic width is wrong")
%   phi             - division of film mass between wake / outer wake
%                      (already a direct physical fraction, unchanged
%                      from the coefficient-level script)
%   htc_wet         - free-stream (wetted film) HTC
%   htc_dry         - wake (dry, vapor-convection) HTC
%
% OUTPUT
% ------
%   - Console table of elasticities (%RMSE change per %parameter change)
%   - Tornado plot ranking all parameters by |elasticity|
%   - Row x parameter elasticity heatmap (is sensitivity flow-regime
%     dependent?) for the top N most sensitive parameters
%   - sensitivityResults.csv / sensitivityResults.mat written to disk
%
% RUNTIME NOTE
% ------------
% Each parameter perturbation requires re-solving the 2-D wall
% conduction problem (ftoglass_ann_optim) at 5 target locations for
% every row, twice (+dp/-dp). With N_rows rows and N_par parameters
% that is N_rows * N_par * 2 * 5 solver calls. For N_rows ~ 30 and
% N_par ~ 15 that's ~4500 calls. Options if this is too slow:
%   1) Set USE_COARSE_MESH = true below to solve the sensitivity study
%      on a coarser (M,N) mesh than your production runs (ranking of
%      sensitivities is generally insensitive to mesh refinement).
%   2) Set ROW_SUBSET below to a smaller set of representative rows.
%   3) Parallelize the row loop inside evalModel with parfor (commented
%      stub included).

import Solvers.Obstruction.WallHT.ftoglass_ann_optim

%% ---------------- User-configurable sensitivity settings -----------------
dp              = 0.1;     % fractional perturbation size (5%) for all params
USE_COARSE_MESH = false;     % coarser wall-conduction mesh -> much faster
COARSE_MN       = 40;       % M=N=this if USE_COARSE_MESH
ROW_SUBSET      = [];       % e.g. 1:10 to use only first 10 rows; [] = all rows
N_TOP_HEATMAP   = 11;       % number of top-ranked params to show in heatmap (all of them, since the list is short now)
OUTDIR          = pwd;      % where to write sensitivityResults.csv/.mat

%% ---------------- Load data (mirrors empiricalModelTest.m) ---------------
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
geom.nu_f    = mu_f/geom.rho_f;
geom.PR_v    = fluid.PRANDTLG;
geom.h_fg    = fluid.HFG;
geom.T_inf   = fluid.TSAT;
geom.T_amb   = 22.5 + 273.15;

% Wall / channel geometry
geom.W      = 0.036;
Hchan       = 0.012;
per_        = 2*geom.W + 2*Hchan;
A_c         = geom.W*Hchan;
geom.D_h    = 4*A_c/per_;
L           = 0.525;
L_obs       = 0.42;
geom.D_obs  = 0.0025;
geom.per    = per_;

geom.x_obs  = L - L_obs;
geom.x_targets = [0.04 0.05 0.06 0.07 0.08];   % TC locations behind obstruction

% Wall HT solver params
if USE_COARSE_MESH
    geom.M = COARSE_MN; geom.N = COARSE_MN;
else
    geom.M = 125; geom.N = 125;
end
geom.W_HT   = geom.W/2;
geom.H_wall = 0.003;
geom.k_wall = 1.2;
geom.htc_amb= 10;

Tmeas = [T.T1 T.T2 T.T3 T.T4 T.T5];   % N_rows x 5, measured TC temps

%% ---------------- Baseline parameter set (fixed phi = 0.6) ---------------
par0.mdot_f = T.mdot_f;
par0.Re_f   = T.Re_f;
par0.SWR    = T.SWR;
par0.tau_w  = T.tau_w;
par0.Re_v   = T.Re_v;
par0.q_eff  = T.q_eff;

par0.phi    = 0.6 * ones(height(T),1);   % fixed phi baseline

% Fixed baseline closure-form coefficients (NOT perturbed in this script
% -- only the resulting physical quantities are perturbed, via the
% scale_* factors below).
par0.Wad_C  = 1.475;   % W_adiabatic = Wad_C * Re_f^Wad_a * SWR^Wad_b * 2/1000
par0.Wad_a  = 0.085;
par0.Wad_b  = 0.22;

par0.yplus  = 5;        % base-film-thickness y+ assumption (-> htc_wet)

par0.DB_C   = 0.023;    % htc_dry = k_vapor/D_h * DB_C * Re_v^DB_m * PR_v^DB_n
par0.DB_m   = 0.8;
par0.DB_n   = 0.4;

par0.wr_p1  = -5;        % wake_ratio = max(0, wr_p1 + wr_p2*Bo + ...
par0.wr_p2  = 0.0058;    %              (wr_p3*Bo^wr_p4)/(wr_p5^wr_p4 + Bo^wr_p4))
par0.wr_p3  = 6.3325;
par0.wr_p4  = 1.6617;
par0.wr_p5  = 0.1719;

% Direct multiplicative scale factors on the physical quantities
% themselves -- this is what gets perturbed for the "Closure Relation"
% category in this version of the script.
par0.scale_Wad     = 1;   % scales W_adiabatic directly
par0.scale_dryW    = 1;   % scales dry_wake_W directly (post Bo calc)
par0.scale_htc_wet = 1;   % scales the free-stream (wetted) HTC directly
par0.scale_htc_dry = 1;   % scales the wake (dry) HTC directly

%% ---------------- Baseline model evaluation ------------------------------
fprintf('Running baseline model (%d rows)...\n', height(T));
[Tcalc0, RMSE0, resid0, RMSE_row0] = evalModel(par0, geom, Tmeas);
fprintf('Baseline RMSE = %.4f C\n\n', RMSE0);

%% ---------------- Parameter perturbation list -----------------------------
% field, category, display label, perturbation kind ('vector' scales a
% per-row column, 'scalar' scales a single closure coefficient)
paramList = {
    'mdot_f',       'Input/Experimental', 'mdot_f (film mass flow)',            'vector'
    'Re_f',         'Input/Experimental', 'Re_f (film Reynolds number)',        'vector'
    'SWR',          'Input/Experimental', 'SWR (shear-to-weight ratio)',        'vector'
    'tau_w',        'Input/Experimental', 'tau_w (wall shear stress)',          'vector'
    'Re_v',         'Input/Experimental', 'Re_v (vapor Reynolds number)',       'vector'
    'q_eff',        'Input/Experimental', 'q_eff (effective heat flux)',        'vector'
    'phi',          'Closure Relation',   'phi (mass division fraction)',       'vector'
    'scale_Wad',     'Closure Relation',   'Adiabatic wake width',              'scalar'
    'scale_dryW',    'Closure Relation',   'Dry wake width',                    'scalar'
    'scale_htc_wet', 'Closure Relation',   'Free-stream (wet) HTC',             'scalar'
    'scale_htc_dry', 'Closure Relation',   'Wake (dry) HTC',                    'scalar'
};
n_par = size(paramList, 1);
N_rows = height(T);

%% ---------------- OAT finite-difference sensitivity sweep -----------------
Elasticity   = nan(n_par, 1);
dRMSE_plus   = nan(n_par, 1);
dRMSE_minus  = nan(n_par, 1);
RowElast     = nan(N_rows, n_par);  % per-row elasticity (using per-row RMSE)

fprintf('Running OAT sensitivity sweep on %d parameters...\n', n_par);
for p = 1:n_par
    field = paramList{p,1};
    kind  = paramList{p,4};

    par_plus  = par0;
    par_minus = par0;
    par_plus.(field)  = par0.(field) * (1 + dp);
    par_minus.(field) = par0.(field) * (1 - dp);

    [~, RMSE_plus,  ~, RMSE_row_plus]  = evalModel(par_plus,  geom, Tmeas);
    [~, RMSE_minus, ~, RMSE_row_minus] = evalModel(par_minus, geom, Tmeas);

    relChange_plus  = (RMSE_plus  - RMSE0) / RMSE0;
    relChange_minus = (RMSE_minus - RMSE0) / RMSE0;

    dRMSE_plus(p)  = relChange_plus  * 100;   % percent
    dRMSE_minus(p) = relChange_minus * 100;   % percent
    Elasticity(p)  = (relChange_plus - relChange_minus) / (2*dp);

    % per-row elasticity (does sensitivity depend on flow condition?)
    % (eps guard avoids Inf/NaN for the rare row with RMSE_row0 ~ 0)
    rowRelChange_plus  = (RMSE_row_plus  - RMSE_row0) ./ max(RMSE_row0, 1e-6);
    rowRelChange_minus = (RMSE_row_minus - RMSE_row0) ./ max(RMSE_row0, 1e-6);
    RowElast(:,p) = (rowRelChange_plus - rowRelChange_minus) / (2*dp);

    fprintf('  [%2d/%2d] %-14s (%-19s) elasticity = %+7.3f  (RMSE %+.2f%% / %+.2f%%)\n', ...
        p, n_par, field, kind, Elasticity(p), dRMSE_plus(p), dRMSE_minus(p));
end

%% ---------------- Results table -------------------------------------------
Results = table(paramList(:,1), paramList(:,2), paramList(:,3), ...
    Elasticity, dRMSE_plus, dRMSE_minus, ...
    'VariableNames', {'Parameter','Category','Label','Elasticity','pctRMSE_plus','pctRMSE_minus'});
[~, ord] = sort(abs(Results.Elasticity), 'descend', 'MissingPlacement','last');
Results = Results(ord, :);

disp(' ');
disp('=== Sensitivity ranking (by |elasticity|) ===');
disp(Results(:, {'Parameter','Category','Elasticity','pctRMSE_plus','pctRMSE_minus'}));

writetable(Results, fullfile(OUTDIR, 'sensitivityResults_physical.csv'));
save(fullfile(OUTDIR, 'sensitivityResults_physical.mat'), 'Results', 'RowElast', 'paramList', ...
     'RMSE0', 'RMSE_row0', 'dp', 'par0');

%% ---------------- Tornado plot (Original Alpha Shading) -------------------
figure('Name','Sensitivity Tornado Plot (Physical Quantities)','Color','w','Position',[100 100 800 550]);
hold on;

% Category colors
catColors = containers.Map({'Input/Experimental','Closure Relation'}, ...
                            {[0.20 0.45 0.75], [0.85 0.35 0.10]});

nR = height(Results);
yPos = nR:-1:1;
barWidth = 0.6;

for i = 1:nR
    c = catColors(Results.Category{i});
    
    valPlus  = Results.pctRMSE_plus(i);
    valMinus = Results.pctRMSE_minus(i);
    
    % Draw the bar with the larger magnitude FIRST so the smaller bar 
    % sits on top rather than being completely hidden underneath.
    if abs(valPlus) >= abs(valMinus)
        barh(yPos(i), valPlus,  barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
        barh(yPos(i), valMinus, barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.35);
    else
        barh(yPos(i), valMinus, barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.35);
        barh(yPos(i), valPlus,  barWidth, 'FaceColor', c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
    end
end

% Reference line at zero
xline(0, 'k-', 'LineWidth', 1.2);

% Formatting axes
yticks(1:nR);
yticklabels(flip(Results.Label));
xlabel(sprintf('Change in RMSE [%%] for a %.0f%% parameter perturbation', dp*100));
title('Sensitivity of Wall-Temperature RMSE to Inputs and Physical Closure Quantities');
grid on; 
box on;

% Legend setup
h1 = patch(NaN, NaN, catColors('Input/Experimental'), 'EdgeColor', 'k', 'FaceAlpha', 0.85);
h2 = patch(NaN, NaN, catColors('Closure Relation'),   'EdgeColor', 'k', 'FaceAlpha', 0.85);
lgd = legend([h1 h2], {'Input / Experimental', 'Closure Relation'}, 'Location', 'southeast');

% Legend subtitle for parameter directions
title(lgd, 'Dark: +10% | Light: -10%', 'FontWeight', 'normal', 'FontSize', 9);

set(gca, 'FontSize', 10, 'YDir', 'normal');
%% ---------------- Row x parameter elasticity heatmap -----------------------
topN = min(N_TOP_HEATMAP, n_par);
topFields = Results.Parameter(1:topN);
[~, colIdx] = ismember(topFields, paramList(:,1));

figure('Name','Row-wise Sensitivity Heatmap (Physical Quantities)','Color','w','Position',[150 150 850 600]);
imagesc(RowElast(:, colIdx));
try
    colormap(gca, redbluecmap());   % diverging colormap, see local fn below
catch
    colormap(gca, parula);
end
colorbar;
xticks(1:topN);
xticklabels(paramList(colIdx,3));
xtickangle(35);
yticks(1:N_rows);
yticklabels(T.Param);
set(gca,'FontSize',8,'TickLabelInterpreter','none');
xlabel('Parameter');
ylabel('Test Run');
title(sprintf('Row-wise Elasticity (top %d most sensitive parameters)', topN));

fprintf('\nDone. Results written to:\n  %s\n  %s\n', ...
    fullfile(OUTDIR,'sensitivityResults_physical.csv'), fullfile(OUTDIR,'sensitivityResults_physical.mat'));

%% =========================== Local functions ===============================

function [T_calc, RMSE, resid, RMSE_row] = evalModel(par, geom, Tmeas)
    % Parameterized re-implementation of the empiricalModelTest.m forward
    % chain (fixed-phi version). Any field of `par` can be perturbed by
    % the caller to run a sensitivity study without touching the rest of
    % the pipeline.
    %
    % par fields (see baseline block above for definitions):
    %   mdot_f, Re_f, SWR, tau_w, Re_v, q_eff   (N_rows x 1)
    %   phi                                     (N_rows x 1, or scalar)
    %   Wad_C, Wad_a, Wad_b, yplus, DB_C, DB_m, DB_n, wr_p1..wr_p5
    %                                           (scalars, fixed baseline
    %                                            closure-form coefficients
    %                                            -- NOT perturbed here)
    %   scale_Wad, scale_dryW, scale_htc_wet, scale_htc_dry
    %                                           (scalars, multiplicative
    %                                            perturbation applied
    %                                            directly to the physical
    %                                            quantity -- THIS is what
    %                                            gets perturbed here)
    %
    % geom fields: rho_f, nu_f, k_film, k_vapor, PR_v, D_h, W, D_obs,
    %              per, h_fg, x_targets, M, N, W_HT, H_wall, k_wall,
    %              x_obs, htc_amb, T_inf, T_amb

    import Solvers.Obstruction.WallHT.ftoglass_ann_optim

    n = numel(par.mdot_f);
    n_x = numel(geom.x_targets);

    % --- Adiabatic wake width closure, then direct scale perturbation ---
    W_adiabatic = par.Wad_C .* par.Re_f.^par.Wad_a .* par.SWR.^par.Wad_b * 2/1000;   % [m], N x 1
    W_adiabatic = W_adiabatic .* par.scale_Wad;

    % --- Wet-wall (free-stream) HTC, then direct scale perturbation ---
    u_star     = sqrt(par.tau_w / geom.rho_f);                 % N x 1
    base_thick = par.yplus ./ u_star * geom.nu_f;              % N x 1
    htc_wet    = geom.k_film ./ base_thick;                    % N x 1
    htc_wet    = htc_wet .* par.scale_htc_wet;

    % --- Dry-wall (wake) HTC (Dittus-Boelter), then direct scale perturbation ---
    htc_dry = geom.k_vapor/geom.D_h * par.DB_C .* par.Re_v.^par.DB_m * geom.PR_v.^par.DB_n;  % N x 1
    htc_dry = htc_dry .* par.scale_htc_dry;

    % --- Mass division at the obstruction / dry wake width vs Bo ---
    mdot_wake  = par.mdot_f * (geom.D_obs/geom.per) .* par.phi;                     % N x 1
    scalarPart = (par.q_eff .* W_adiabatic) ./ (mdot_wake .* geom.h_fg);            % N x 1
    Bo = scalarPart .* geom.x_targets;                                              % N x n_x

    wake_ratio = max(0, par.wr_p1 + par.wr_p2 .* Bo + ...
                 (par.wr_p3 .* Bo.^par.wr_p4) ./ (par.wr_p5.^par.wr_p4 + Bo.^par.wr_p4));
    dry_wake_W = W_adiabatic .* wake_ratio;    % N x n_x
    dry_wake_W = dry_wake_W .* par.scale_dryW; % direct scale perturbation
    w_free     = geom.W - dry_wake_W;          % N x n_x

    % --- Solve 2-D wall conduction at each target location ---
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
    RMSE_row = sqrt(mean(resid.^2, 2));         % N x 1, per-run RMSE (across 5 TCs)
    RMSE     = sqrt(mean(resid(:).^2));         % scalar, RMSE across all runs & TCs
end

function cmap = redbluecmap()
    % Small diverging colormap (blue-white-red) with no toolbox dependency,
    % used for the row-wise elasticity heatmap so negative/positive
    % sensitivities are visually distinct. Falls back to parula above if
    % this errors for any reason.
    n = 128;
    top = [linspace(1,0.7,n)', linspace(0,0,n)', linspace(0,0,n)'];
    bot = [linspace(0,0,n)', linspace(0,0,n)', linspace(0.7,1,n)'];
    cmap = [flipud(bot); top];
end
