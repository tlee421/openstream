%close all;
% clearvars

import Inputs.*
import Solvers.*
import Solvers.Mixture.*

%Inputs
 inputSet = InputSet( ...
             modelFilePath         = './inputs/models.inp',  modelID    = 'MFVALS', ...
             optionsFilePath       = './inputs/options.inp', optionsID  = 'DEFAULT', ...
             geometryFilePath      = './inputs/geom.inp',    geometryID = 'MFVAL_dev', ...
             bcFilePath            = './inputs/bc_Regime3_c.inp', ...
             sessionParentDir      = fullfile(pwd,'outputs'), ...
             overwriteSessionFiles = true, ...
             LOGMODE               = 'BOTH');
        
%inputSet = InputSet( ...
%            modelFilePath         = './inputs/models.inp',  modelID    = 'BARCS', ...
%            optionsFilePath       = './inputs/options.inp', optionsID  = 'STEADY', ...
%            geometryFilePath      = './inputs/geom.inp',    geometryID = 'BARC', ...
%            bcFilePath            = './inputs/bc_barc.inp', ...
%            sessionParentDir      = fullfile(pwd,'outputs'), ...
%            overwriteSessionFiles = true, ...
%            LOGMODE               = 'BOTH');        % Use LOGTOFILEONLY for (almost) no print ouptut
%         
% Create the mixture solver
mixSolver = MixtureSolver(inputSet);

% mixSolver is initialized at construction. Here, it is explicitly initialized for clarity.
mixSolver.initializeSolver(); 

% Solve (does not accept any argument)
mixSolver.solve();

% Axial and temporal plots
mixSolver.plotz(1);
mixSolver.plotz(mixSolver.NTIME);
%mixSolver.plott(mixSolver.NZ,'solveMode','STEADY')
%mixSolver.plott(mixSolver.NZ)
%mixSolver.plott([1 floor(mixSolver.NZ./8.*(2:8))])

% Save results
%mixSolver.saveResults(saveFormat="MAT");

%return

%% Three-field

import Solvers.ThreeField.*

% Create the three-field solver
tfSolver = ThreeFieldSolver(inputSet, mixSolver);                                                        % Plot initialization state (dep and ent set to 0 but calculated value is plotted instead)

% Plot the axial distribution of mass flow rates
tfSolver.plotz('display','W');

% Plot the axial distribution of velocities
tfSolver.plotz('display','U');