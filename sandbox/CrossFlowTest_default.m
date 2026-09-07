clearvars
import Inputs.*
import Solvers.*

warning off
% Construct the input set
inputSet = InputSet( ...
    modelFilePath          = 'inputs/models.inp', ...
    modelID                = 'TUTORIAL1', ...
    optionsFilePath        = 'inputs/options.inp', ...
    optionsID              = 'DEFAULT', ...
    geometryFilePath       = 'inputs/geom.inp', ...
    geometryID             = 'CROSSFLOWTEST', ...
    bcFilePath             = 'inputs/crossflowtest.inp', ...
    sessionParentDir       = fullfile(pwd,'outputs'), ...
    overwriteSessionFiles = true, ...
    LOGMODE                = 'BOTH');

warning on

import Solvers.Mixture.*

% Construct the mixture solver object and solve the flow field.
mixSolver = MixtureSolver(inputSet);
mixSolver.solve();

%% Three-field

import Solvers.ThreeField.*

% Construct the three-field solver object and solve the flow fields.
tfSolver = ThreeFieldSolver(inputSet, mixSolver);
tfSolver.solve();
%%
% Plot field mass flow rates and velocities.
tfSolver.plotz('display', "U");
