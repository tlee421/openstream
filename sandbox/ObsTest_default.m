%% Obstructed Solver

warning off
%Create inputSet object
inputSet= Inputs.InputSet( ...
         modelFilePath         = './inputs/models.inp',  modelID    = 'MFVALS', ...
         optionsFilePath       = './inputs/options.inp', optionsID  = 'OBS', ...
         geometryFilePath      = './inputs/geom.inp',    geometryID = 'MFVAL_obs', ...
         bcFilePath            = './inputs/bc_mfval_obs.inp', ...
         obsFilePath = './inputs/obstruction.inp', obsIDs = {'LARGEOBS'}, ...
         sessionParentDir      = fullfile(pwd,'outputs'), ...
         overwriteSessionFiles = true, ...
         LOGMODE               = 'BOTH');
warning on 

%Initialize obstruction solver
obsSolver = Solvers.Obstruction.ObstructionSolver(inputSet, 'THREEFIELD');

%Solve obstruction solver to compute hydrodynamics
obsSolver.solve();

%Solve wall heat transfer
obsSolver.solveWallHT();

%%
%Plot the three field model results for the three tracks merged together
obsSolver.plotz('display','THICK')
%%
%Plot the wall temperature distribution
obsSolver.plotWallTemperature();
