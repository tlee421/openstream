clear vars

%Import Inputs and Solvers classes
import Inputs.*
import Solvers.*

%turn of warnings during initial setup
warning off
%Create input set 
%Inputs
 inputSet = InputSet( ...
             modelFilePath         = './inputs/models.inp',  modelID    = 'MFVALS', ...
             optionsFilePath       = './inputs/options.inp', optionsID  = 'DEFAULT', ...
<<<<<<< HEAD
             geometryFilePath      = './inputs/geom.inp',    geometryID = 'MFVAL_dev', ...
             bcFilePath            = './inputs/bc_Regime3_c.inp', ...
=======
             geometryFilePath      = './inputs/geom.inp',    geometryID = 'MFVAL_pre_obs', ...
             bcFilePath            = './inputs/bc_mfval_pre_obs.inp', ...
>>>>>>> origin/lee
             sessionParentDir      = fullfile(pwd,'outputs'), ...
             overwriteSessionFiles = true, ...
             LOGMODE               = 'BOTH');
%Turn warnings back on 
warning on

%import mixture solver class
import Solvers.Mixture.*
% Create the mixture solver
mixSolver = MixtureSolver(inputSet);
%initialize the mixture solver
mixSolver.initializeSolver();
%solve the mixture solver
mixSolver.solve();
%plot the mixture solver mass flow rate results 
mixSolver.plotz('display','W');

%import three-field solver class
import Solvers.ThreeField.*
% Create the three-field solver
tfSolver = ThreeFieldSolver(inputSet, mixSolver);                                                        
%Solve the three-field solver
tfSolver.solve();
%Plot the three-field solver mass flow rate results
<<<<<<< HEAD
tfSolver.plotz('display','W', 'oafLine', false);
=======
tfSolver.plotz('display','THICK');
tfSolver.plotz('display','U');
tfSolver.plotz('display','FME');
tfSolver.plotz('display','FME', 'solveMode','NULL');
>>>>>>> origin/lee

