classdef DRYWAKEWIDTH
    %WAKEWIDTH Enumeration of dry wake width models for obstruction 
    %solver
    %
    % This class defines the available :attr:`Inputs.Model.DRYWAKEWIDTH`
    % models for calculating the width of the dry wake within the 
    % obstruction solver. The width of the dry wake is utilized in solveWallHT()
    % by determining where the transition from dry to wet HTC lies
    %
    % Models:
    %
    % - EMPIRICAL     — Emperically determined wake width as a function of wake boiling number. Data was gathered in MFVAL large obstruction  
    % - TWOTRACK      — Assumes dry wake width is equal to width of hydrodynamic wake (wake track wall)
    %
    
    

    enumeration
        EMPIRICAL          % Empirical wake width model
        TWOTRACK           % Two-track wake width model
        
    end
end