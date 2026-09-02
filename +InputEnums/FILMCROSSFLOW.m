classdef FILMCROSSFLOW
    %MOMENTFILM Enumeration of models for film flow between adjacent walls
    %
    % This class defines the available :attr:`Inputs.Model.FILMCROSSFLOW` models
    % for solving the rate of film mass transfer between walls in the three
    % field solver. 
    %
    % Models:
    %
    % - NONE          — No film crossflow between walls is considered
    % - FTPOTENTIAL   — Rate of film mass transfer is based on difference in film thickness between adjacent walls 
   

    enumeration
        NONE            % no film crossflow
        FTPOTENTIAL     % Rate of film crossflow is based on differences in film thickness
    end
end