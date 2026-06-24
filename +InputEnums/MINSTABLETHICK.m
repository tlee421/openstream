classdef MINSTABLETHICK
    %WAKEWIDTH Enumeration of minimum stable film thickness models, used for detecting dryou 
    %
    %
    % This class defines the available :attr:`Inputs.Model.MINSTABLETHICK`
    % models for film dryout detection.
    %
    %
    % Models:
    %
    % - CHUN            — Chun et. al (2002) "Development of the critical film thickness correlation  for an advanced annular film mechanistic dryout model applicable to MARS code"
     

    enumeration
        CHUN            % Chun et. al (2002) 
    end
end