classdef WALLSHAPE
    %WAKEWIDTH Enumeration of wall shapes for conduction modeling
    %solver
    %
    % This class defines the available :attr:`Inputs.Model.WALLSHAPE`
    % models for wall heat transfer and conduction
    % Currently implemented for solving lateral conduction within the obstruction solver
    %
    % Models:
    %
    % - PLANAR          — Assumes wall has planar geometry for conduction models
    % - CYLINDRICAL     — Not yet available, to be implemented in the future 

    enumeration
        PLANAR            % Wall has planar geometry
    end
end