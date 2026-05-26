classdef THINFILMFRIC
    %THINFILMFRIC Enumeration of thin film friction models
    %
    % This class defines the available :attr:`Inputs.Model.THINFILMFRIC`
    % models for calculating wall friction under thin film conditions, used
    % by the three-field and four-field solvers.
    %
    % Models:
    %
    % - LAMINAR   — Laminar flow model (:math:`16/\mathrm{Re}_{film}`)
    % - TURBULENT — Turbulent flow model

    enumeration
        LAMINAR              % Laminar thin film friction model
        TURBULENT            % Turbulent thin film friction model
        TRANSITION           % User defined transition from laminar to turbulent thin film friction models
        TRACE                % Wall friction model for annular flows utilized in TRACE 
    end
end