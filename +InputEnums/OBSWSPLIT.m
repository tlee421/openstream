classdef OBSWSPLIT
    %OBSWSPLIT Enumeration of film mass division at obstruction models
    %
    % This class defines the available :attr:`Inputs.Model.OBSWSPLIT`
    % models for calculating the ratio of mass that enters the wake in the near obstruction field
    %
    % Models:
    %
    % - SCALAR   — Ratio of film mass incipient upon obstruction that makes its way into the wake. 
    %              Valid range of this closure is between 0-1, with 0 indicating no mass travels around obs.
    %              A value of 1 represents all incipient film mass travling into the wake, as if the obs isnt there at all.
    %
    % - BASEFILM — Assumes all basefilm travels into the wake region. Basefilm thickness is computed based on BASEYPLUS closure relation value.
    %              A nusselt velocity profile is assumed to compute the mass flow rate of basefilm at a thickness specified by BASEYPLUS.

    enumeration
        SCALAR                 % Film mass division is a user defined scalar
        BASEFILM               % Basefilm travels into wake model
    end
end