classdef ENTRAINMENT
    %ENTRAINMENT Enumeration of droplet entrainment models
    %
    % This class defines the available :attr:`Inputs.Model.ENTRAINMENT` models
    % for calculating the droplet entrainment mass flux, used by the three-field
    % and four-field solvers.
    %
    % NOTE: Entrainment and deposition correlations are coupled.
    % It is recommended to use the same model for both entrainment
    % and deposition mass flux to ensure consistency.
    %
    % Models:
    %
    % - NONE          — No entrainment
    % - GOVAN         — Hewitt and Govan correlation (:cite:t:`hewitt1990phenomenological`)
    % - OKAWA2003     — Okawa et al. correlation (:cite:t:`okawa2003`)
    % - OKAWA2004     — Okawa et al. correlation (:cite:t:`OKAWA2004`)
    % - OKAWA2004MOD  — Modified Okawa (2004) model from (:cite:t:`ADAMSSON20112843`)
    % - OKAWAGEN      — Generic Okawa-based model using user-defined :attr:`Inputs.Model.OKAWACOEFS`

    enumeration
        NONE                 % No entrainment
        GOVAN                % Hewitt and Govan (1990)
        OKAWA2003            % Okawa et al. (2003)
        OKAWA2004            % Okawa et al. (2004)
        OKAWA2004MOD         % Modified Okawa et al. (2004) from Adamsson and Le Corre (2011)
        OKAWAMFVAL           % Modified Okawa et al. (2004) with Rodarte (2015) data from MFVAL test facility
        OKAWACT              % Modified Okawa et al. (2004) based on fit from Ciancolini and Thome Correlation (2012)
        OKAWARD              % Modified Okawa et al. (2004) based on fit for refrigerants from Rodarte (2015)
        OKAWAGEN             % Generic Okawa model
    end
end