classdef WAKEWIDTH
    %WAKEWIDTH Enumeration of wake track width models for obstruction
    %solver
    %
    % This class defines the available :attr:`Inputs.Model.WAKEWIDTH`
    % models for calculating the width of wake track in the two track
    % obstruction solver.
    %
    % Models:
    %
    % - OBSWIDTH     — Assumes width of wake is equal to projected perimeter of the obstruction
    % - LARGEOBS     — Uses closure relation developed for large obstruction wake width within the MFVAL facility

    enumeration
        OBSWIDTH          % Wake width is equal to obs projected perimeter
        LARGEOBS          % MFVAL large obstruction wake width closure model
    end
end