classdef SolverMode < uint16
%SOLVERMODE Enumeration class of possible solver modes
%
%   Solvers can be set up for a new solution, or a continuation of an
%   existing solution, or to solve a subset of the domain.
    enumeration
        NEW                    (0)
        CONTINUE               (1)
        SUBSET                 (2)
    end
end
