classdef SOLVER
%SOLVER Name of available solver
%   Defines the solver name for selecting which solver type to use
%
    enumeration
        MIXTURE
        TWOFLUID
        THREEFIELD
        FOURFIELD
        OBSTRUCTION
    end

    methods

        function solverPath = solverPath(enum)
            import InputEnums.SOLVER

            switch (enum)
                case SOLVER.MIXTURE     
                    solverPath = "Solvers.Mixture.MixtureSolver";
                case SOLVER.TWOFLUID
                    solverPath = "Solvers.TwoFluid.TwoFluidSolver";
                case SOLVER.THREEFIELD
                    solverPath = "Solvers.ThreeField.ThreeFieldSolver";
                case SOLVER.FOURFIELD
                    solverPath = "Solvers.FourField.FourFieldSolver";
                case SOLVER.OBSTRUCTION
                    solverPath = "Solvers.Obstruction.ObstructionSolver";
                otherwise
                    error("SOLVER:invalidSolver", "Solver type not recognized");
            end
        end

    end

end
