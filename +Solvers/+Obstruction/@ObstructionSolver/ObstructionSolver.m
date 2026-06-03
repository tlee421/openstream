classdef ObstructionSolver < Solvers.AbstractSolver
    %OBSTRUCTIONSOLVER Summary of this class goes here
    %
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        
        NZ           (1,1) double  {mustBeNumeric}                          = 0         % [-] Number of axial steps
        NTIME        (1,1) double  {mustBeNumeric}                          = 0         % [-] Number of time steps
        TIME         (:,1) double  {mustBeNumeric}                          = 0         % [s] Time series
        DT           (1,1) double  {mustBeNumeric}                          = 0         % [s] Time step size
        Z            (:,1) double  {mustBeNumeric}                          = 1.        % [m] Elevation
        DZ           (1,1) double  {mustBeNumeric}                          = 0         % [m] Axial step size

        fluid       {isa(fluid,'Inputs.FluidProperties')}
        boundaryConditions
        
        solutionSets (:,1) Solvers.Obstruction.SolutionSet
        axialBounds
        axialBoundIndices
        spanBounds

     end

     properties (SetAccess = protected)
        originalInputset
        inputSet
        STATE                                                               = Solvers.SolverState.UNSOLVED
        solverType
        SOLVERMODE                                                          = Solvers.SolverMode.NEW
     end

    methods
        function obsSolver = ObstructionSolver(inputSet, solverType)
            %OBSTRUCTIONSOLVER Creates a Obstruction solver 
            %  
            %  Detailed explanation goes here

            arguments
                inputSet            {isa(inputSet,'Inputs.InputSet')}
                solverType          InputEnums.SOLVER
            end

            import Solvers.Obstruction.*

            % Call abstract class constructor to initialize solver inputSet
            obsSolver = obsSolver@Solvers.AbstractSolver(inputSet);

            % Check if obstructions are nonempty in inputSet
            if isempty(inputSet.obs)
                error("MIXTURESOLVER:NoObstructionError", "No obstructions are in the inputSet");
            end

            % Save original inputSet
            obsSolver.originalInputset = inputSet;

            % For now, assume only one obstruction is given
            obs = inputSet.obs(1);

            % Wall perimeter
            wallPerim = obsSolver.inputSet.geometry.PERIM(obs.WALL);

            % Place the obstruction on the wall; determine bounds
            [axialBounds, spanBounds] = obs.placeObsOnWall( ...
                                                obsSolver.inputSet.geometry.LENGTH, ...
                                                wallPerim);
            
            obsSolver.axialBounds = [0 axialBounds obsSolver.inputSet.geometry.LENGTH];
            obsSolver.spanBounds = [0 spanBounds wallPerim];

            % Split inputSet into axial segments
            [obsSolver.inputSet, obsSolver.axialBoundIndices] = obsSolver.inputSet.SplitByAxialPosition(obsSolver.axialBounds);

            % Split last/3rd inputSet by tracks
            obsSolver.inputSet(end).SplitBySpanPosition(obs.WALL, obsSolver.spanBounds);

            % Store solver type
            obsSolver.solverType = solverType;
            
            
            
        end
        
        function initializeSolver(obsSolver)
        %INITIALIZESOLVER Initialize solver using the stored inputSet
        %
        % 
            
            import Solvers.Obstruction.SolutionSet
            import Solvers.SolverState            

            % set STATE to INITIALIZED
            obsSolver.STATE = SolverState.INITIALIZED;

            % Log 
            
        end

        function solve(obsSolver)
        %SOLVE
        %
        % Solve the problem

            import Solvers.Obstruction.SolutionSet
            import Solvers.SolverState

            % Solve un-obstructed mixture solver
            mixSolver_unObs = Solvers.Mixture.MixtureSolver(obsSolver.originalInputset);
            mixSolver_unObs.solve();

            % Split the solution into three parts: 
            
            %%
            % 1. Before the obstruction
            obsSolver.solutionSets(1) = SolutionSet( ...
                                            "NONOBS", ...
                                            obsSolver.inputSet(1), ...
                                            obsSolver.axialBounds(1:2));
            % Subset of mixSolver
            subset_StartIdx = obsSolver.axialBoundIndices(1,1);
            subset_length = diff(obsSolver.axialBoundIndices(:,1))+1;
            mixSolver(1) = mixSolver_unObs.subset(subset_StartIdx, subset_length, "inputSet", obsSolver.inputSet(1));
            
            disp(obsSolver.solverType.solverPath()); 
            % Solver for the current step
            stepSolver(1) = Solvers.(obsSolver.solverType.solverPath())(obsSolver.inputSet(1), mixSolver);
            
            % Solve
            stepSolver(1).solve();


            % Store stepSolver to solutionset
            obsSolver.solutionSets(1).setSolver(stepSolver(1));
            

            %%
            % 2. Along the obstruction
            
            % Update inputSet P
            for tIdx = 1:length(obsSolver.inputSet(2).bc)
                obsSolver.inputSet(2).bc(tIdx).setProperty("PRESSURE", stepSolver(1).mixSolver.mixture(1).P(end));
                obsSolver.inputSet(2).bc(tIdx).setProperty("HIN", stepSolver(1).mixSolver.mixture(1).H(end));
            end
            obsSolver.solutionSets(2) = SolutionSet( ...
                                            "OBS", ...
                                            obsSolver.inputSet(2), ...
                                          obsSolver.axialBounds(2:3));
            
            % Subset of mixSolver
            subset_StartIdx = obsSolver.axialBoundIndices(1,2);
            subset_length = diff(obsSolver.axialBoundIndices(:,2))+1;
            mixSolver(2) = mixSolver_unObs.subset(subset_StartIdx, subset_length, "inputSet", obsSolver.inputSet(2));

            % % Set initSolver to true
            % mixSolver(2) = Solvers.Mixture.MixtureSolver(obsSolver.inputSet(2), "solverMode","CONTINUE", "prevSolution", stepSolver(1).mixSolver);
            % 
            % % Solve mixSolver
            % mixSolver(2).solve();
            
            % Create mixture solver if another solver is used
            if obsSolver.solverType ~= InputEnums.SOLVER.MIXTURE                

                % Create solver with previous solution
                stepSolver(2) = Solvers.(obsSolver.solverType.solverPath())(obsSolver.inputSet(2), mixSolver(2),  "solverMode","CONTINUE",  "prevSolution", stepSolver(1));

                % Solve
                stepSolver(2).solve();
                
            else
                % Solve
                %stepSolver(2) = mixSolver(2);
            end

            % Store stepSolver to solutionset
            obsSolver.solutionSets(2).setSolver(stepSolver(2));
            %return

            %%
            %   3. After the obstruction
            % Update inputSet P
            for tIdx = 1:length(obsSolver.inputSet(3).bc)
                obsSolver.inputSet(3).bc(tIdx).setProperty("PRESSURE", stepSolver(2).mixSolver.mixture(1).P(end));
                obsSolver.inputSet(3).bc(tIdx).setProperty("HIN", stepSolver(2).mixSolver.mixture(1).H(end));
            end
            obsSolver.solutionSets(3) = SolutionSet( ...
                                            "NONOBS", ...
                                            obsSolver.inputSet(3), ...
                                            obsSolver.axialBounds(3:4));

            % Subset of mixSolver for post-obstruction segment (use un-obstructed solution grid)
            subset_StartIdx = obsSolver.axialBoundIndices(1,3);
            subset_length = diff(obsSolver.axialBoundIndices(:,3))+1;
            mixSolver(3) = mixSolver_unObs.subset(subset_StartIdx, subset_length, "inputSet", obsSolver.inputSet(3));

            % Create mixture solver if another solver is used
            if obsSolver.solverType ~= InputEnums.SOLVER.MIXTURE                

                % Create solver without init
                stepSolver(3) = Solvers.(obsSolver.solverType.solverPath())(obsSolver.inputSet(3), mixSolver(3),  "solverMode","CONTINUE",  "prevSolution", stepSolver(2));

                % Solve
                stepSolver(3).solve();
                
            else
                % Solve
                stepSolver(3) = mixSolver(3);
            end

            % Store stepSolver to solutionset
            obsSolver.solutionSets(3).setSolver(stepSolver(3));

        end

        function plotter = plott(obsSolver, zIdx,opts)
            warning('plott not yet available within obstruction branch')

        end

        function plotter = plotz(obsSolver, tIdx, opts)
            %PLOTZ Plots spatial (axial distribution) of a selected parameter for each obstruction solver region.
            %
            % Generates axial plots of selected parameters at a specified
            % time indicies
            % The obstruction domain is divided into three regions: before,
            % on, and after the obstruction. This method
            % delegates to each region's solver plotz method.
            %
            % Inputs:
            %
            % - obsSolver     — :class:`Solvers.Obstruction.ObstructionSolver` object
            % - tIdx          — Time index or indices (vector of positive integers)
            % - opt.display   — Parameter(s) to display (delegated to each solver)
            % - opt.label     — Corresponding labels for display parameters
            % - opt.unit      — Units for each parameter
            % - opt.field     — Field to plot (delegated to each solver)
            % - opts.solveMode    — Solve mode: 'REAL' or 'NULL'
            % - opts.wall         — Wall index(es) to plot
            % - opts.zIdx         — Axial indices to include in the plot
            % - opts.annular      — Logical flag to plot annular two-phase flow parameters only from the onset of annular flow
            % - opts.unitTemp     — Temperature unit: 'K' or 'C'
            % - opts.arrangement  — Plot arrangement: 'flow', 'vertical', or 'horizontal'
            % - opts.resize       — Resize factor for plot scaling
            %
            % Notes:
            %
            % - Plotting is delegated to each region's solver (before, on, after obstruction)
            % - Each region's solver determines which parameters can be plotted
            % - Figure titles indicate the obstruction region being plotted

            arguments
                obsSolver
                tIdx              (:,1) double {mustBeInteger,mustBePositive}                                         = []
                opts.display      {matlab.system.mustBeMember(opts.display,{'HFLUX','W','WL','U','THICK','FWE','FME','DME','ALL','VR'})} = {'HFLUX','W','U'}
                opts.solveMode    {mustBeMember(opts.solveMode,{'REAL','NULL'})}                                      = 'REAL'
                opts.wall         (1,:) double {mustBeVector,mustBeInteger,mustBePositive}                            = 1
                opts.zIdx         (:,1) double {mustBeVector,mustBeInteger,mustBePositive}                            = []
                opts.obstructions (1,1) logical                                                                       = false
                opts.annular      (1,1) logical                                                                       = true
                opts.unitTemp     {mustBeMember(opts.unitTemp,{'K','C'})}                                             = 'K'
                opts.arrangement  {mustBeMember(opts.arrangement,{'flow','vertical','horizontal'})}                   = 'flow'
                opts.resize       (1,1) double {mustBeNonnegative}                                                    = 0
                opts.nearWall     (1,1) logical                                                                       = false
            end
            
            %Grab display parameter
            displayParam = opts.display;

            %Surpress Figure display
            set(groot, 'DefaultFigureVisible', 'off')
            %Plot the three field mass results for the three tracks
            obsSolver.solutionSets(1).solver.plotz('display',displayParam);
            obsSolver.solutionSets(2).solver.plotz('display',displayParam, 'oafLine', false, 'annular', false);
            obsSolver.solutionSets(3).solver.plotz('display',displayParam,'oafLine', false, 'annular',false);
            %Restore Figure display
            set(groot, 'DefaultFigureVisible', 'on');
            
            %Find and sort figures to get the last 4
            allFigures = findall(groot, 'Type', 'figure');
            % Filter out invalid/deleted figure handles
            allFigures = allFigures(isvalid(allFigures));
            % Guard: exit early if no valid figures found
            if isempty(allFigures)
                warning('No valid figures found.');
                return;
            end
            figNumbers = [allFigures.Number];
            [~, sortIdx] = sort(figNumbers, 'descend');
            targetFigures = allFigures(sortIdx(end:-1:max(end-3,1)));


            % 1. COLLECT AXIS AND STORE PROPERTIES
            targetAxes = [];
            for i = 1:length(targetFigures)
                ax = findall(targetFigures(i), 'Type', 'axes');
                targetAxes = [targetAxes; ax];
            end
            % Scrape style from first plot
            refAx  = targetAxes(1);
            refFig = ancestor(refAx, 'figure');            
            % Figure-level position
            figPos = get(refFig, 'Position');            
            % Axes-level appearance
            styleFields.LineWidth          = get(refAx, 'LineWidth');
            styleFields.FontName           = get(refAx, 'FontName');
            styleFields.FontSize           = get(refAx, 'FontSize');
            styleFields.FontWeight         = get(refAx, 'FontWeight');
            styleFields.Box                = get(refAx, 'Box');
            styleFields.Color              = get(refAx, 'Color');
            styleFields.XColor             = get(refAx, 'XColor');
            styleFields.YColor             = get(refAx, 'YColor');
            styleFields.XScale             = get(refAx, 'XScale');
            styleFields.YScale             = get(refAx, 'YScale');
            styleFields.TickDir            = get(refAx, 'TickDir');
            styleFields.TickLength         = get(refAx, 'TickLength');           
            % Grid
            styleFields.XGrid              = get(refAx, 'XGrid');
            styleFields.YGrid              = get(refAx, 'YGrid');
            styleFields.GridAlpha          = get(refAx, 'GridAlpha');
            styleFields.GridColor          = get(refAx, 'GridColor');
            styleFields.GridLineStyle      = get(refAx, 'GridLineStyle');
            styleFields.XMinorGrid         = get(refAx, 'XMinorGrid');
            styleFields.YMinorGrid         = get(refAx, 'YMinorGrid');
            styleFields.MinorGridAlpha     = get(refAx, 'MinorGridAlpha');
            styleFields.MinorGridColor     = get(refAx, 'MinorGridColor');
            styleFields.MinorGridLineStyle = get(refAx, 'MinorGridLineStyle');          
            % Axes position/size within figure (normalized units)
            oldUnits = get(refAx, 'Units');
            set(refAx, 'Units', 'normalized');
            axPos  = get(refAx, 'Position');
            axOPos = get(refAx, 'OuterPosition');
            set(refAx, 'Units', oldUnits);
            

            % 2. COLLECT METADATA (titles, labels, axis limits)
            allTitles  = {};
            allXLabels = {};
            allYLabels = {};
            allXLims   = [];
            allYLims   = [];
            
            for i = 1:length(targetAxes)
                ax = targetAxes(i);
            
                t = get(get(ax, 'Title'),  'String');
                x = get(get(ax, 'XLabel'), 'String');
                y = get(get(ax, 'YLabel'), 'String');
            
                if ~isempty(t), allTitles{end+1}  = t; end
                if ~isempty(x), allXLabels{end+1} = x; end
                if ~isempty(y), allYLabels{end+1} = y; end
            
                allXLims = [allXLims; get(ax, 'XLim')];
                allYLims = [allYLims; get(ax, 'YLim')];
            end
            

            % 3. CREATE MASTER FIGURE
            masterFig = figure('Color', get(refFig, 'Color'), 'Position', figPos);
            masterAx  = axes(masterFig);
            set(masterAx, 'Units', 'normalized', 'Position', axPos, 'OuterPosition', axOPos);
            hold(masterAx, 'on');
            % Apply all scraped axis properties
            fieldNames = fieldnames(styleFields);
            for f = 1:length(fieldNames)
                try
                    set(masterAx, fieldNames{f}, styleFields.(fieldNames{f}));
                catch
                    % Skip anything that doesn't apply to axes
                end
            end
            

            % 4. COPY AND ADD LINES FROM FIGURES
            segmentNames = {'Pre-obstruction', 'Along obstruction', 'Free stream', 'Wake'};
            
            legendHandles = [];
            legendLabels  = {};
            seenLabels    = {};
            filmHandles   = containers.Map('KeyType', 'int32', 'ValueType', 'any');
            
            % --- First pass: collect colors already used by non-film lines ---
            takenColors = [];
            for i = 1:length(targetAxes)
                linesInAxis = findobj(targetAxes(i), 'Type', 'line');
                for j = 1:length(linesInAxis)
                    origName = strtrim(get(linesInAxis(j), 'DisplayName'));
                    isFilm   = ~isempty(regexpi(origName, '^film$', 'once'));
                    if ~isFilm
                        takenColors = [takenColors; get(linesInAxis(j), 'Color')];
                    end
                end
            end
            takenColors = unique(takenColors, 'rows');
            
            % --- Pick 4 colors for film segments that don't clash with taken colors ---
            candidateColors = [
                0.00, 0.45, 0.74;
                0.85, 0.33, 0.10;
                0.93, 0.69, 0.13;
                0.49, 0.18, 0.56;
                0.47, 0.67, 0.19;
                0.30, 0.75, 0.93;
                0.64, 0.08, 0.18;
                0.07, 0.62, 0.62;
                0.20, 0.20, 0.80;
                0.80, 0.20, 0.60;
            ];
            
            filmColors = [];
            for k = 1:size(candidateColors, 1)
                c = candidateColors(k,:);
                if ~isempty(takenColors)
                    dists = sqrt(sum((takenColors - c).^2, 2));
                    if min(dists) < 0.3
                        continue;
                    end
                end
                filmColors = [filmColors; c];
                if size(filmColors, 1) == 4
                    break;
                end
            end
            
            % Fallback: pad with lines(4) if not enough distinct colors found
            if size(filmColors, 1) < 4
                filmColors = [filmColors; lines(4 - size(filmColors, 1))];
            end
            
            % --- Scrape film marker style from first film line found ---
            filmMarker          = 'o';
            filmMarkerSize      = 6;
            filmMarkerFaceColor = 'none';
            filmMarkerEdgeColor = 'auto';
            filmLineStyle       = '-';
            
            found = false;
            for i = 1:length(targetAxes)
                if found, break; end
                linesInAxis = findobj(targetAxes(i), 'Type', 'line');
                for j = 1:length(linesInAxis)
                    origName = strtrim(get(linesInAxis(j), 'DisplayName'));
                    if ~isempty(regexpi(origName, '^film$', 'once'))
                        filmMarker          = get(linesInAxis(j), 'Marker');
                        filmMarkerSize      = get(linesInAxis(j), 'MarkerSize');
                        filmMarkerFaceColor = get(linesInAxis(j), 'MarkerFaceColor');
                        filmMarkerEdgeColor = get(linesInAxis(j), 'MarkerEdgeColor');
                        filmLineStyle       = get(linesInAxis(j), 'LineStyle');
                        found = true;
                        break;
                    end
                end
            end
            
            % --- Second pass: copy lines ---
            for i = 1:length(targetAxes)
                linesInAxis = findobj(targetAxes(i), 'Type', 'line');
            
                for j = 1:length(linesInAxis)
                    newLine  = copyobj(linesInAxis(j), masterAx);
                    origName = strtrim(get(linesInAxis(j), 'DisplayName'));
                    isFilm   = ~isempty(regexpi(origName, '^film$', 'once'));
            
                    if isFilm
                        set(newLine, ...
                            'Color',           filmColors(i,:), ...
                            'DisplayName',     '', ...
                            'LineStyle',       filmLineStyle, ...
                            'Marker',          filmMarker, ...
                            'MarkerSize',      filmMarkerSize, ...
                            'MarkerFaceColor', filmMarkerFaceColor, ...
                            'MarkerEdgeColor', filmMarkerEdgeColor);
            
                        if ~isKey(filmHandles, int32(i))
                            filmHandles(int32(i)) = newLine;
                        end
                    else
                        set(newLine, 'DisplayName', '');
            
                        if ~isempty(origName) && ~any(strcmp(seenLabels, origName))
                            seenLabels{end+1}   = origName;
                            legendHandles       = [legendHandles; newLine];
                            legendLabels{end+1} = origName;
                        end
                    end
                end
            end

            %Add OAF dashed line
            % Find the first x position where Film or Drop data begins in segment 1
            onsetX = [];
            linesInSeg1 = findobj(targetAxes(1), 'Type', 'line');
            
            for j = 1:length(linesInSeg1)
                origName = strtrim(get(linesInSeg1(j), 'DisplayName'));
                isFilmOrDrop = ~isempty(regexpi(origName, '^(film|drop)$', 'once'));
            
                if isFilmOrDrop
                    xData = get(linesInSeg1(j), 'XData');
                    yData = get(linesInSeg1(j), 'YData');
            
                    % Find first index where value is non-zero and non-NaN
                    validIdx = find(~isnan(yData) & yData > 0, 1, 'first');
            
                    if ~isempty(validIdx)
                        candidate = xData(validIdx);
                        if isempty(onsetX) || candidate < onsetX
                            onsetX = candidate;  % take the earliest across film and drop
                        end
                    end
                end
            end
            
            % Draw the OAF line if found
            if ~isempty(onsetX)
                yLimsCurrent = [min(allYLims(:,1)), max(allYLims(:,2))];
                line(masterAx, [onsetX onsetX], yLimsCurrent, ...
                    'Color',     [1 0 0], ...
                    'LineStyle', '--', ...
                    'LineWidth', 1.5, ...
                    'DisplayName', '');
            end
            
            % Append film legend entries in segment order
            for i = 1:length(segmentNames)
                if isKey(filmHandles, int32(i))
                    legendHandles       = [legendHandles; filmHandles(int32(i))];
                    legendLabels{end+1} = ['Film - ' segmentNames{i}];
                end
            end
            
            hold(masterAx, 'off');
            

            % 5. APPLY METADATA TO MASTER AXES          
            % Title
            uniqueTitles = unique(allTitles, 'stable');
            if numel(uniqueTitles) == 1
                title(masterAx, uniqueTitles{1});
            elseif numel(uniqueTitles) > 1
                title(masterAx, strjoin(uniqueTitles, ' | '));
            end
            
            % X label
            uniqueXLabels = unique(allXLabels, 'stable');
            if numel(uniqueXLabels) == 1
                xlabel(masterAx, uniqueXLabels{1});
            elseif numel(uniqueXLabels) > 1
                xlabel(masterAx, strjoin(uniqueXLabels, ' / '));
            end
            
            % Y label
            uniqueYLabels = unique(allYLabels, 'stable');
            if numel(uniqueYLabels) == 1
                ylabel(masterAx, uniqueYLabels{1});
            elseif numel(uniqueYLabels) > 1
                ylabel(masterAx, strjoin(uniqueYLabels, ' / '));
            end
            
            % Axis limits spanning all source axes
            if ~isempty(allXLims)
                xlim(masterAx, [min(allXLims(:,1)), max(allXLims(:,2))]);
            end
            if ~isempty(allYLims)
                ylim(masterAx, [min(allYLims(:,1)), max(allYLims(:,2))]);
            end
            
            % Legend
            legend(masterAx, legendHandles, legendLabels, 'Location', 'best');
            
            % 6. CLEAN UP WORKSPACE
            close(targetFigures);
        end

        function plotzt(obsSolver, opt)
            %PLOTZT Plots 2D time/elevation distributions for each obstruction solver region.
            %
            % Generates surface plots of obstruction-related parameters over time 
            % and axial (elevation) positions. The obstruction domain is divided 
            % into three regions: before, on, and after the obstruction. This method
            % delegates to each region's solver plotzt method.
            %
            % Inputs:
            %
            % - obsSolver     — :class:`Solvers.Obstruction.ObstructionSolver` object
            % - opt.display   — Parameter(s) to display (delegated to each solver)
            % - opt.label     — Corresponding labels for display parameters
            % - opt.unit      — Units for each parameter
            % - opt.field     — Field to plot (delegated to each solver)
            % - opt.wall      — Wall index(es) to plot
            % - opt.zIdx      — Axial indices (must include at least 2)
            % - opt.tIdx      — Time indices (must include at least 2)
            % - opt.reverseTime — Logical flag to reverse time axis
            % - opt.shading   — Surface shading style: 'faceted', 'flat', or 'interp'
            % - opt.view      — View angle for 3D plot [azimuth elevation]
            %
            % Notes:
            %
            % - Plotting is delegated to each region's solver (before, on, after obstruction)
            % - Each region's solver determines which parameters can be plotted
            % - Figure titles indicate the obstruction region being plotted

            arguments
                obsSolver
                opt.display      {mustBeA(opt.display,{'cell','char'})}                   = {'HFLUX','W','U','H','X','VF'}
                opt.label        {mustBeA(opt.label,{'cell','char'})}                     = {}
                opt.unit         {mustBeA(opt.unit,{'cell','char'})}                      = {}
                opt.field        {mustBeA(opt.field,{'cell','char'})}                     = 'mixture'
                opt.solveMode    {mustBeMember(opt.solveMode,{'REAL','NULL'})}            = 'REAL'
                opt.wall         (1,:) double {mustBeVector,mustBeInteger,mustBePositive} = 1
                opt.zIdx         (:,1) double {mustBeVector,mustBeInteger,mustBePositive} = []
                opt.tIdx         (:,1) double {mustBeVector,mustBeInteger,mustBePositive} = []
                opt.reverseTime  (1,1) logical                                            = false
                opt.shading      {mustBeMember(opt.shading,{'faceted','flat','interp'})}  = 'interp'
                opt.view         (1,2) double                                             = [0 90]
            end

            % Process options to match solver plotzt interface
            solverOpt = opt;
            if ~isempty(opt.display) && ~iscell(opt.display)
                solverOpt.display = {opt.display};
            end
            if ~isempty(opt.field) && ~iscell(opt.field)
                solverOpt.field = {opt.field};
            end

            % Region names for display
            regionNames = {'Before Obstruction', 'At Obstruction', 'After Obstruction'};

            % Plot each region's solver data
            for regionIdx = 1:length(obsSolver.solutionSets)
                if isempty(obsSolver.solutionSets(regionIdx).solver)
                    obsSolver.log(sprintf('Warning: No solver data available for region %d (%s)\n', ...
                        regionIdx, regionNames{regionIdx}));
                    continue
                end

                % Create a new figure for this region
                regionSolver = obsSolver.solutionSets(regionIdx).solver;
                
                % Call the stepSolver's plotzt method with the provided options
                try
                    regionSolver.plotzt(solverOpt);
                    
                    % Update figure titles to indicate region
                    currentFig = gcf;
                    if ~isempty(currentFig.Name)
                        currentFig.Name = sprintf('%s - %s', regionNames{regionIdx}, currentFig.Name);
                    end
                    
                catch ME
                    obsSolver.log(sprintf('Error plotting region %d (%s): %s\n', ...
                        regionIdx, regionNames{regionIdx}, ME.message));
                end
            end
        end

         function saveResults(ffSolver, opts)
        %SAVERESULTS
        %
        arguments
            ffSolver
            opts.saveFormat {mustBeMember(opts.saveFormat,["MAT"])}   = "MAT"
        end
            session = ffSolver.inputSet.session;
            switch opts.saveFormat
                case "MAT"
                    results = struct( ...
                                'Z', ffSolver.Z, ...
                                'TIME', ffSolver.TIME, ...
                                'sessionName', session.name, ...
                                'boundaryConditions', ffSolver.boundaryConditions, ...
                                'filmInit', struct(ffSolver.filmInit), ...
                                'dropInit', struct(ffSolver.dropInit), ...
                                'film', struct(ffSolver.film), ...
                                'base', struct([ffSolver.film.base]), ...
                                'wave', struct([ffSolver.film.wave]), ...
                                'drop', struct(ffSolver.drop));
                    
                    if isfolder(session.directory)
                        save( ...
                            fullfile( ...
                                session.directory, strcat(session.name,'.mat')), ...
                            '-struct', ...
                            "results", ...
                            "-mat" ...
                        );
                    else
                        error('%s does not exist. Check Session.log.LOGMODE. Try session.makeSessionDirectory()');
                    end
            end


        end

    end
end
