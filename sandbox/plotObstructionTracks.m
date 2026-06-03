function plotObstructionTracks(obsSolver, segmentIdx, tIdx, zTarget, varargin)
%PLOTOBSTRUCTIONTRACKS Film-only plotting for two post-obstruction tracks
%   plotObstructionTracks(obsSolver)
%plots segment 3, tIdx=1, middle z; plots film W and U for the two
%tracks side-by-side (columns = tracks, rows = axial & time series)
%   plotObstructionTracks(obsSolver, segmentIdx, tIdx, zTarget)
%       - specify segment, time index and axial z for the time series
%
% Notes:
% - This helper uses only film fields (film.W and film.U). No droplet
%   fields are accessed or plotted.
% - The two tracks are assumed to be at WALL and WALL+1 where WALL is the
%   original obstructing wall index stored in obsSolver.originalInputset.obs(1).WALL

    % defaults
    if nargin < 2 || isempty(segmentIdx)
        segmentIdx = 3;
    end
    if nargin < 3 || isempty(tIdx)
        tIdx = 1;
    end

    % Retrieve solution set
    try
        solSet = obsSolver.solutionSets(segmentIdx);
    catch
        error('plotObstructionTracks:InvalidSegment', 'Segment %d does not exist in obsSolver.solutionSets', segmentIdx);
    end
    solver = solSet.solver;

    % require film
    if ~isprop(solver, 'film')
        error('plotObstructionTracks:NoFilmField', 'Selected segment solver does not have film fields. Use a ThreeField solver for this helper.');
    end

    % original obstructing wall index -> tracks at WALL and WALL+1
    wallID = obsSolver.originalInputset.obs(1).WALL;
    trackCol1 = wallID;
    trackCol2 = wallID + 1;

    % axial coords
    z = solver.Z;

    % clamp tIdx
    tIdx = max(1, min(tIdx, length(solver.TIME)));

    % extract film snapshot for axial plots
    film = solver.film(tIdx);

    % basic size checks
    nwalls = size(film.W,2);
    if trackCol1 > nwalls || trackCol2 > nwalls
        error('plotObstructionTracks:TrackColumns', 'Computed track columns (%d,%d) exceed available walls (%d). Inspect inputSet.geometry.PERIM.', trackCol1, trackCol2, nwalls);
    end

    % axial arrays (film-only)
    track1_W_axial = film.W(:, trackCol1);
    track2_W_axial = film.W(:, trackCol2);

    track1_U_axial = film.U(:, trackCol1);
    track2_U_axial = film.U(:, trackCol2);
   
    % Plot: 2 rows x 2 cols -> columns are tracks
    figure('Name', sprintf('Obstruction tracks (film-only): segment %d', segmentIdx));
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % Top-left: Mass Flow track 1
    ax11 = nexttile(1);
    plot(ax11, z, track1_W_axial, '-o', 'DisplayName', sprintf('W (film) track %d', trackCol1)); hold(ax11,'on');
    ylabel(ax11,'W [kg/s]');
    xlabel(ax11,'Axial position z [m]');
    title(ax11, sprintf('Track %d (non-wake) — axial (t=%0.3f s)', trackCol1, solver.TIME(tIdx)));
    legend(ax11,'Location','best');
    grid(ax11,'on');

    % Top-right: Mass Flow track 2
    ax12 = nexttile(2);
    plot(ax12, z, track2_W_axial, '-o', 'DisplayName', sprintf('W (film) track %d', trackCol2)); hold(ax12,'on');
    ylabel(ax12,'W [kg/s]');
    xlabel(ax12,'Axial position z [m]');
    title(ax12, sprintf('Track %d (wake) — axial (t=%0.3f s)', trackCol2, solver.TIME(tIdx)));
    legend(ax12,'Location','best');
    grid(ax12,'on');

    %{
    %Bottom-left: Velcotiy series track 1 
    ax21 = nexttile(3);
    plot(ax21, z, track1_U_axial, '-o', 'DisplayName', sprintf('U (film) track %d', trackCol1)); hold(ax21,'on');
    ylabel(ax21,'U [m/s]');
    xlabel(ax21,'Axial position z [m]');
    title(ax21, sprintf('Track %d (non-wake) — axial (t=%0.3f s)', trackCol1, solver.TIME(tIdx)));
    legend(ax21,'Location','best');
    grid(ax21,'on');

    % Bottom-right: time series track 2 at zIdx
    ax22 = nexttile(4);
    plot(ax22, z, track2_U_axial, '-o', 'DisplayName', sprintf('U (film) track %d', trackCol2)); hold(ax22,'on');
    ylabel(ax22,'U [m/s]');
    xlabel(ax22,'Axial position z [m]');
    title(ax22, sprintf('Track %d (wake) — axial (t=%0.3f s)', trackCol2, solver.TIME(tIdx)));
    legend(ax22,'Location','best');
    grid(ax22,'on');

    % link time series x axis
    linkaxes([ax21, ax22],'x');
    
    %}
end
