classdef (Abstract) AbstractField < matlab.mixin.Copyable
    %ABSTRACTFIELD Base class for field definitions across solvers
    %
    % This abstract class defines shared properties and methods for all field-type
    % classes used in solver implementations. It provides mechanisms for transient
    % data extraction, plotting, memoization, struct conversion, and property copying.

    properties (Abstract=true, SetAccess={?Solvers.AbstractSolver, ?Solvers.AbstractField})

        NZ           (1,1) double  {mustBeNumeric}                         % Number of axial steps [-]
        NTIME        (1,1) double  {mustBeNumeric}                         % Number of time steps [-]
        TIME         (1,1) double  {mustBeNumeric}                         % Time series [s]
        DT           (1,1) double  {mustBeNumeric}                         % Time step size [s]
        TIDX         (1,1) double  {mustBeNumeric}                         % Time step index [-]
        Z            (:,1) double  {mustBeNumeric}                         % Elevation [m]
        ITR          (1,1) struct                                          % Iteration properties

    end

    properties (Access = protected)

        %memoizedFunctions = dictionary();
        % Currently using containers.Map() for MATLAB version
        memoizedFunctions = containers.Map();                              % Memoization store for function handles

    end

    properties (SetAccess = protected, Hidden)

        flowProperties (:,:) cell = {'W','U','H','ITR'}                    % Flow properties used for copying

    end

    methods

        function absField = AbstractField()
            %ABSTRACTFIELD Constructor for AbstractField class

        end

        function paramData = transient(obj, param, subobj, opt)
            %TRANSIENT Extracts transient distribution array for a given parameter
            %
            % Supports wall-wise and axial/time slicing with optional
            % subobject access.

            arguments
                obj
                param         (1,1) string {mustBeTextScalar}
                subobj                                                                 = []
                opt.wall      (1,:) double {mustBeVector,mustBeInteger,mustBePositive} = []
                opt.zIdx      (:,1) double {mustBeVector,mustBeInteger,mustBePositive} = 1:obj(1).NZ
                opt.tIdx      (:,1) double {mustBeVector,mustBeInteger,mustBePositive} = 1:length(obj)
            end

            % Read parameter
            p = split(param,'.');
            if isscalar(p)
                if isempty(subobj)
                    paramData = cell2mat(arrayfun(@(x) x.(param),obj(opt.tIdx),'uni',0));
                else
                    paramData = cell2mat(arrayfun(@(x,y) x.(param)(y),obj(opt.tIdx),subobj(opt.tIdx),'uni',0));
                end
            else
                if isempty(subobj)
                    paramData = cell2mat(arrayfun(@(x) x.(p{1}).(p{2}),obj(opt.tIdx),'uni',0));
                else
                    paramData = cell2mat(arrayfun(@(x,y) x.(p{1}).(p{2})(y),obj(opt.tIdx),subobj(opt.tIdx),'uni',0));
                end
            end

            % Size parameter based on options
            NWALL = size(paramData,2)/length(opt.tIdx);
            if isempty(opt.wall), opt.wall = 1:NWALL; end
            if NWALL > 1
                idx = cell2mat(arrayfun(@(n) [n:NWALL:size(paramData,2)],opt.wall,'uni',0));
                paramData = paramData(:,idx);                              % Wall discretization
            end
            paramData = paramData(opt.zIdx,:);                             % Keep relevant axial length
            paramData = reshape(paramData,length(opt.zIdx),length(opt.tIdx),length(opt.wall));

            % Special case for single elevation
            if isscalar(opt.zIdx)
                paramData = permute(paramData,[3 2 1]);
            end
        end

        function ax = plotzt(obj, param ,ylabelText ,ylabelUnit ,k ,opt, subobj, annular,time)
            %PLOTZT 2D Generates a 2D space-time distribution plot for a given parameter

            % Supports unit conversion, pre-annular flow masking, and
            % customizable view options.

            if nargin < 8, annular = false; end
            if nargin < 9, time = []; end

            z     = obj(1).Z(opt.zIdx);                                    % [m]
            if isempty(time), time = [obj(opt.tIdx).TIME];  end            % [s]
            if opt.reverseTime
                time = time -time(end);
            end

            % Load data
            try
                paramData = obj.transient(param,       'wall',k,'zIdx',opt.zIdx,'tIdx',opt.tIdx);
            catch
                paramData = obj.transient(param,subobj,'wall',k,'zIdx',opt.zIdx,'tIdx',opt.tIdx);
            end

            % Adjust for temperature unit
            if strcmp(ylabelUnit,'C')
                paramData = paramData-273.15;
            end

            % Remove pre-annular flow region
            if annular
                try
                    OAFIDX = arrayfun(@(x) x.OAFIDX,obj);
                catch
                    OAFIDX = arrayfun(@(x) x.mix.OAFIDX,obj);
                end
                OAFIDX = OAFIDX - opt.zIdx(1) + 1;
                for k = 1:length(OAFIDX)
                    paramData(1:OAFIDX(k),k) = nan;
                end
            end

            ax = nexttile; hold on; grid on; title(ylabelText)
            [t_mesh,z_mesh] = meshgrid(time,z);
            surf(z_mesh,t_mesh,paramData,'edgeColor','none');
            xlabel('Axial position [m]'); xlim([min(z)       max(z)]);
            ylabel('Time [s]')          ; ylim([min(time) max(time)]);
            cb = colorbar(); cb.Label.String = [ylabelText ' [' ylabelUnit ']']; cb.Label.FontSize = 14;
            set(gca,'fontSize',14)
            shading(opt.shading)
            view(opt.view);
        end

        function out = memoizeFunction(obj, methodStr, methodHandle, varargin)
            %MEMOIZEDMETHOD Registers and retrieves memoized function handles
            %
            % Avoids redundant computation by caching results.
            %
            % Adapted from https://stackoverflow.com/a/75037451

            %For the first call with a particular method, create and
            %memoize a function handle view of the method
            if ~isConfigured(obj.memoizedFunctions) || ~obj.memoizedFunctions.isKey(methodStr)
                fn_method = @(varargin)methodHandle(varargin{:});
                fn = memoize(fn_method);
                obj.memoizedFunctions(methodStr) = fn;
            end

            %For all calls, get the store function handle out of
            %storage, and use it.
            fn = obj.memoizedFunctions(methodStr);
            out = fn(varargin{:});
        end

        function out = struct(obj)
            %STRUCT Converts field object to a structured array
            %
            % Includes TIME, ITR, and flow properties.

            flowProps = obj.flowProperties;
            for i = length(obj):-1:1

                % Add `TIME` and `ITR` by default
                outElement = struct('TIME', obj(i).TIME, ...
                    'ITR', obj(i).ITR);

                % Add flow properties as specified
                for flowPropIdx = 1:length(flowProps)

                    % Name of flow property
                    flowProp = flowProps{flowPropIdx};

                    % Set flowProp as new field
                    outElement.(flowProp) = obj(i).(flowProp);

                end

                % add outElement to out
                out(i) = outElement;
            end
        end
        
        function copyFlowProperties(srcObj, targetObj, opts)
            %COPYFLOWPROPERTIES Copies flow properties from source to target field object
            %
            % Supports full or partial copying depending on 'opts.copyMode'.
            %
            % Copy modes:
            %   "rest"     - Partial copy, preserving inlet conditions (default)
            %   "full"     - Full copy of all flow properties
            %   "continue" - Copies last spatial element of src to first of target
            %
            % Deprecated:
            %   opts.all (logical) - Use opts.copyMode="full" instead

            arguments
                srcObj
                targetObj (1,:) Solvers.AbstractField
                opts.copyMode (1,1) string {mustBeMember(opts.copyMode, {'full','rest','continue'})} = "rest"
                opts.all      (1,1) logical = false  % deprecated, use copyMode="full"
            end

            % Resolve legacy opts.all into copyMode
            if opts.all
                warning( ...
                    'AbstractField:copyFlowProperties:deprecatedOption', ...
                    'opts.all is deprecated and will be removed in a future version. Use copyMode="full" instead.' ...
                    );
                opts.copyMode = "full";
            end

            % TODO: add type check for srcObj and targetObj

            for i = 1:length(targetObj)

                % Mesh check (not required for "continue" as src and target
                % may intentionally have different meshes)
                if opts.copyMode == "full" || opts.copyMode == "rest"
                    if srcObj.Z ~= targetObj(i).Z
                        classType = class(srcObj);
                        throw( ...
                            MException( ...
                                'AbstractFieldError:copyFlowPropertiesError', ...
                                sprintf('Source and target objects (%s) have mismatched spatial meshes', classType) ...
                                ) ...
                            );
                    end
                end

                % Copy properties
                propNames = srcObj.flowProperties;
                for j = 1:length(propNames)

                    if opts.copyMode == "full"
                        % Full copy
                        targetObj(i).(propNames{j}) = srcObj.(propNames{j});

                    elseif opts.copyMode == "rest"
                        % Partial copy, preserving inlet conditions (index 1)
                        if isstruct(targetObj(i).(propNames{j})) && isscalar(targetObj(i).(propNames{j}))
                            % Scalar structs are copied per field
                            structFields = fieldnames(targetObj(i).(propNames{j}));
                            for ii = 1:length(structFields)
                                targetObj(i).(propNames{j}).(structFields{ii})(2:end) = ...
                                    srcObj.(propNames{j}).(structFields{ii})(2:end);
                            end
                        else
                            % Non-scalar properties copied as a vector
                            targetObj(i).(propNames{j})(2:end) = srcObj.(propNames{j})(2:end);
                        end

                    elseif opts.copyMode == "continue"
                        % Delegate to overridable hook method to allow subclass-specific
                        % behaviour (e.g. wall-splitting) without modifying this function
                        
                        % Scalar structs are copied per field
                        if isstruct(targetObj(i).(propNames{j})) && isscalar(targetObj(i).(propNames{j}))
                            structFields = fieldnames(targetObj(i).(propNames{j}));
                            for ii = 1:length(structFields)
                                targetObj(i).(propNames{j}).(structFields{ii})(1) = ...
                                    srcObj(i).(propNames{j}).(structFields{ii})(end);
                            end
                        % Non-scalar properties are copied as a vector
                        else
                            % If property size shows different num. of walls, 
                            % look at inputset.obs for hints, FOR NOW
                            % TODO: if there are more than 1 obstruction,
                            %       major changes will be needed.
                            if size(targetObj(i).(propNames{j}), 2) ~= size(srcObj(i).(propNames{j}), 2)

                                % TODO: Check if an obs exists
                                
                                % Determine obstruction wall id
                                wallID = srcObj.inputSet.obs(1).WALL;

                                % source value
                                srcVal = srcObj(i).(propNames{j})(end,:);

                                % Split srcVal at wallID to 2
                                % ex. if wallID ==1 , targetVal(:,[1,2])
                                % will correspond to srcVal(:,1)
                                if propNames{j} == 'W'
                                    massSplitRatio = srcObj.inputSet.model.OBSWSPLITRATIO;
                                    massSplitRatio = massSplitRatio.'./sum(massSplitRatio);
                                    
                                    targetVal = [srcVal(:,1:wallID-1), repmat(srcVal(:,wallID),1,2).*massSplitRatio, srcVal(:,wallID+1:end)];
                                else
                                    targetVal = [srcVal(:,1:wallID-1), repmat(srcVal(:,wallID),1,2), srcVal(:,wallID+1:end)];
                                end

                                % Assign targetVal
                                targetObj(i).(propNames{j})(1,:) = targetVal;                                

                                
                            else
                                if isobject(srcObj(i).(propNames{j}))
                                    % copy flow properties
                                    srcObj(i).(propNames{j})(end).copyFlowProperties(targetObj(i).(propNames{j})(1),"copyMode","continue");
                                else
                                    % simply copy if same size
                                    targetObj(i).(propNames{j})(1,:) = srcObj(i).(propNames{j})(end,:);
                                end
                            end
                        end
                    end
                end
            end
        end

    end

end

