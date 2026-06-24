classdef Obstruction < Inputs.Input
    %OBSTRUCTION Reuseable wall obstruction definition
    %   Detailed explanation goes here
    
    properties (SetAccess=?Inputs.Input)
        
        ID         (1,1) string  {mustBeTextScalar}                                                % Channel ID
        SHAPE      (1,1) InputEnums.OBSSHAPE                                                       % Base geom. shape
        WALL       (1,1) double  {mustBePositive,mustBeNonempty}           = 1                     % Wall index as defined by Geometry
        RLOCATION  (2,1) double  {mustBePositive,mustBeNonempty}           = [0.5 0.5]             % Relative location on wall (axial, span)

        % CIRCLE GEOM
        DIAMETER   (1,1) double  {mustBePositive,mustBeNonempty}           = 1                     % Diameter [m]

        % TRIANGLE GEOM
        SIDELENGTH (3,1) double  {mustBePositive,mustBeNonempty}           = [1 1 1].*1E-3         % Side lengths [m] (from pointing vertex, CW)
        POINTANGLE (1,1) double  {mustBeNumeric,mustBeNonempty}            = 0                     % Point direction against flow direction (0) [deg]

        %WALL PROPERTIES (for heat transfer model)
        WALLSHAPE  (1,1) InputEnums.WALLSHAPE                              = 'PLANAR'              % Wall Geometry shape, either 'PLANAR' or 'CYLINDRICAL'
        WALLTHICK  (1,1) double  {mustBePositive,mustBeNonempty}           = 0.003                 % Thickness of wall [m]

        % Calculated properties
        LENGTH     (1,1) double  {mustBePositive,mustBeNonempty}           = 1                     % Axial length [m]
        WIDTH      (1,1) double  {mustBePositive,mustBeNonempty}           = 1                     % Spanwise width [m]
        AREA       (1,1) double  {mustBePositive,mustBeNonempty}           = 1                     % Coolant area [m^2] 
        PERIM      (1,:) double  {mustBePositive,mustBeNonempty}           = 1                     % Perimeters [m]
        NWALL      (1,1) double  {mustBePositive,mustBeNonempty}           = 1                     % Number of walls [-] TODO: consider removing this property
        
    end

    methods
        function obj = Obstruction(filePath,obstructionID)
            %MODEL Construct an instance of this class
            %
            %   Detailed explanation goes here
            arguments
                filePath = ""
                obstructionID = ""
            end

            % Call superclass constructor to parse file and select
            % specified optionsID using "ID" key
            obj = obj@Inputs.Input(filePath, 'ID', obstructionID);

            % Return default value if empty inputs are given
            if strlength(filePath) == 0
                obj.ID = "DEFAULT";
                return
            end
            
            %
            % List of immutable obj property names
            objPropnames = obj.listInputProperties("exclude",{"LENGTH", "WDITH", "AREA", "PERIM", "NWALL"});
            %objPropnames = obj.listInputProperties()
            
            % Array of fieldnames using default values
            defaultValueFieldNames = string().empty();
            defaultValues = {};
            
            % Iterate through obj property names
            for idx = 1:length(objPropnames)

                % Retrieve idx-th item in objPropnames
                objPropname = objPropnames(idx);

                % Check if the objPropname entry is specified, and if the
                % default value should be used
                [isSpecified, useDefault, defaultValue] = obj.validateInputEntry(objPropname,id=obstructionID);

                if ~useDefault

                    % Take care of special cases
                    % Function handles provided in string format cannot be
                    % automatically cast to a function_handle. Here, a
                    % validation is first performed to detect restricted
                    % keywords, then converted.
                    if isa(obj.(objPropname), "function_handle") && isstring(obj.inputStruct.(objPropname))

                        % Check for insecure keywords in function handle
                        obj.validateFunctionHandleInput(obj.inputStruct.(objPropname))    ;
                        obj.(objPropname) = ...
                            str2func(obj.inputStruct.(objPropname));
                    else
                        obj.(objPropname) = ...
                            upper(obj.inputStruct.(objPropname));
                    end
                elseif useDefault
                    defaultValueFieldNames(end+1) = objPropname;
                    defaultValues{end+1} = defaultValue;
                end

                if isSpecified
                    % Remove objPropname from inputStruct
                    obj.inputStruct = rmfield(obj.inputStruct, objPropname);
                end
            end

            % Default value used warning
            if ~isempty(defaultValueFieldNames)
                defaultValueWarningString = obj.defaultValueUsedReport(defaultValueFieldNames, defaultValues);
                if nargout == 0
                    warning('Model:defaultValueUsedWarning', ...
                        sprintf('%s\n',defaultValueWarningString));
                else
                    w = struct('warnID', 'Model:defaultValueUsedWarning', ...
                        'msg', defaultValueWarningString);
                    if isempty(obj.warnings)
                        obj.warnings = w;
                    else
                        obj.warnings(end+1) = w;
                    end
                end
            end

            % If extra fields in obj.inputStruct remain, warn user
            remainingInputStructFields = fieldnames(obj.inputStruct);
            if ~isempty(remainingInputStructFields)
                warning( ...
                    '%s: These entries were not used: \n\t %s ', ...
                    upper(class(obj)), sprintf('%s ',remainingInputStructFields{:}) ...
                    );
                obj.extra = obj.inputStruct;
            end

            % Remove dynamic property inputStruct
            inputStructProp = obj.findprop('inputStruct');
            delete(inputStructProp)
            
            % Calculate Length, Area, Perim, NWall
            switch (obj.SHAPE)
                case InputEnums.OBSSHAPE.CIRCLE
                    obj.LENGTH = obj.DIAMETER;
                    obj.WIDTH = obj.DIAMETER;
                    obj.AREA = 0.25.*pi.*obj.DIAMETER.^2;
                    obj.PERIM = pi.*obj.DIAMETER;
                    obj.NWALL = 1;
                otherwise
                    error("TWOPHASESOLVER:errorObstruction:invalidShape", ...
                        "Obstruction shape invalid or not yet implemented");
            end

        end

        function dh = HDIAM(obj)
            % HDIAM Hydraulic diameter
            %
            %   TODO: Verfiy this is correct and typical for flow obs.
            
            dh = 4*obj.AREA/sum(obj.PERIM);

        end

        function [axialBounds, spanBounds] = placeObsOnWall(obj, wallLength, wallWidth)

            % Assume RLOCATION refers to the symmetrical width/length of
            % the obstruction
            axialBounds = wallLength.*obj.RLOCATION(1)+(obj.LENGTH./2).*[-1, 1];
            axialBounds(axialBounds<0) = 0;
            spanBounds = wallWidth.*obj.RLOCATION(2)+(obj.WIDTH./2).*[-1, 1];

        end

        

    end
    
    methods (Static)
        function writeInputFile(filePathName, ID, varargin)
            Inputs.Input.writeInputFile(filePathName, "a+", "ID", ID, varargin{:});
        end
    end

end
