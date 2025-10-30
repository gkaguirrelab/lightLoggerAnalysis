% Object to support reading AVI video files using python routines

classdef videoIOWrapper < handle

    properties (Constant)

        myConstant = 8;
    end

    % Private properties
    properties (GetAccess=private)
        utility_library

    end

    % Calling function can see, but not modify
    properties (SetAccess=private)

        Name
        Path
        NumFrames
        Duration
        FrameRate
        Width
        Height

    end

    % These may be modified after object creation
    properties (SetAccess=public)

        % Verbosity
        verbose = false;

    end

    methods

        % Constructor
        function obj = videoIOWrapper(videoFileName,options)

            arguments
                videoFileName
                options.ioAction char = 'read'
            end

            % Import the Python helper library 
            obj.utility_library = import_pyfile(getpref("lightLoggerAnalysis", "Pi_util_path")); 

            % Store the name and path to the video
            [filepath,name,ext] = fileparts(videoFileName);
            obj.Name = [name,ext];
            obj.Path = filepath;

            % Store the number of frames in the video 
            obj.NumFrames =  obj.utility_library.inspect_video_frame_count(videoFileName);

            % Retrieve the FPS of the video
            obj.FrameRate =  obj.utility_library.inspect_video_FPS(videoFileName); 

            % Set the duration of the video 
            obj.Duration = (obj.NumFrames / FrameRate); 

            % Retrieve the frame size of the video 
            frame_size = cell(obj.utility_library.inspect_video_framesize(videoFileName)); 
            rows = double(frame_size{1});
            cols = double(frame_size{2});
            
            % Set the dimensions of the video
            obj.Height = rows; 
            obj.Width = cols; 

        end

        % Required methds
        frame = read(obj,frameNum)

    end
end