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
        Width
        Height
        full_video_path  
        filename
        last_frame_read = 0; 
        last_frame_written = 0; 
        mode;

    end

    % These may be modified after object creation
    properties (SetAccess=public)

        % Verbosity
        verbose = false;
        
        % Public so that we can set it when writing
        FrameRate

    end

    methods

        % Constructor
        function obj = videoIOWrapper(videoFileName,varargin)

            % Using an old-fashioned input parser as the arguments block
            % syntax fails under Matlab 2025a.
            p = inputParser; p.KeepUnmatched = false;            
            p.addParameter('ioAction','read',@ischar);
            p.parse(varargin{:})
            options.ioAction = p.Results.ioAction;

            % Import the Python helper library 
            obj.utility_library = import_pyfile(getpref("lightLoggerAnalysis", "video_io_util_path")); 

            % Store the name and path to the video
            [filepath,name,ext] = fileparts(videoFileName);
            obj.Name = [name,ext];
            obj.Path = filepath;
            obj.full_video_path = videoFileName; 
            obj.filename = name; 
            
            % Store the mode this was opened with 
            obj.mode = options.ioAction; 

            switch(options.ioAction)
                % In the case of read, load in some information about the video 
                case "read"
                    % Store the number of frames in the video 
                    obj.NumFrames =  double(obj.utility_library.inspect_video_frame_count(videoFileName));

                    % Retrieve the FPS of the video
                    obj.FrameRate = double(obj.utility_library.inspect_video_FPS(videoFileName)); 

                    % Set the duration of the video 
                    obj.Duration = (obj.NumFrames / obj.FrameRate); 

                    % Retrieve the frame size of the video 
                    frame_size = cell(obj.utility_library.inspect_video_framesize(videoFileName)); 
                    rows = double(frame_size{1});
                    cols = double(frame_size{2});
                    
                    % Set the dimensions of the video
                    obj.Height = rows; 
                    obj.Width = cols; 

                case "write"

                otherwise 
                    error("Unsupported action")

            end 

        end

        % Required methds
        % FileIO 
        open(obj); 
        close(obj);

        % video I/O 
        frame = readFrame(obj,frameNum, options)
        writeVideo(obj,frameNum, options)

        % Delete/clear overload 
        delete(obj); 

        

    end
end