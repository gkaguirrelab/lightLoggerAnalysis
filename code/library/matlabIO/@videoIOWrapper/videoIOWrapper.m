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
        temporary_reading_hdf5_filepath
        current_reading_color_mode 
        read_ahead_buffer
        read_ahead_buffer_size;
        buffer_start_frame
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
            p.addParameter('readAheadBufferSize',1000);
            p.parse(varargin{:})
            options.ioAction = p.Results.ioAction; 
            options.readAheadBufferSize = p.Results.readAheadBufferSize;

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

            % Store the read ahead buffer size for making reads faster 
            obj.read_ahead_buffer_size = options.readAheadBufferSize;

            switch(options.ioAction)
                % In the case of read, load in some information about the video 
                case "read"
                    % First, we need to convert the video to an HDF5 file. This is 
                    % because using the Python wrapper directly on the .avi is slow 
                    % because of the MATLAB interface. We also have to open/close 
                    % the file per frame, which is slow. With HDF5, we will 
                    % generate it quickly with Python, then access it quickly 
                    % with MATLAB, removing it when the object closes 
                    % We will do this conversion when the first frame is read 
                    obj.temporary_reading_hdf5_filepath = fullfile(filepath, name+"_temp.hdf5");

                    % Store the number of frames in the video 
                    obj.NumFrames = double(obj.utility_library.inspect_video_frame_count(videoFileName));

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