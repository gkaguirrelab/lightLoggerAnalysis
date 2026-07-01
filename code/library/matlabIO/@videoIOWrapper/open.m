function open(obj)
% Prepare filesystem state for buffered video writing.
%
% Syntax:
%   open(obj)
%
% Description:
%   This method creates the temporary frame directory used by
%   `videoIOWrapper` when it is operating in write mode. Each subsequent
%   call to `writeVideo` saves a numbered image into that directory, and
%   the final AVI is assembled from those images when `close` is called.
%   In read mode the method intentionally does nothing because the HDF5
%   read buffer is created lazily on the first `readFrame` call.
%
% Inputs:
%   obj                      - `videoIOWrapper` instance.
%
% Outputs:
%   None.
%
% Examples:
%{
    writer = videoIOWrapper("output.avi", "ioAction", "write");
    open(writer);
%}

    arguments 
        obj 
    end     

    % If we are reading and not writing, just do nothing 
    if(~strcmp(obj.mode, 'write'))
        return; 
    end     

    % Construct the path to the temporary file where we will write 
    % out frames to disk 
    if(strcmp(obj.mode, 'write'))
        temp_dir_path = fullfile(obj.Path, obj.filename); 
        
        if ~exist(temp_dir_path, 'dir')
            mkdir(temp_dir_path);
        end
    end
end
