function delete(obj)
% Delete overload that ensures temporary video resources are cleaned up.
%
% Syntax:
%   delete(obj)
%
% Description:
%   MATLAB calls this method when the wrapper object is cleared or goes out
%   of scope. The implementation delegates to `close` so both explicit and
%   implicit teardown paths share the same cleanup logic for HDF5 buffers
%   and temporary frame directories.
%
% Inputs:
%   obj                      - `videoIOWrapper` instance being destroyed.
%
% Outputs:
%   None.
%
% Examples:
%{
    delete(reader);
%}

    obj.close();
end 
