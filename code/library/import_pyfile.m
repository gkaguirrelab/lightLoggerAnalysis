function module = import_pyfile(path)
% Import a Python module directly from an explicit file path.
%
% Syntax:
%   module = import_pyfile(path)
%
% Description:
%   This function loads a Python source file by absolute path instead of by
%   module basename. It creates a unique synthetic module name for every
%   import so multiple different files with the same basename can coexist
%   in one MATLAB/Python session, and it temporarily ensures the source
%   directory is present on `py.sys.path` so sibling imports inside that
%   module still resolve correctly. This avoids the common failure mode in
%   which Python reuses a cached module with the right name but from the
%   wrong project directory.
%
% Inputs:
%   path                     - String or char vector. Full path to the
%                              Python `.py` file to import.
%
% Outputs:
%   module                   - Python module object corresponding to the
%                              requested file.
%
% Examples:
%{
    module = import_pyfile("/path/to/world_util.py");
%}

    arguments
        path {mustBeText}; % The path to the Python file we want to import
    end

    % Resolve the incoming path so we load this exact file rather than
    % whichever same-named module Python may already have cached. Importing
    % by basename (py.importlib.import_module) returns whatever module is
    % already cached in sys.modules under that name, which silently loads the
    % wrong file when several projects define e.g. world_util.py / Pi_util.py.
    path = string(path);
    [target_dir, filename, ext] = fileparts(path);
    assert(ext == ".py", "import_pyfile:ExpectedPythonFile", ...
        "Expected a .py file, got %s", path);

    % Give each load a unique synthetic module name so multiple different
    % files sharing a basename can coexist in one MATLAB/Python session.
    module_name = matlab.lang.makeValidName(filename + "_" + string(char(java.util.UUID.randomUUID)));

    % Some utility files import siblings from the same folder, so ensure the
    % target directory is importable while we execute this module.
    target_dir = char(target_dir);
    sys_path = py.sys.path;
    path_count = py.getattr(sys_path, 'count');
    if int64(path_count(target_dir)) == 0
        path_insert = py.getattr(sys_path, 'insert');
        path_insert(int32(0), target_dir);
    end

    % Import the module directly from the requested file path.
    spec = py.importlib.util.spec_from_file_location(module_name, char(path));
    module = py.importlib.util.module_from_spec(spec);
    setitem = py.getattr(py.sys.modules, '__setitem__');
    setitem(module_name, module);
    spec.loader.exec_module(module);

    % Return the imported module
    return ;


end
