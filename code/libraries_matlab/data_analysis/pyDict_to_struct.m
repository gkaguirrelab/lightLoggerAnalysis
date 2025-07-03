function s = pyDict_to_struct(d)

keys = d.keys;
    s = struct();
    for i = 1:double(py.len(keys))
        key = string(keys{i});
        val = d{keys{i}};
        
        % If nested dict, convert recursively
        if isa(val, 'py.dict')
            s.(key) = pyDict_to_struct(val);
        % If numpy array, convert to double
        elseif isa(val, 'py.numpy.ndarray')
            s.(key) = double(val);
        else
            s.(key) = val;
        end
    end
end

%% SCRIPT TO USE
file_handle = py.open('/Users/sophiamirabal/Documents/MATLAB/test.pkl', 'rb');
pkl = py.importlib.import_module('pickle');
chunks_as_py = pkl.load(file_handle);
file_handle.close()

chunks = cell(chunks_as_py);
chunks = cellfun(@(x) chunk_dict_to_matlab(x), chunks, 'UniformOutput', false);


