function s = load_h5_struct(filename, groupname)
% Recursively load a group from HDF5 as a MATLAB struct
%
% filename: path to the .h5 file
% groupname: path to the group (must start with '/'), e.g., '/Params'

info = h5info(filename, groupname);
s = struct();

% read all datasets in this group
for i = 1:length(info.Datasets)
    dname = info.Datasets(i).Name;           % dataset name, e.g., 'g'
    dataset_path = [groupname '/' dname];    % full HDF5 path
    data = h5read(filename, dataset_path);   % read dataset
    
    % convert string datasets (uint8 arrays) to MATLAB char
    if ischar(data)
        s.(dname) = data;
    elseif isnumeric(data) || islogical(data)
        s.(dname) = data;
    elseif iscell(data) && all(cellfun(@isnumeric,data))
        % sometimes strings stored as cell arrays of uint8
        s.(dname) = char(cell2mat(data)');
    else
        s.(dname) = data;
    end
end

% recursively read subgroups
for i = 1:length(info.Groups)
    gname = info.Groups(i).Name;                 % full path, e.g., '/Params/field'
    [~, fieldname] = fileparts(gname);          % extract the last part as field name
    s.(fieldname) = load_h5_struct(filename, gname);
end
end