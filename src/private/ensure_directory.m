function ensure_directory(folder_path)
%ENSURE_DIRECTORY Create a directory when it does not already exist.
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end
end
