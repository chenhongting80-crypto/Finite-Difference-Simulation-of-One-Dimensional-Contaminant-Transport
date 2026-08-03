function save_figure(fig, filename)
%SAVE_FIGURE Export a MATLAB figure with a compatibility fallback.
[parent_folder, ~, ~] = fileparts(filename);
ensure_directory(parent_folder);

try
    exportgraphics(fig, filename, 'Resolution', 220);
catch
    saveas(fig, filename);
end
close(fig);
end
