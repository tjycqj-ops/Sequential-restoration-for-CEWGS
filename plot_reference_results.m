%% Plot the archived reference results without running the optimization

cewgs_plot_reference_entry(fileparts(mfilename('fullpath')));

function cewgs_plot_reference_entry(package_root)
original_matlab_path = path;
path_cleanup = onCleanup(@()path(original_matlab_path));
addpath(fullfile(package_root,'scripts'));
source_csv = fullfile(package_root,'reference_results', ...
    'four_case_low_pressure_source.csv');
output_dir = fullfile(package_root,'reference_figures');
plot_four_case_source_csv(source_csv,output_dir);
clear path_cleanup
end
