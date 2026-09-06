%% Run all four CEWGS comparison cases and generate result figures

cewgs_run_comparison_entry(fileparts(mfilename('fullpath')),false);

function cewgs_run_comparison_entry(package_root,skip_plot)
original_matlab_path = path;
path_cleanup = onCleanup(@()path(original_matlab_path));
addpath(fullfile(package_root,'scripts'));
run_four_case_comparison_impl(package_root,skip_plot,false);
clear path_cleanup
end
