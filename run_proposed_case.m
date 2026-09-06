%% Run the proposed CEWGS restoration method
% This entry runs the proposed pressure-driven transient model only.

cewgs_run_proposed_entry(fileparts(mfilename('fullpath')));

function cewgs_run_proposed_entry(package_root)
original_matlab_path = path;
path_cleanup = onCleanup(@()path(original_matlab_path));
addpath(fullfile(package_root,'scripts'));
run_four_case_comparison_impl(package_root,true,false,{'proposed'});
clear path_cleanup
end
