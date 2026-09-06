%% Validate the public package without invoking YALMIP or Gurobi

package_root = fileparts(mfilename('fullpath'));
required_files = {
    fullfile(package_root,'model','main_cewgs_restoration.m')
    fullfile(package_root,'scripts','run_four_case_comparison_impl.m')
    fullfile(package_root,'data','cewgs_33_20_15_test_system.mat')
    fullfile(package_root,'reference_results','four_case_low_pressure_source.csv')
    fullfile(package_root,'reference_results','four_case_low_pressure_summary.csv')
    };
assert(all(cellfun(@(f)exist(f,'file')==2,required_files)), ...
    'The public package is incomplete.');
original_matlab_path = path;
path_cleanup = onCleanup(@()path(original_matlab_path));
addpath(fullfile(package_root,'model'),fullfile(package_root,'tests'));
test_gas_source_pressure_boundaries;
fprintf('PASS: public package structure and test-system data are valid.\n');
clear path_cleanup
