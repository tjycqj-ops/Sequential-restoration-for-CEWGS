function run_four_case_comparison_impl(package_root,skip_plot,validate_only,selected_case_keys,environment_overrides)
%RUN_FOUR_CASE_COMPARISON_IMPL Execute the four CEWGS comparison cases.

if nargin < 2
    skip_plot = false;
end
if nargin < 3
    validate_only = false;
end
if nargin < 4
    selected_case_keys = [];
end
if nargin < 5
    environment_overrides = cell(0,2);
end
validateattributes(package_root,{'char','string'},{'scalartext'});
validateattributes(skip_plot,{'logical','numeric'},{'scalar'});
validateattributes(validate_only,{'logical','numeric'},{'scalar'});
assert(iscell(environment_overrides)&&size(environment_overrides,2)==2, ...
    'environment_overrides must be an N-by-2 cell array of environment name/value pairs.');
package_root = char(package_root);
skip_plot = logical(skip_plot);
validate_only = logical(validate_only);

model_dir = fullfile(package_root,'model');
scripts_dir = fullfile(package_root,'scripts');
model_file = fullfile(model_dir,'main_cewgs_restoration.m');
case_file = fullfile(package_root,'data', ...
    'cewgs_33_20_15_test_system.mat');
run_outputs_root = fullfile(package_root,'run_outputs');
run_figures_root = fullfile(package_root,'run_figures');
batch_id = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
output_root = fullfile(run_outputs_root,batch_id);
figure_root = fullfile(run_figures_root,batch_id);

original_matlab_path = path;
path_cleanup = onCleanup(@()path(original_matlab_path));
addpath(model_dir,scripts_dir);
setup_cewgs_dependencies(package_root);

assert(exist(model_file,'file')==2,'Missing model file: %s',model_file);
assert(exist(case_file,'file')==2,'Missing case-data file: %s',case_file);
compressor_data=load(case_file,'GDS_compessor');
assert(isfield(compressor_data,'GDS_compessor')&& ...
    all(abs(compressor_data.GDS_compessor(:,5)-0.32)<1e-12), ...
    'All comparison cases require compressor sigma = 0.32 in column 5.');
assert(exist('sdpvar','file')==2,'YALMIP is not available on the MATLAB path.');
assert(~isempty(which('gurobi')),'The Gurobi MATLAB interface is not available.');
if ~exist(run_outputs_root,'dir'), mkdir(run_outputs_root); end
if ~exist(run_figures_root,'dir'), mkdir(run_figures_root); end
if ~exist(output_root,'dir'), mkdir(output_root); end
if ~exist(figure_root,'dir'), mkdir(figure_root); end

if validate_only
    fprintf('PASS: direct MATLAB optimization entry is configured correctly.\n');
    fprintf('Model: %s\nCase data: %s\n',model_file,case_file);
    return
end

common = {
    'CEWGS_STORAGE_MODEL','dynamic_storage';
    'CEWGS_INITIAL_WATER_PRESSURE_SCALE','0.4';
    'CEWGS_INITIAL_GAS_PRESSURE_SCALE','0.4';
    'CEWGS_GAS_SOURCE_PRESSURE_BAR','20';
    'CEWGS_GAS_NODE1_SOURCE_CAPACITY_KCFH','40';
    'CEWGS_GAS_CUSTOMER_MIN_PRESSURE_BAR','10';
    'CEWGS_GAS_FLOW_SOURCE_NODE','8';
    'CEWGS_GAS_FLOW_SOURCE_RATE_KCFH','10';
    'CEWGS_GAS_FLOW_SOURCE_MODE','upper_bound';
    'CEWGS_GAS_SOC_PENALTY_WEIGHT','0.1';
    'CEWGS_ENFORCE_GAS_NODE_MIN_PRESSURE','1';
    'CEWGS_GAS_CUSTOMER_INACTIVE_PRESSURE_UPPER','0';
    'CEWGS_WATER_SOC_PENALTY_WEIGHT','0.001';
    % Use the pipe-specific capacities from the case-data workbook.  The
    % temporary shared 0.46 m^3/s cap reduced the source-connected capacity
    % from 0.9 to 0.46 m^3/s and is not part of the matched comparison.
    'CEWGS_WATER_PIPE_FLOW_CAPACITY_M3S','';
    'CEWGS_CASE_FILE',case_file;
    'CEWGS_TEST_STEPS','16';
    'CEWGS_EXECUTE_STEPS','16';
    'CEWGS_ROLLING_HORIZON','3';
    'CEWGS_LOOKAHEAD_START_STEP','1';
    'CEWGS_BINARY_MOVE_BLOCKING','0';
    'CEWGS_TOPOLOGY_MOVE_BLOCKING','1';
    'CEWGS_FLOW_DIRECTION_MOVE_BLOCKING','1';
    'CEWGS_WATER_REGION_MOVE_BLOCKING','0';
    'CEWGS_GT_INACTIVE_PRESSURE_UPPER','0';
    'CEWGS_SCA_ROLLING_WARM_START','0';
    'CEWGS_SCA_WARM_START_MODE','none';
    'CEWGS_SCA_WARM_START_WEIGHT','1';
    'CEWGS_FIXED_PRESSURE_FORMULATION','compact_exact';
    'CEWGS_PRESSURE_DRIVEN_FORMULATION','compact_exact';
    'CEWGS_TOPOLOGY_FORMULATION','strengthened';
    'CEWGS_TOPOLOGY_SOLUTION_MODE','direct';
    'CEWGS_SCA_MIP_WARM_START','1';
    'CEWGS_SCA_INTEGER_STABILIZATION','none';
    'CEWGS_THREADS','8';
    'CEWGS_SOLVE_TIME_LIMIT','';
    'CEWGS_ROLLING_WINDOW_TIME_LIMIT','';
    'CEWGS_CAPTURE_SOLVER_TRACE','1';
    'CEWGS_WATER_SOURCE_HEAD','20.6';
    'CEWGS_FIX_WATER_SOURCE_HEAD','1';
    'CEWGS_SWITCH_ACTION_COST','0.1';
    'CEWGS_TIE_SWITCH_MULTIPLIER','2';
    'CEWGS_WT_CAPACITY_BOOST_MW','0';
    'CEWGS_WT_CAPACITY_BOOST_UNIT','1';
    'CEWGS_ELECTRIC_LOAD_REDUCTION_MW','0';
    'CEWGS_DIAGNOSTIC_PLAN_FILE','';
    'CEWGS_DIAGNOSTIC_FIX_SCOPE','topology_and_pumps';
    'CEWGS_FORCE_FULL_WATER_STEP','0';
    'CEWGS_FORCE_FULL_WATER_NODES','';
    'CEWGS_PUMP_POWER_SCA_MODE','full';
    'CEWGS_PUMP_SPEED_MODEL','projection';
    'CEWGS_ENFORCE_PUMP_PHYSICS','1';
    'CEWGS_SCA_PENALTY_START','1';
    'CEWGS_SCA_FLOW_REF_TOL','0.01';
    'CEWGS_SCA_PUMP_REF_TOL','1';
    'CEWGS_SCA_PUMP_TOL','0.1';
    'CEWGS_PUMP_POWER_TOL_KW','0.1';
    'CEWGS_PIPE_REVERSE_SCA','0';
    'CEWGS_RPP_RELAX_TIME_LIMIT','';
    'CEWGS_RPP_PROJECTION_TIME_LIMIT','';
    'CEWGS_RPP_POLISH_TIME_LIMIT','';
    'CEWGS_MCB_MAX_ITERATIONS','';
    'CEWGS_MCB_RELATIVE_GAP','';
    'CEWGS_MCB_RELAX_TIME_LIMIT','';
    'CEWGS_MCB_INTEGER_TIME_LIMIT',''
    };

cases = struct( ...
    'key',{'proposed','fixed_pressure','water_steady','gas_steady'}, ...
    'horizon',{'3','2','3','3'}, ...
    'service',{'pressure_driven','fixed_pressure','pressure_driven','pressure_driven'}, ...
    'water',{'transient','transient','steady_state','transient'}, ...
    'gas',{'transient','transient','transient','steady_state'}, ...
    'water_tank',{'dynamic_storage','dynamic_storage','no_storage','dynamic_storage'}, ...
    'gas_tank',{'dynamic_storage','dynamic_storage','dynamic_storage','no_storage'}, ...
    'water_region',{'0','0','0','0'}, ...
    'mip_gap',{'','0.05','',''});
if ~isempty(selected_case_keys)
    selected_case_keys = string(selected_case_keys);
    keep = ismember(string({cases.key}),selected_case_keys);
    assert(any(keep),'No requested case key matches the configured cases.');
    cases = cases(keep);
end
batch_result_files = cell(1,numel(cases));

environment_names = unique([common(:,1); environment_overrides(:,1); {
    'CEWGS_SERVICE_MODEL';'CEWGS_WATER_NETWORK_MODEL';
    'CEWGS_GAS_NETWORK_MODEL';'CEWGS_WATER_TANK_MODEL';
    'CEWGS_GAS_TANK_MODEL';'CEWGS_MIP_GAP';'CEWGS_OUTPUT_DIR'}]);
old_values = cellfun(@getenv,environment_names,'UniformOutput',false);
environment_cleanup = onCleanup(@()restore_environment(environment_names,old_values));

for case_index = 1:numel(cases)
    apply_environment(common);
    apply_environment(environment_overrides);
    cfg = cases(case_index);
    case_output = fullfile(output_root,cfg.key);
    if ~exist(case_output,'dir'), mkdir(case_output); end
    setenv('CEWGS_SERVICE_MODEL',cfg.service);
    setenv('CEWGS_ROLLING_HORIZON',cfg.horizon);
    horizon_override_index=find(strcmp(environment_overrides(:,1), ...
        'CEWGS_ROLLING_HORIZON'),1,'last');
    if ~isempty(horizon_override_index)
        setenv('CEWGS_ROLLING_HORIZON', ...
            environment_overrides{horizon_override_index,2});
    end
    setenv('CEWGS_WATER_NETWORK_MODEL',cfg.water);
    setenv('CEWGS_GAS_NETWORK_MODEL',cfg.gas);
    setenv('CEWGS_WATER_TANK_MODEL',cfg.water_tank);
    setenv('CEWGS_GAS_TANK_MODEL',cfg.gas_tank);
    setenv('CEWGS_WATER_REGION_MOVE_BLOCKING',cfg.water_region);
    water_region_override_index=find(strcmp(environment_overrides(:,1), ...
        'CEWGS_WATER_REGION_MOVE_BLOCKING'),1,'last');
    if ~isempty(water_region_override_index)
        setenv('CEWGS_WATER_REGION_MOVE_BLOCKING', ...
            environment_overrides{water_region_override_index,2});
    end
    setenv('CEWGS_MIP_GAP',cfg.mip_gap);
    mip_gap_override_index=find(strcmp(environment_overrides(:,1),'CEWGS_MIP_GAP'),1,'last');
    if ~isempty(mip_gap_override_index)
        setenv('CEWGS_MIP_GAP',environment_overrides{mip_gap_override_index,2});
    end
    setenv('CEWGS_OUTPUT_DIR',case_output);

    gap_label = getenv('CEWGS_MIP_GAP');
    if isempty(gap_label), gap_label = 'Gurobi default'; end
    horizon_label=getenv('CEWGS_ROLLING_HORIZON');
    if isempty(horizon_label), horizon_label='model default'; end
    water_region_label=getenv('CEWGS_WATER_REGION_MOVE_BLOCKING');
    if isempty(water_region_label), water_region_label='model default'; end
    fprintf('\n[%d/%d] Running %s: H=%s, water-region linked=%s, MIPGap=%s\n', ...
        case_index,numel(cases),cfg.key,horizon_label,water_region_label,gap_label);
    try
        execute_model_script(model_file);
    catch model_error
        error('Case %s failed: %s',cfg.key,model_error.message);
    end
    listing = dir(fullfile(case_output,'results_*.mat'));
    assert(~isempty(listing),'Case %s did not create a result MAT file.',cfg.key);
    [~,latest_index] = max([listing.datenum]);
    batch_result_files{case_index} = fullfile(listing(latest_index).folder,listing(latest_index).name);
end

if ~skip_plot
    analyze_four_case_results(output_root,case_file,figure_root, ...
        str2double({cases.horizon}));
    analyze_computation_time(output_root,figure_root);
end

batch_manifest = struct( ...
    'batch_id',batch_id, ...
    'created_at',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'result_root',output_root, ...
    'figure_root',figure_root, ...
    'result_files',{batch_result_files}, ...
    'figures_generated',~skip_plot, ...
    'configuration',['Mixed rolling horizons: proposed H=3, fixed-pressure ', ...
        'H=2, steady-water H=3, steady-gas H=3; compressor sigma=0.32 for ', ...
        'all cases; topology and pipe-direction binaries linked across the ', ...
        'prediction steps; other binaries independent; fixed-pressure ', ...
        'MIPGap=0.05']);
if numel(cases)==4
    manifest_name='latest_completed_batch.mat';
else
    manifest_name=['latest_' strjoin({cases.key},'_') '_batch.mat'];
end
save(fullfile(run_outputs_root,manifest_name),'batch_manifest');

fprintf('\nCompleted %d case(s).\nResults: %s\n',numel(cases),output_root);
if ~skip_plot, fprintf('Figures: %s\n',figure_root); end
clear environment_cleanup path_cleanup
end

function apply_environment(pairs)
for pair_index = 1:size(pairs,1)
    setenv(pairs{pair_index,1},pairs{pair_index,2});
end
end

function restore_environment(names,values)
for name_index = 1:numel(names)
    setenv(names{name_index},values{name_index});
end
end

function execute_model_script(model_file)
% Isolate the legacy script's CLEAR command from the user's base workspace.
run(model_file);
end


