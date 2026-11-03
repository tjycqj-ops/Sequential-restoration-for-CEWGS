function [common,cases] = four_case_configuration(package_root,overrides)
%FOUR_CASE_CONFIGURATION Shared numerical settings and model-only case differences.
if nargin<2,overrides=cell(0,2);end
assert(iscell(overrides)&&size(overrides,2)==2, ...
    'Overrides must be an N-by-2 cell array.');
case_file=fullfile(package_root,'data','cewgs_33_20_15_test_system.mat');
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
    % Use the pipe-specific capacities in the packaged system data.
    'CEWGS_WATER_PIPE_FLOW_CAPACITY_M3S','';
    'CEWGS_CASE_FILE',case_file;
    'CEWGS_WIND_FILE',fullfile(package_root,'data','wind_uncertainty.txt');
    'CEWGS_TEST_STEPS','16';
    'CEWGS_EXECUTE_STEPS','16';
    'CEWGS_ROLLING_HORIZON','3';
    'CEWGS_MIP_GAP','0.0001';
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
    'CEWGS_SCA_ITER_MAX','40';
    'CEWGS_SCA_FLOW_REF_TOL','0.01';
    'CEWGS_SCA_PUMP_REF_TOL','1';
    'CEWGS_SCA_PUMP_TOL','0.1';
    'CEWGS_PUMP_POWER_TOL_KW','0.1';
    'CEWGS_PIPE_REVERSE_SCA','0';
    };

cases = struct( ...
    'key',{'proposed','fixed_pressure','water_steady','gas_steady'}, ...
    'service',{'pressure_driven','fixed_pressure','pressure_driven','pressure_driven'}, ...
    'water',{'transient','transient','steady_state','transient'}, ...
    'gas',{'transient','transient','transient','steady_state'}, ...
    'water_tank',{'dynamic_storage','dynamic_storage','no_storage','dynamic_storage'}, ...
    'gas_tank',{'dynamic_storage','dynamic_storage','dynamic_storage','no_storage'}, ...
    'label',{'Proposed','Fixed-pressure','Steady-state water','Steady-state gas'});

for k=1:size(overrides,1)
    idx=find(strcmp(common(:,1),overrides{k,1}));
    assert(isscalar(idx),'Unknown shared setting: %s',overrides{k,1});
    validateattributes(overrides{k,2},{'char','string'},{'scalartext'});
    common{idx,2}=char(overrides{k,2});
end
steps=str2double(common{strcmp(common(:,1),'CEWGS_TEST_STEPS'),2});
horizon=str2double(common{strcmp(common(:,1),'CEWGS_ROLLING_HORIZON'),2});
executed=str2double(common{strcmp(common(:,1),'CEWGS_EXECUTE_STEPS'),2});
gap=str2double(common{strcmp(common(:,1),'CEWGS_MIP_GAP'),2});
assert(isfinite(steps)&&steps>=1&&steps<=32&&steps==floor(steps), ...
    'Total steps must be an integer from 1 to 32.');
assert(isfinite(horizon)&&horizon>=1&&horizon<=steps&&horizon==floor(horizon), ...
    'Rolling horizon must be an integer from 1 to total steps.');
assert(isfinite(executed)&&executed>=1&&executed<=steps&&executed==floor(executed), ...
    'Executed steps must be an integer from 1 to total steps.');
assert(isfinite(gap)&&gap>=0&&gap<1,'MIPGap must be a finite number in [0,1).');
end
