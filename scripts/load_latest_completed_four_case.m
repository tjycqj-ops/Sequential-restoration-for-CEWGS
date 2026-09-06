function latest=load_latest_completed_four_case(package_root)
%LOAD_LATEST_COMPLETED_FOUR_CASE Resolve the newest complete four-case batch.
% The manifest is written only after all four optimization cases finish, so
% plotting scripts never consume a partially written in-progress batch.

if nargin<1||isempty(package_root)
    package_root=fileparts(fileparts(mfilename('fullpath')));
end
manifest_file=fullfile(package_root,'run_outputs','latest_completed_batch.mat');
assert(exist(manifest_file,'file')==2, ...
    'No completed four-case batch manifest was found: %s',manifest_file);

M=load(manifest_file,'batch_manifest');
assert(isfield(M,'batch_manifest'),'Invalid completed-batch manifest: %s',manifest_file);
batch_manifest=M.batch_manifest;
required_manifest_fields={'batch_id','result_root','figure_root','result_files'};
assert(all(isfield(batch_manifest,required_manifest_fields)), ...
    'Completed-batch manifest fields are incomplete: %s',manifest_file);

result_root=fullfile(package_root,'run_outputs',char(batch_manifest.batch_id));
figure_root=fullfile(package_root,'run_figures',char(batch_manifest.batch_id));
raw_files=batch_manifest.result_files;
assert(iscell(raw_files)&&numel(raw_files)==4, ...
    'The latest completed batch must contain exactly four result files.');
case_keys={'proposed','fixed_pressure','water_steady','gas_steady'};
case_labels={'Proposed','Fixed-pressure','Steady-state water','Steady-state gas'};
ordered_files=cell(1,4);
case_files=strings(1,4);
horizons=zeros(1,4);
step_hours=zeros(1,4);

for raw_index=1:numel(raw_files)
    result_file=char(raw_files{raw_index});
    if exist(result_file,'file')~=2
        path_parts=strsplit(strrep(result_file,'\','/'),'/');
        result_file=fullfile(result_root,path_parts{end-1},path_parts{end});
    end
    assert(exist(result_file,'file')==2,'Missing result file: %s',result_file);
    info=load(result_file,'service_model','water_network_model', ...
        'gas_network_model','case_file','Th','TT','executed_steps','detT', ...
        'solver_problem_r');
    required_result_fields={'service_model','water_network_model', ...
        'gas_network_model','case_file','Th','TT','executed_steps','detT', ...
        'solver_problem_r'};
    assert(all(isfield(info,required_result_fields)), ...
        'Result metadata fields are incomplete: %s',result_file);
    assert(info.executed_steps==info.TT&&all(info.solver_problem_r==0), ...
        'Result is not a completed successful run: %s',result_file);

    if strcmp(info.service_model,'fixed_pressure')
        slot=2;
    elseif strcmp(info.water_network_model,'steady_state')&& ...
            strcmp(info.gas_network_model,'transient')
        slot=3;
    elseif strcmp(info.water_network_model,'transient')&& ...
            strcmp(info.gas_network_model,'steady_state')
        slot=4;
    elseif strcmp(info.service_model,'pressure_driven')&& ...
            strcmp(info.water_network_model,'transient')&& ...
            strcmp(info.gas_network_model,'transient')
        slot=1;
    else
        error('Unrecognized comparison case in result file: %s',result_file);
    end
    assert(isempty(ordered_files{slot}), ...
        'Duplicate %s result in completed batch.',case_keys{slot});
    ordered_files{slot}=result_file;
    case_files(slot)=string(info.case_file);
    horizons(slot)=info.Th;
    step_hours(slot)=info.detT;
end

assert(all(~cellfun(@isempty,ordered_files)), ...
    'The latest completed batch does not contain all four comparison cases.');
assert(all(case_files==case_files(1)), ...
    'The latest four results do not use the same case-data file.');
assert(all(abs(step_hours-step_hours(1))<1e-12), ...
    'The latest four results do not use the same execution time step.');
case_file=char(case_files(1));
if exist(case_file,'file')~=2
    [~,case_name,case_extension]=fileparts(strrep(case_file,'\','/'));
    case_file=fullfile(package_root,'input',[case_name case_extension]);
end
assert(exist(case_file,'file')==2,'Missing case-data file: %s',case_file);

latest=struct('batch_id',batch_manifest.batch_id, ...
    'manifest_file',manifest_file,'result_root',result_root, ...
    'figure_root',figure_root,'result_files',{ordered_files}, ...
    'case_file',case_file,'case_keys',{case_keys},'case_labels',{case_labels}, ...
    'horizons',horizons,'step_hours',step_hours(1));
end
