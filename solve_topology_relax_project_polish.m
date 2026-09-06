function [sol,trace,plan]=solve_topology_relax_project_polish( ...
    constraints,objective,ops,y_active,a_active,context,options)
%SOLVE_TOPOLOGY_RELAX_PROJECT_POLISH Relax only y/a, graph-project them,
% then fix them and polish every other integer in the original coupled model.

started=tic;
if nargin<7,options=struct;end
options=defaults(options);
[original,recovery,diagnostic]=export(constraints,objective,ops);
if ~isempty(diagnostic)&&diagnostic.problem~=0
    error('TopologyRPP:Export','YALMIP export failed: %s',diagnostic.info);
end
original=normalize_model(original);
used=recovery.used_variables(:);
y_cols=variable_columns(y_active,used);
a_cols=variable_columns(a_active,used);
topology_cols=unique([y_cols(:);a_cols(:)]);
assert(all(original.vtype(topology_cols)=='B'|original.vtype(topology_cols)=='I'), ...
    'TopologyRPP:VariableType','The selected node/line variables are not binary in the original model.');

non_topology=setdiff((1:numel(original.obj))',topology_cols);
non_topology_vtype=original.vtype(non_topology); % Contract: untouched discrete recourse.
relaxed=original;
relaxed.vtype(topology_cols)='C';
relaxed.lb(topology_cols)=max(relaxed.lb(topology_cols),0);
relaxed.ub(topology_cols)=min(relaxed.ub(topology_cols),1);
assert(isequal(relaxed.vtype(non_topology),non_topology_vtype), ...
    'TopologyRPP:UnexpectedRelaxation','A non-topology variable was relaxed.');

relax_params=solver_params(options.relax_time_limit,options);
relax_tic=tic;relax_result=gurobi(relaxed,relax_params);relax_seconds=toc(relax_tic);
if ~isfield(relax_result,'x')
    [sol,trace,plan]=failure('TOPOLOGY_RELAXATION_FAILED',relax_result,started,relax_seconds);
    return
end

y_rel=reshape(relax_result.x(y_cols),size(y_active));
a_rel=reshape(relax_result.x(a_cols),size(a_active));
move_blocked=(numel(unique(y_cols))<numel(y_cols))||(numel(unique(a_cols))<numel(a_cols));
if move_blocked
    % The rolling model represents every look-ahead column with the same
    % YALMIP binary.  Project that unique decision once, then replicate it.
    y_projection_reference=mean(y_rel,2);
    a_projection_reference=mean(a_rel,2);
else
    y_projection_reference=y_rel;
    a_projection_reference=a_rel;
end
force_on_nodes=false(size(y_rel));
if isfield(context,'protected')
    for group_id=1:numel(context.protected)
        protected_cols=variable_columns(context.protected(group_id).variable,used);
        protected_values=reshape(relax_result.x(protected_cols),size(context.protected(group_id).variable));
        protected_buses=context.protected(group_id).buses(:);
        assert(size(protected_values,1)==numel(protected_buses),'TopologyRPP:ProtectedSize');
        for tau=1:size(protected_values,2)
            force_on_nodes(protected_buses(protected_values(:,tau)>0.5),tau)=true;
        end
    end
end
if move_blocked,force_on_nodes=any(force_on_nodes,2);end
projection_options=struct('switch_budget',context.switch_budget, ...
    'commodity_bound',context.commodity_bound,'threads',options.threads, ...
    'seed',options.seed,'time_limit',options.projection_time_limit, ...
    'node_weight',options.node_weight,'line_weight',options.line_weight, ...
    'switch_cost',options.projection_switch_cost,'force_on_nodes',force_on_nodes);
if isfield(context,'node_benefit'),projection_options.node_benefit=context.node_benefit;end
projection=project_restoration_topology(context.line,context.roots,y_projection_reference,a_projection_reference, ...
    context.y_previous,context.a_previous,projection_options);
if isempty(projection.y)
    [sol,trace,plan]=failure('TOPOLOGY_PROJECTION_FAILED',projection.raw,started,relax_seconds);
    trace.projection=projection;
    return
end
if move_blocked&&size(y_active,2)>1
    projection.y=repmat(projection.y,1,size(y_active,2));
    projection.a=repmat(projection.a,1,size(a_active,2));
    projection.F=repmat(projection.F,1,size(a_active,2));
    projection.root_supply=repmat(projection.root_supply,1,size(a_active,2));
end

projected_by_col=nan(numel(original.obj),1);
projected_by_col=place_values(projected_by_col,y_cols,projection.y);
projected_by_col=place_values(projected_by_col,a_cols,projection.a);
projected_values=projected_by_col(topology_cols);
assert(all(isfinite(projected_values))&&all(ismember(projected_values,[0,1])), ...
    'TopologyRPP:ProjectionValues','Projected topology is not binary.');

fixed=original;
fixed.lb(topology_cols)=projected_values;
fixed.ub(topology_cols)=projected_values;
fixed.vtype(topology_cols)='C';
fixed.start=relax_result.x;
fixed.start(topology_cols)=projected_values;
polish_params=solver_params(options.polish_time_limit,options);
polish_tic=tic;polish_result=gurobi(fixed,polish_params);polish_seconds=toc(polish_tic);
if ~isfield(polish_result,'x')
    sol=struct('problem',1,'info',['RPP polish failed: ',polish_result.status], ...
        'solvertime',toc(started),'solverinput',struct('model',fixed), ...
        'solveroutput',struct('result',polish_result));
else
    assert(numel(used)<=numel(polish_result.x),'TopologyRPP:RecoverySize','Solver vector is shorter than the YALMIP recovery map.');
    assign(recover(used),polish_result.x(1:numel(used)));
    sol=struct('problem',0,'info',sprintf('RPP relax=%s, project=%s, polish=%s', ...
        relax_result.status,projection.status,polish_result.status), ...
        'solvertime',toc(started),'solverinput',struct('model',fixed), ...
        'solveroutput',struct('result',polish_result,'relax_result',relax_result, ...
        'projection',projection));
end
plan=struct('y',projection.y,'a',projection.a,'y_relaxed',y_rel,'a_relaxed',a_rel);
trace=struct('relax_status',relax_result.status,'polish_status',polish_result.status, ...
    'relax_seconds',relax_seconds,'projection_seconds',projection.wall_time, ...
    'polish_seconds',polish_seconds,'total_seconds',toc(started), ...
    'relaxed_objective',field_or(relax_result,'objval',NaN), ...
    'polished_objective',field_or(polish_result,'objval',NaN), ...
    'projection_distance',projection.distance,'topology_variables',numel(topology_cols), ...
    'remaining_integer_variables',sum(original.vtype(non_topology)=='B'|original.vtype(non_topology)=='I'), ...
    'move_blocked',move_blocked,'forced_on_nodes',sum(force_on_nodes,'all'), ...
    'y_relaxed',y_rel,'a_relaxed',a_rel,'projected_y',projection.y, ...
    'projected_a',projection.a,'projection',projection);
end

function cols=variable_columns(variable,used)
cols=zeros(numel(variable),1);
for k=1:numel(variable)
    id=getvariables(variable(k));
    assert(isscalar(id),'TopologyRPP:VariableMap','A topology element does not map to one YALMIP variable.');
    [found,cols(k)]=ismember(id,used);
    assert(found,'TopologyRPP:VariableMap','A topology variable is absent from the solver export.');
end
cols=reshape(cols,size(variable));
end

function target=place_values(target,cols,values)
for k=1:numel(cols)
    col=cols(k);value=round(values(k));
    if isfinite(target(col))
        assert(target(col)==value,'TopologyRPP:MoveBlocking','Repeated topology variable received inconsistent projected values.');
    else
        target(col)=value;
    end
end
end

function model=normalize_model(model)
model.obj=model.obj(:);model.lb=model.lb(:);model.ub=model.ub(:);
model.vtype=model.vtype(:);model.rhs=model.rhs(:);model.sense=model.sense(:);
if ~isfield(model,'objcon'),model.objcon=0;end
end

function params=solver_params(time_limit,options)
params=struct('OutputFlag',options.output_flag,'Threads',options.threads, ...
    'Seed',options.seed,'MIPGap',options.mip_gap,'TimeLimit',time_limit, ...
    'NonConvex',0);
end

function value=field_or(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end

function [sol,trace,plan]=failure(reason,result,started,relax_seconds)
status=field_or(result,'status','UNKNOWN');
sol=struct('problem',1,'info',[reason,': ',status],'solvertime',toc(started), ...
    'solverinput',struct,'solveroutput',struct('result',result));
trace=struct('relax_status',status,'polish_status','NOT_RUN', ...
    'relax_seconds',relax_seconds,'projection_seconds',0,'polish_seconds',0, ...
    'total_seconds',toc(started),'relaxed_objective',NaN,'polished_objective',NaN, ...
    'projection_distance',NaN,'topology_variables',0,'remaining_integer_variables',0);
plan=struct('y',[],'a',[],'y_relaxed',[],'a_relaxed',[]);
end

function options=defaults(options)
defaults=struct('relax_time_limit',60,'projection_time_limit',10, ...
    'polish_time_limit',60,'mip_gap',0.05,'threads',1,'seed',0, ...
    'output_flag',0,'node_weight',1,'line_weight',1,'projection_switch_cost',1e-4, ...
    'integral_tolerance',1e-5);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(options,names{k}),options.(names{k})=defaults.(names{k});end
end
end
