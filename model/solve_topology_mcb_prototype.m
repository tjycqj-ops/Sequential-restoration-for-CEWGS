function [sol,trace,plan]=solve_topology_mcb_prototype( ...
    constraints,objective,ops,y_active,a_active,context,options)
%SOLVE_TOPOLOGY_MCB_PROTOTYPE PROTOTYPE — do not use as production solver.
% Modified combinatorial Benders experiment for one move-blocked rolling
% window.  The graph master contains y/a and theta.  Fixed-topology full
% mixed-integer recourse supplies incumbents; its continuous relaxation
% supplies value-function cuts through topology-fixing-row duals.

started=tic;
if nargin<7,options=struct;end
options=defaults(options);
[original,recovery,diagnostic]=export(constraints,objective,ops);
if ~isempty(diagnostic)&&diagnostic.problem~=0
    error('TopologyMCB:Export','YALMIP export failed: %s',diagnostic.info);
end
original=normalize_model(original);
used=recovery.used_variables(:);
y_cols=variable_columns(y_active,used);a_cols=variable_columns(a_active,used);
assert(all(all(y_cols==y_cols(:,1)))&&all(all(a_cols==a_cols(:,1))), ...
    'TopologyMCB:Scope','Prototype requires move-blocked y/a within the look-ahead window.');
y_model_cols=y_cols(:,1);a_model_cols=a_cols(:,1);
topology_cols=[y_model_cols;a_model_cols];
assert(numel(unique(topology_cols))==numel(topology_cols));

global_lb=options.global_lower_bound;
master=build_master(context,global_lb);
graph_seed=optimize_restoration_topology_graph(context.graph, ...
    struct('threads',options.threads,'time_limit',options.master_time_limit));
if isempty(graph_seed.y)
    [sol,trace,plan]=failure('GRAPH_SEED_FAILED',graph_seed.raw,started);
    return
end
candidate=[graph_seed.y;graph_seed.a];
best_obj=Inf;best_x=[];best_candidate=[];cuts=struct('alpha',{},'beta',{},'kind',{});
rows=struct([]);termination='ITERATION_LIMIT';

for iteration=1:options.max_iterations
    if iteration>1
        next_master=solve_master(master,options);
        if ~isfield(next_master,'x'),termination=['MASTER_',next_master.status];break;end
        candidate=round(next_master.x(master.map.z));
        lb_before=next_master.objval;
    else
        lb_before=global_lb;
    end
    y_candidate=candidate(1:numel(y_model_cols));
    a_candidate=candidate(numel(y_model_cols)+(1:numel(a_model_cols)));

    relax_tic=tic;
    [relax_result,fix_rows]=solve_relaxed_recourse(original,topology_cols,candidate,options);
    relax_seconds=toc(relax_tic);relaxed_value=Inf;dual_cut_added=false;
    if isfield(relax_result,'x')
        relaxed_value=relax_result.objval;
        if isfield(relax_result,'pi')&&numel(relax_result.pi)>=max(fix_rows)
            beta=relax_result.pi(fix_rows);
            alpha=relaxed_value-beta'*candidate;
            [master,cuts]=add_value_cut(master,cuts,alpha,beta,'RELAXED_DUAL');
            dual_cut_added=true;
        else
            point_value=field_or(relax_result,'objbound',relaxed_value);
            [master,cuts]=add_point_cut(master,cuts,candidate,point_value,global_lb,'RELAXED_POINT');
        end
    end

    integer_tic=tic;
    integer_result=solve_integer_recourse(original,topology_cols,candidate,relax_result,options);
    integer_seconds=toc(integer_tic);integer_value=Inf;integer_bound=Inf;
    if isfield(integer_result,'x')
        integer_value=integer_result.objval;
        integer_bound=field_or(integer_result,'objbound',integer_value);
        if integer_value<best_obj
            best_obj=integer_value;best_x=integer_result.x;best_candidate=candidate;
        end
        cut_value=min(integer_value,max(global_lb,integer_bound));
        [master,cuts]=add_point_cut(master,cuts,candidate,cut_value,global_lb,'INTEGER_POINT');
    elseif strcmp(field_or(integer_result,'status','UNKNOWN'),'INFEASIBLE')
        [master,cuts]=add_no_good(master,cuts,candidate,'INTEGER_INFEASIBLE');
    end

    cut_lower_at_candidate=global_lb;
    for cut_id=1:numel(cuts)
        if isfinite(cuts(cut_id).alpha)
            cut_lower_at_candidate=max(cut_lower_at_candidate, ...
                cuts(cut_id).alpha+cuts(cut_id).beta'*candidate);
        end
    end
    if isfinite(integer_value)
        cut_violation=max(cut_lower_at_candidate-integer_value,0);
        assert(cut_violation<=1e-4,'TopologyMCB:InvalidCut', ...
            'A generated value cut exceeds verified integer recourse by %.6g.',cut_violation);
    else
        cut_violation=NaN;
    end

    post=solve_master(master,options);
    if isfield(post,'x'),lb_after=post.objval;else,lb_after=lb_before;end
    if isfinite(best_obj)
        relative_gap=max(best_obj-lb_after,0)/max(1,abs(best_obj));
    else
        relative_gap=Inf;
    end
    row=struct('iteration',iteration,'candidate_y',y_candidate,'candidate_a',a_candidate, ...
        'master_lb_before',lb_before,'relaxed_value',relaxed_value, ...
        'integer_value',integer_value,'integer_bound',integer_bound, ...
        'best_ub',best_obj,'master_lb_after',lb_after,'relative_gap',relative_gap, ...
        'relax_status',field_or(relax_result,'status','UNKNOWN'), ...
        'integer_status',field_or(integer_result,'status','UNKNOWN'), ...
        'relax_seconds',relax_seconds,'integer_seconds',integer_seconds, ...
        'dual_cut_added',dual_cut_added,'cut_count',numel(cuts), ...
        'cut_lower_at_candidate',cut_lower_at_candidate,'cut_violation',cut_violation);
    if isempty(rows),rows=row;else,rows(end+1)=row;end %#ok<AGROW>
    fprintf('[MCB-PROTOTYPE] iter=%d LB=%.6g relax=%.6g integer=%.6g UB=%.6g gap=%.3f%% cuts=%d\n', ...
        iteration,lb_after,relaxed_value,integer_value,best_obj,100*relative_gap,numel(cuts));
    if relative_gap<=options.relative_gap
        termination='GAP';break
    end
    if iteration<options.max_iterations
        next_master=post; %#ok<NASGU> Used at the start of the next iteration.
    end
end

if isempty(best_x)
    [sol,trace,plan]=failure('NO_INTEGER_RECOURSE',struct('status',termination),started);
    trace.iterations=rows;trace.cuts=cuts;return
end
assign(recover(used),best_x(1:numel(used)));
sol=struct('problem',0,'info',sprintf('MCB prototype: %s, UB %.6g',termination,best_obj), ...
    'solvertime',toc(started),'solverinput',struct('model',original), ...
    'solveroutput',struct('best_objective',best_obj,'termination',termination));
plan=struct('y',repmat(best_candidate(1:numel(y_model_cols)),1,size(y_active,2)), ...
    'a',repmat(best_candidate(numel(y_model_cols)+(1:numel(a_model_cols))),1,size(a_active,2)));
trace=struct('prototype',true,'termination',termination,'iterations',rows,'cuts',cuts, ...
    'best_objective',best_obj,'best_y',plan.y,'best_a',plan.a, ...
    'total_seconds',toc(started),'graph_seed',graph_seed, ...
    'topology_variables',numel(topology_cols), ...
    'remaining_integer_variables',sum(original.vtype=='B'|original.vtype=='I')-numel(topology_cols));
end

function master=build_master(c,theta_lb)
line=c.line;roots=c.roots(:);y_prev=round(c.y_previous(:));a_prev=round(c.a_previous(:));
from=line(:,2);to=line(:,3);nb=numel(y_prev);ne=numel(a_prev);nr=numel(roots);
ny=nb;na=ne;nf=ne;nR=nr;iy=1:ny;ia=ny+(1:na);iF=ny+na+(1:nf);
iR=ny+na+nf+(1:nR);itheta=ny+na+nf+nR+1;N=itheta;
lb=-inf(N,1);ub=inf(N,1);vtype=repmat('C',N,1);obj=zeros(N,1);obj(itheta)=1;
lb([iy,ia])=0;ub([iy,ia])=1;vtype([iy,ia])='B';lb(iR)=0;lb(itheta)=theta_lb;
A=sparse(0,N);rhs=[];sense='';
    function addrow(indices,values,b,relation)
        row=sparse(1,N);row(indices)=values;A(end+1,:)=row;
        rhs(end+1,1)=b;sense(end+1,1)=relation;
    end
for root=roots',addrow(iy(root),1,1,'=');end
for node=find(y_prev)',addrow(iy(node),-1,-1,'<');end
for ell=1:ne
    if a_prev(ell)>0,addrow(ia(ell),-1,-1,'<');end
    addrow([ia(ell),iy(from(ell))],[1,-1],0,'<');
    addrow([ia(ell),iy(to(ell))],[1,-1],0,'<');
    addrow([iF(ell),ia(ell)],[1,-c.commodity_bound],0,'<');
    addrow([iF(ell),ia(ell)],[-1,-c.commodity_bound],0,'<');
end
addrow(ia,ones(1,ne),c.switch_budget+sum(a_prev),'<');
addrow([ia,iy],[ones(1,ne),-ones(1,nb)],-nr,'=');
inc_from=sparse(from,(1:ne)',1,nb,ne);inc_to=sparse(to,(1:ne)',1,nb,ne);balance=inc_from-inc_to;
for node=setdiff((1:nb)',roots)'
    nz=find(balance(node,:)~=0);addrow([iF(nz),iy(node)],[full(balance(node,nz)),1],0,'=');
end
for root_id=1:nr
    node=roots(root_id);nz=find(balance(node,:)~=0);
    addrow([iF(nz),iR(root_id)],[full(balance(node,nz)),-1],0,'=');
end
master=struct('A',A,'rhs',rhs,'sense',sense,'obj',obj,'modelsense','min', ...
    'lb',lb,'ub',ub,'vtype',vtype,'map',struct('y',iy,'a',ia,'z',[iy,ia],'theta',itheta));
end

function result=solve_master(master,options)
model=rmfield(master,'map');
params=struct('OutputFlag',0,'Threads',options.threads,'Seed',options.seed, ...
    'MIPGap',0,'TimeLimit',options.master_time_limit);
result=gurobi(model,params);
end

function [result,fix_rows]=solve_relaxed_recourse(original,cols,values,options)
model=original;model.vtype(:)='C';
[model,fix_rows]=append_fix_rows(model,cols,values);
params=struct('OutputFlag',0,'Threads',options.threads,'Seed',options.seed, ...
    'TimeLimit',options.relax_time_limit,'NonConvex',0,'QCPDual',1);
result=gurobi(model,params);
end

function result=solve_integer_recourse(original,cols,values,relax_result,options)
model=original;model.lb(cols)=values;model.ub(cols)=values;model.vtype(cols)='C';
if isfield(relax_result,'x'),model.start=relax_result.x;end
params=struct('OutputFlag',0,'Threads',options.threads,'Seed',options.seed, ...
    'TimeLimit',options.integer_time_limit,'MIPGap',options.mip_gap, ...
    'NonConvex',0,'DualReductions',0);
result=gurobi(model,params);
end

function [model,rows]=append_fix_rows(model,cols,values)
nvar=numel(model.obj);nfix=numel(cols);base=size(model.A,1);
extra=sparse(1:nfix,cols,1,nfix,nvar);
model.A=[model.A;extra];model.rhs=[model.rhs;values(:)];
model.sense=[model.sense;repmat('=',nfix,1)];rows=base+(1:nfix);
end

function [master,cuts]=add_value_cut(master,cuts,alpha,beta,kind)
row=sparse(1,numel(master.obj));row(master.map.z)=beta(:)';row(master.map.theta)=-1;
master.A(end+1,:)=row;master.rhs(end+1,1)=-alpha;master.sense(end+1,1)='<';
cuts(end+1)=struct('alpha',alpha,'beta',beta(:),'kind',kind);
end

function [master,cuts]=add_point_cut(master,cuts,z,value,lower_bound,kind)
if ~isfinite(value),return;end
scale=max(value-lower_bound,0);beta=scale*(2*z(:)-1);
alpha=value-scale*sum(z);
[master,cuts]=add_value_cut(master,cuts,alpha,beta,kind);
end

function [master,cuts]=add_no_good(master,cuts,z,kind)
beta=2*z(:)-1;row=sparse(1,numel(master.obj));row(master.map.z)=beta';
master.A(end+1,:)=row;master.rhs(end+1,1)=sum(z)-1;master.sense(end+1,1)='<';
cuts(end+1)=struct('alpha',NaN,'beta',beta,'kind',kind);
end

function cols=variable_columns(variable,used)
cols=zeros(numel(variable),1);
for k=1:numel(variable)
    id=getvariables(variable(k));assert(isscalar(id),'TopologyMCB:VariableMap');
    [found,cols(k)]=ismember(id,used);assert(found,'TopologyMCB:VariableMap');
end
cols=reshape(cols,size(variable));
end

function model=normalize_model(model)
model.obj=model.obj(:);model.lb=model.lb(:);model.ub=model.ub(:);
model.vtype=model.vtype(:);model.rhs=model.rhs(:);model.sense=model.sense(:);
if ~isfield(model,'objcon'),model.objcon=0;end
end

function value=field_or(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end

function [sol,trace,plan]=failure(reason,result,started)
status=field_or(result,'status','UNKNOWN');
sol=struct('problem',1,'info',[reason,': ',status],'solvertime',toc(started), ...
    'solverinput',struct,'solveroutput',struct('result',result));
trace=struct('prototype',true,'termination',reason,'iterations',struct([]), ...
    'cuts',struct([]),'total_seconds',toc(started));plan=struct('y',[],'a',[]);
end

function options=defaults(options)
defaults=struct('max_iterations',3,'relative_gap',0.05,'global_lower_bound',0, ...
    'master_time_limit',10,'relax_time_limit',30,'integer_time_limit',30, ...
    'mip_gap',0.05,'threads',1,'seed',0);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(options,names{k}),options.(names{k})=defaults.(names{k});end
end
end
