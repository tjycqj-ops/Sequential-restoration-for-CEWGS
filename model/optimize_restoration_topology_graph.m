function out=optimize_restoration_topology_graph(context,options)
%OPTIMIZE_RESTORATION_TOPOLOGY_GRAPH Direct rooted-forest topology MILP.
% The graph layer maximizes a service-potential surrogate subject to the
% actual switching budget.  It does not consume a fractional topology and
% therefore is not a rounding/projection method.  The returned topology is
% intended to be fixed and physically checked by the full CEWGS dispatch.

started=tic;
if nargin<2,options=struct;end
options=defaults(options);
required={'line','roots','y_previous','a_previous','customer_nodes','priority', ...
    'pump_buses','water_pump_paths','compressor_buses','gas_demand', ...
    'gas_base_capacity','compressor_capacity_gain','switch_budget'};
for k=1:numel(required)
    assert(isfield(context,required{k}),'TopologyGraph:MissingField', ...
        'Missing context field %s.',required{k});
end

line=context.line;roots=context.roots(:);
if size(line,2)>=3,from=line(:,2);to=line(:,3);else,from=line(:,1);to=line(:,2);end
from=from(:);to=to(:);ne=numel(from);
y_prev=round(context.y_previous(:));a_prev=round(context.a_previous(:));
nb=numel(y_prev);nr=numel(roots);
customer_nodes=context.customer_nodes(:);priority=context.priority(:);
pump_buses=context.pump_buses(:);compressor_buses=context.compressor_buses(:);
water_paths=context.water_pump_paths;
gas_demand=context.gas_demand(:);compressor_gain=context.compressor_capacity_gain(:);
np=numel(pump_buses);nw=size(water_paths,1);nc=numel(compressor_buses);ng=numel(gas_demand);
assert(numel(a_prev)==ne&&numel(priority)==numel(customer_nodes));
assert(size(water_paths,2)==np&&nw==numel(priority));
assert(ng==numel(priority)&&numel(compressor_gain)==nc);
assert(all(roots>=1&roots<=nb)&&all(customer_nodes>=1&customer_nodes<=nb));
assert(all(pump_buses>=1&pump_buses<=nb)&&all(compressor_buses>=1&compressor_buses<=nb));

cursor=0;
iy=cursor+(1:nb);cursor=cursor+nb;
ia=cursor+(1:ne);cursor=cursor+ne;
iF=cursor+(1:ne);cursor=cursor+ne;
iR=cursor+(1:nr);cursor=cursor+nr;
ie=cursor+(1:numel(priority));cursor=cursor+numel(priority);
ixp=cursor+(1:np);cursor=cursor+np;
iw=cursor+(1:nw);cursor=cursor+nw;
ixc=cursor+(1:nc);cursor=cursor+nc;
ig=cursor+(1:ng);cursor=cursor+ng;N=cursor;

obj=zeros(N,1);
obj(ia)=options.switch_cost;
obj(ie)=-priority;
obj(iw)=-priority;
obj(ig)=-priority;
lb=-inf(N,1);ub=inf(N,1);vtype=repmat('C',N,1);
bounded=[iy,ia,ie,ixp,iw,ixc,ig];lb(bounded)=0;ub(bounded)=1;
vtype([iy,ia])='B';lb(iR)=0;

A=sparse(0,N);rhs=zeros(0,1);sense=char(zeros(0,1));
    function addrow(indices,values,b,relation)
        row=sparse(1,N);row(indices)=values;
        A(end+1,:)=row;rhs(end+1,1)=b;sense(end+1,1)=relation;
    end

for root=roots',addrow(iy(root),1,1,'=');end
for node=find(y_prev)',addrow(iy(node),-1,-1,'<');end
for ell=1:ne
    if a_prev(ell)>0,addrow(ia(ell),-1,-1,'<');end
    addrow([ia(ell),iy(from(ell))],[1,-1],0,'<');
    addrow([ia(ell),iy(to(ell))],[1,-1],0,'<');
    addrow([iF(ell),ia(ell)],[1,-options.commodity_bound],0,'<');
    addrow([iF(ell),ia(ell)],[-1,-options.commodity_bound],0,'<');
end
addrow(ia,ones(1,ne),context.switch_budget+sum(a_prev),'<');
addrow([ia,iy],[ones(1,ne),-ones(1,nb)],-nr,'=');

incidence_from=sparse(from,(1:ne)',1,nb,ne);
incidence_to=sparse(to,(1:ne)',1,nb,ne);
source_balance=incidence_from-incidence_to;
nonroots=setdiff((1:nb)',roots);
for node=nonroots'
    nz=find(source_balance(node,:)~=0);
    addrow([iF(nz),iy(node)],[full(source_balance(node,nz)),1],0,'=');
end
for root_id=1:nr
    node=roots(root_id);nz=find(source_balance(node,:)~=0);
    addrow([iF(nz),iR(root_id)],[full(source_balance(node,nz)),-1],0,'=');
end

for k=1:numel(priority)
    addrow([ie(k),iy(customer_nodes(k))],[1,-1],0,'<');
end
for pump=1:np
    addrow([ixp(pump),iy(pump_buses(pump))],[1,-1],0,'<');
end
for customer=1:nw
    for pump=find(water_paths(customer,:)>0.5)
        addrow([iw(customer),ixp(pump)],[1,-1],0,'<');
    end
end
for compressor=1:nc
    addrow([ixc(compressor),iy(compressor_buses(compressor))],[1,-1],0,'<');
end
addrow([ig,ixc],[gas_demand',-compressor_gain'],context.gas_base_capacity,'<');

model=struct('A',A,'rhs',rhs,'sense',sense,'obj',obj,'modelsense','min', ...
    'lb',lb,'ub',ub,'vtype',vtype);
params=struct('OutputFlag',options.output_flag,'Threads',options.threads, ...
    'Seed',options.seed,'MIPGap',0,'TimeLimit',options.time_limit);
raw=gurobi(model,params);
out=struct('status',raw.status,'y',[],'a',[],'F',[],'root_supply',[], ...
    'electric_service',[],'pump_available',[],'water_service',[], ...
    'compressor_available',[],'gas_service',[],'objective',Inf, ...
    'graph_residual',Inf,'wall_time',toc(started),'raw',raw);
if ~isfield(raw,'x'),return;end
out.y=round(raw.x(iy));out.a=round(raw.x(ia));out.F=raw.x(iF);
out.root_supply=raw.x(iR);out.electric_service=raw.x(ie);
out.pump_available=raw.x(ixp);out.water_service=raw.x(iw);
out.compressor_available=raw.x(ixc);out.gas_service=raw.x(ig);
out.objective=raw.objval;
endpoint_res=max([out.a-out.y(from);out.a-out.y(to);0]);
flow_res=max([abs(source_balance(nonroots,:)*out.F+out.y(nonroots));0]);
root_res=max([abs(source_balance(roots,:)*out.F-out.root_supply);0]);
tree_res=abs(sum(out.a)-sum(out.y)+nr);
switch_res=max(sum(out.a-a_prev)-context.switch_budget,0);
monotone_res=max([y_prev-out.y;a_prev-out.a;0]);
out.graph_residual=max([endpoint_res;flow_res;root_res;tree_res;switch_res;monotone_res]);
out.wall_time=toc(started);
end

function options=defaults(options)
defaults=struct('switch_cost',1e-3,'commodity_bound',165, ...
    'threads',1,'seed',0,'time_limit',10,'output_flag',0);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(options,names{k}),options.(names{k})=defaults.(names{k});end
end
assert(options.commodity_bound>0&&options.time_limit>0&&options.switch_cost>=0);
end
