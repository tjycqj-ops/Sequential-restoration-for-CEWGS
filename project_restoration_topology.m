function out=project_restoration_topology(line,roots,y_rel,a_rel,y_prev,a_prev,options)
%PROJECT_RESTORATION_TOPOLOGY Project a fractional restoration plan onto a
% rooted radial forest.  This projection contains only topology variables;
% all coupled electric-water-gas decisions are re-optimized in the polish.

started=tic;
if nargin<7,options=struct;end
options=defaults(options);
roots=roots(:);y_prev=round(y_prev(:));a_prev=round(a_prev(:));
if size(line,2)>=3,from=line(:,2);to=line(:,3);else,from=line(:,1);to=line(:,2);end
from=from(:);to=to(:);
[nb,T]=size(y_rel);[ne,Ta]=size(a_rel);
assert(T==Ta&&ne==numel(from),'TopologyRPP:ProjectionSize','Fractional topology dimensions do not match the graph.');
assert(numel(y_prev)==nb&&numel(a_prev)==ne,'TopologyRPP:PreviousSize','Previous topology dimensions do not match the graph.');
assert(all(roots>=1&roots<=nb)&&numel(unique(roots))==numel(roots),'TopologyRPP:Roots','Invalid restoration roots.');

nr=numel(roots);nonroots=setdiff((1:nb)',roots);
ny=nb*T;na=ne*T;nf=ne*T;nR=nr*T;
iy=reshape(1:ny,nb,T);ia=reshape(ny+(1:na),ne,T);
iF=reshape(ny+na+(1:nf),ne,T);iR=reshape(ny+na+nf+(1:nR),nr,T);
N=ny+na+nf+nR;

obj=zeros(N,1);
node_weight=expand_option(options.node_weight,nb,T,'node_weight');
node_benefit=expand_option(options.node_benefit,nb,T,'node_benefit');
obj(iy(:))=node_weight(:).*(1-2*min(max(y_rel(:),0),1))-node_benefit(:);
obj(ia(:))=options.line_weight*(1-2*min(max(a_rel(:),0),1))+options.switch_cost;
lb=-inf(N,1);ub=inf(N,1);vtype=repmat('C',N,1);
lb(iy(:))=0;ub(iy(:))=1;vtype(iy(:))='B';
lb(ia(:))=0;ub(ia(:))=1;vtype(ia(:))='B';
lb(iR(:))=0;

incidence_from=sparse(from,(1:ne)',1,nb,ne);
incidence_to=sparse(to,(1:ne)',1,nb,ne);
source_balance=incidence_from-incidence_to;

A=sparse(0,N);rhs=zeros(0,1);sense=char(zeros(0,1));
    function addrow(indices,values,b,relation)
        row=sparse(1,N);row(indices)=values;
        A(end+1,:)=row;rhs(end+1,1)=b;sense(end+1,1)=relation;
    end

for tau=1:T
    for root=roots'
        addrow(iy(root,tau),1,1,'=');
    end
    if ~isempty(options.force_on_nodes)
        forced=find(options.force_on_nodes(:,min(tau,size(options.force_on_nodes,2))));
        for node=forced',addrow(iy(node,tau),1,1,'=');end
    end
    for ell=1:ne
        addrow([ia(ell,tau),iy(from(ell),tau)],[1,-1],0,'<');
        addrow([ia(ell,tau),iy(to(ell),tau)],[1,-1],0,'<');
        addrow([iF(ell,tau),ia(ell,tau)],[1,-options.commodity_bound],0,'<');
        addrow([iF(ell,tau),ia(ell,tau)],[-1,-options.commodity_bound],0,'<');
    end
    if tau==1
        previous=a_prev;
    else
        previous_indices=ia(:,tau-1);
    end
    if tau==1
        for ell=1:ne,addrow(ia(ell,tau),-1,-previous(ell),'<');end
        addrow(ia(:,tau),ones(1,ne),options.switch_budget+sum(previous),'<');
    else
        for ell=1:ne,addrow([ia(ell,tau),previous_indices(ell)],[-1,1],0,'<');end
        addrow([ia(:,tau);previous_indices],[ones(1,ne),-ones(1,ne)],options.switch_budget,'<');
    end
    addrow([ia(:,tau);iy(:,tau)],[ones(1,ne),-ones(1,nb)],-nr,'=');
    for row_id=1:numel(nonroots)
        node=nonroots(row_id);
        nz=find(source_balance(node,:)~=0);
        addrow([iF(nz,tau);iy(node,tau)],[full(source_balance(node,nz)),1],0,'=');
    end
    for root_id=1:nr
        node=roots(root_id);nz=find(source_balance(node,:)~=0);
        addrow([iF(nz,tau);iR(root_id,tau)],[full(source_balance(node,nz)),-1],0,'=');
    end
end

model=struct('A',A,'rhs',rhs,'sense',sense,'obj',obj,'modelsense','min', ...
    'lb',lb,'ub',ub,'vtype',vtype);
params=struct('OutputFlag',0,'Threads',options.threads,'Seed',options.seed, ...
    'MIPGap',0,'TimeLimit',options.time_limit);
raw=gurobi(model,params);
out=struct('status',raw.status,'y',[],'a',[],'F',[],'root_supply',[], ...
    'objective',Inf,'distance',Inf,'wall_time',toc(started),'raw',raw);
if ~isfield(raw,'x'),return;end
out.y=round(reshape(raw.x(iy(:)),nb,T));
out.a=round(reshape(raw.x(ia(:)),ne,T));
out.F=reshape(raw.x(iF(:)),ne,T);
out.root_supply=reshape(raw.x(iR(:)),nr,T);
out.objective=raw.objval;
out.distance=options.node_weight*sum(abs(out.y-y_rel),'all')+ ...
    options.line_weight*sum(abs(out.a-a_rel),'all');
out.wall_time=toc(started);
end

function options=defaults(options)
defaults=struct('switch_budget',6,'commodity_bound',30,'node_weight',1, ...
    'line_weight',1,'switch_cost',1e-4,'threads',1,'seed',0,'time_limit',10, ...
    'force_on_nodes',[],'node_benefit',0);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(options,names{k}),options.(names{k})=defaults.(names{k});end
end
assert(options.switch_budget>=0&&options.commodity_bound>0&&options.time_limit>0);
end

function value=expand_option(value,nrow,ncol,name)
if isscalar(value),value=repmat(value,nrow,ncol);
elseif isvector(value)&&numel(value)==nrow,value=repmat(value(:),1,ncol);
else,assert(isequal(size(value),[nrow,ncol]),'TopologyRPP:OptionSize','%s has the wrong dimensions.',name);
end
end
