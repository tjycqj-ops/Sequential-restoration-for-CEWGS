function paths=restoration_water_pump_paths(pipe,pumps,customers)
%RESTORATION_WATER_PUMP_PATHS Pump incidence on source-customer paths.
paths=zeros(numel(customers),size(pumps,1));
G=graph(pipe(:,2),pipe(:,3),[],max(pipe(:,2:3),[],'all'));
for customer=1:numel(customers)
    nodes=shortestpath(G,1,customers(customer));
    for step=1:numel(nodes)-1
        u=nodes(step);v=nodes(step+1);
        edge=find((pipe(:,2)==u&pipe(:,3)==v)|(pipe(:,2)==v&pipe(:,3)==u),1);
        pump=find(pumps(:,1)==edge,1);
        if ~isempty(pump),paths(customer,pump)=1;end
    end
end
end
