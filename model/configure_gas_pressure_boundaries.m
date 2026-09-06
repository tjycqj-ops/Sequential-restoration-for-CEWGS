function [GDS_node,source_nodes,pressure_source_nodes,uniform_min_nodes] = ...
    configure_gas_pressure_boundaries(GDS_node,GDS_station,flow_source_node,uniform_min_bar)
%CONFIGURE_GAS_PRESSURE_BOUNDARIES Separate pressure and flow source bounds.
% Fixed-pressure source nodes keep their data-table pressure limits. Ordinary
% nodes and a flow-controlled source receive the common minimum-pressure
% override because the flow source pressure is determined by network physics.

validateattributes(GDS_node,{'numeric'},{'2d','nonempty','finite'});
validateattributes(GDS_station,{'numeric'},{'2d','nonempty','finite'});
validateattributes(flow_source_node,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(uniform_min_bar,{'numeric'},{'scalar','positive','finite'});
assert(size(GDS_node,2)>=4,'GDS_node must include node, minimum, and maximum pressure columns.');
assert(size(GDS_station,2)>=2,'GDS_station must include its connected gas-node column.');

source_nodes=GDS_station(:,2);
if flow_source_node>0
    assert(nnz(source_nodes==flow_source_node)==1, ...
        'The flow-controlled source must identify exactly one gas station node.');
end

pressure_source_nodes=setdiff(source_nodes,flow_source_node,'stable');
uniform_min_nodes=setdiff(GDS_node(:,1),pressure_source_nodes,'stable');
uniform_min_rows=ismember(GDS_node(:,1),uniform_min_nodes);
GDS_node(uniform_min_rows,3)=uniform_min_bar;

assert(all(GDS_node(:,4)>=GDS_node(:,3)), ...
    'Configured gas minimum pressure exceeds a node maximum pressure.');
end
