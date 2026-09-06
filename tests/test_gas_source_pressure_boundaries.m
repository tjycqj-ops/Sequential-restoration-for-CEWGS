function test_gas_source_pressure_boundaries
%TEST_GAS_SOURCE_PRESSURE_BOUNDARIES Regression for source-bound conflict.

package_root=fileparts(fileparts(mfilename('fullpath')));
model_dir=fullfile(package_root,'model');
case_file=fullfile(package_root,'data','cewgs_33_20_15_test_system.mat');
addpath(model_dir);
S=load(case_file,'GDS_node','GDS_station');

[nodes,source_nodes,pressure_sources,uniform_nodes]= ...
    configure_gas_pressure_boundaries(S.GDS_node,S.GDS_station,8,10);
row1=find(nodes(:,1)==1);
row8=find(nodes(:,1)==8);
assert(nodes(row8,3)==10, ...
    'Flow source node 8 must use the configured 10 bar minimum pressure.');
assert(nodes(row1,3)==S.GDS_node(row1,3), ...
    'Pressure source node 1 must preserve its data-table minimum pressure.');
assert(ismember(8,uniform_nodes)&&~ismember(8,pressure_sources));
assert(ismember(1,source_nodes)&&ismember(1,pressure_sources));

[nodes_without_flow_source,~,pressure_sources_without_flow_source,~]= ...
    configure_gas_pressure_boundaries(S.GDS_node,S.GDS_station,0,10);
assert(nodes_without_flow_source(row8,3)==S.GDS_node(row8,3), ...
    'Without a flow-source designation, node 8 must preserve its source boundary.');
assert(all(ismember(source_nodes,pressure_sources_without_flow_source)));

fprintf('PASS: gas flow-source and pressure-source boundaries are separated.\n');
end

