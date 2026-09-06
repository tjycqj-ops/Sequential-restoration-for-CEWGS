function cycle_matrix = restoration_cycle_matrix(edge_nodes)
%RESTORATION_CYCLE_MATRIX Fundamental-cycle incidence matrix for an undirected graph.
% Each row identifies one cycle. It is used only to add redundant radiality
% cuts, sum(a(cycle)) <= |cycle|-1, and introduces no new decision variables.

assert(size(edge_nodes,2)==2,'edge_nodes must be an m-by-2 endpoint matrix.');
assert(all(edge_nodes(:)>=1)&&all(edge_nodes(:)==round(edge_nodes(:))), ...
    'Graph endpoints must be positive integer node indices.');

edge_nodes=round(edge_nodes);
edge_count=size(edge_nodes,1);
node_count=max(edge_nodes,[],'all');
parent=1:node_count;
tree_edge=false(edge_count,1);
chord_rows=zeros(edge_count,1);
chord_count=0;

for edge_row=1:edge_count
    u=edge_nodes(edge_row,1);
    v=edge_nodes(edge_row,2);
    root_u=find_root(parent,u);
    root_v=find_root(parent,v);
    if root_u~=root_v
        parent(root_v)=root_u;
        tree_edge(edge_row)=true;
    else
        chord_count=chord_count+1;
        chord_rows(chord_count)=edge_row;
    end
end
chord_rows=chord_rows(1:chord_count);

adjacent_nodes=cell(node_count,1);
adjacent_edges=cell(node_count,1);
for edge_row=find(tree_edge)'
    u=edge_nodes(edge_row,1);
    v=edge_nodes(edge_row,2);
    adjacent_nodes{u}(end+1)=v;
    adjacent_edges{u}(end+1)=edge_row;
    adjacent_nodes{v}(end+1)=u;
    adjacent_edges{v}(end+1)=edge_row;
end

cycle_matrix=zeros(chord_count,edge_count);
for cycle_row=1:chord_count
    chord=chord_rows(cycle_row);
    source=edge_nodes(chord,1);
    target=edge_nodes(chord,2);
    [path_edges,found]=tree_path(source,target,adjacent_nodes,adjacent_edges,node_count);
    assert(found,'The fundamental-cycle tree path could not be reconstructed.');
    cycle_matrix(cycle_row,[path_edges(:);chord])=1;
end
end

function root=find_root(parent,node)
root=node;
while parent(root)~=root
    root=parent(root);
end
end

function [path_edges,found]=tree_path(source,target,adjacent_nodes,adjacent_edges,node_count)
visited=false(node_count,1);
previous_node=zeros(node_count,1);
previous_edge=zeros(node_count,1);
queue=zeros(node_count,1);
queue_head=1;
queue_tail=1;
queue(1)=source;
visited(source)=true;

while queue_head<=queue_tail&&~visited(target)
    current=queue(queue_head);
    queue_head=queue_head+1;
    for neighbor_index=1:numel(adjacent_nodes{current})
        neighbor=adjacent_nodes{current}(neighbor_index);
        if ~visited(neighbor)
            visited(neighbor)=true;
            previous_node(neighbor)=current;
            previous_edge(neighbor)=adjacent_edges{current}(neighbor_index);
            queue_tail=queue_tail+1;
            queue(queue_tail)=neighbor;
        end
    end
end

found=visited(target);
path_edges=zeros(0,1);
if ~found,return;end
current=target;
while current~=source
    path_edges(end+1,1)=previous_edge(current); %#ok<AGROW>
    current=previous_node(current);
end
end
