%Test of the create_node_cell

% cell_struct = create_cell(1,1,0);
% key = 7;
% cell_struct.key = key;
% % 
% % 
% Node1 = create_node_cell(cell_struct);


%%%%%%%%%%%%%
HT = cell(2,1);  
size_HT_p = size(HT{1},2);
disp(size_HT_p);
cel1 = create_cell(1,1,0);
node1 = create_node_cell(cel1);
HT{1}=Linked_list_cell(node1);
% addNode_key(HT{1},9);
