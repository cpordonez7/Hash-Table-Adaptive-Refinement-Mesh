%Test of the linked list cell

%Creates the nodes
cell1 = create_cell(1,1,0);
cell1.key = 1;
Node1 = create_node_cell(cell1);
cell2 = create_cell(2,2,0);
cell2.key = 2;
cell3 = create_cell(3,3,1);
cell3.key = 3;
cell4 = create_cell(4,4,1);
cell4.key = 4;

%Creates the linked list
llcell_Obj = Linked_list_cell(Node1);

%%%%%%%%%%%%%%%%
%Add nll2_Obj.print_list;odes to the linked list
addNode_cell(llcell_Obj,cell2);
addNode_cell(llcell_Obj,cell3);
addNode_cell(llcell_Obj,cell4);

llcell_Obj.print_list;

%Delete a node fron the linked list
% % deleteNode_key(ll2_Obj,4)
% % 
% % ll2_Obj.print_list;
% % 
% % searchNode(ll2_Obj,4)
% % 
% % print_list(ll2_Obj)
% % 
% % ll2_Obj.makeEmpty;
% % 
% % print_list(ll2_Obj)


