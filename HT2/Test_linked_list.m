%Test of the linked list

%Creates the nodes
Node1 = create_node(1);
Node2 = create_node(2);
Node3 = create_node(3);
Node4 = create_node(4);

%Creates the linked list
LL_Obj = Linked_list(Node1);

%Add nodes to the linked list
addNode(LL_Obj,Node2);
addNode(LL_Obj,Node3);
addNode(LL_Obj,Node4);

%Print linkes list
print_list(LL_Obj)


%Delete a node fron the linked list
deleteNode(LL_Obj,Node1)
print_list(LL_Obj)

searchNode(LL_Obj,Node2.value)

