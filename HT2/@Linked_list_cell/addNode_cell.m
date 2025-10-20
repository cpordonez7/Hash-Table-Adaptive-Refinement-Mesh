function [] = addNode_cell(obj,cell_struct)
    %Add a new node to the end of the linked list

    %Creates a new node, input a struct
    New_node = create_node_cell(cell_struct);
    % disp(New_node.cell_struct)

    %Check if the linked list is empty
    if obj.first_node.value == 'Null'
        obj.first_node = New_node;
        obj.last_node = obj.first_node;
        obj.current = obj.first_node;
        obj.head.next = obj.first_node;
    else
        New_node.previous = obj.last_node;
        obj.last_node.next = New_node;
        obj.last_node = New_node;
    end

    %counter the number nodes in the linked list
    obj.head.counter =  obj.head.counter + 1;
end