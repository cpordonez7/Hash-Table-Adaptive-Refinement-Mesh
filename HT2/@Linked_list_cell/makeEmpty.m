function [] = makeEmpty(obj)
    %This function makes empty the linked list 

    if obj.first_node.value == 'Null'
        printf('The list is empty');
    end


    while obj.first_node ~= obj.last_node
        helper_node = obj.first_node;
        obj.first_node = obj.first_node.next;
        % 
        % %Delete current 
        clear helper_node;
    end

    if obj.first_node == obj.last_node
        obj.first_node.value = 'Null';
        obj.first_node.cell_struct.i = 0;
        obj.first_node.cell_struct.j = 0;
        obj.first_node.cell_struct.l = 0;
        obj.first_node.cell_struct.key = 'Null';
        obj.first_node.cell_struct.id = 0;
        obj.first_node.next = 'Null';
        obj.first_node.previous = 'Null';
        obj.last_node = obj.first_node;
    end

    obj.head.next = obj.first_node;
    obj.head.counter = 0;
end