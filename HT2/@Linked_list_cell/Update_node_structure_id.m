function [] = Update_node_structure_id(obj, value, id)
    %This function search for a node in the linked list, return True or
    %False

    %Check if the linked list is empty
    if obj.first_node.value == 'Null'
        fprintf('The list is empty');
    elseif obj.last_node.value == value 
        obj.last_node.cell_struct.id = id;
    else        
        obj.current = obj.first_node;
        while (obj.current.value ~= value && obj.current ~= obj.last_node)
            obj.current = obj.current.next;
        end

        if obj.current.value == value
            obj.current.cell_struct.id = id;
        end
    end
end