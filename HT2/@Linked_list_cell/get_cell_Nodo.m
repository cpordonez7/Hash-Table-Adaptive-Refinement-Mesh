function [Cell_struct] = get_cell_Nodo(obj, value)
    %This function search for a node in the linked list and returns the
    %structure-type cell of that node

    %Check if the linked list is empty
    if obj.first_node.value == 'Null'
        fprintf('The list is empty');
        % Cell_struct = obj.first_node.cell_struct;
    else
        obj.current = obj.first_node;
        while (obj.current.value ~= value && obj.current ~= obj.last_node)
            obj.current = obj.current.next;
        end

        if obj.current.value == value
            % bool = Obj.current; 
            Cell_struct = obj.current.cell_struct;
        else
            fprintf('The cell is not in the linked list');
            % bool = 'False';
        end
    end
end