function [] = print_list(obj)
    %Print the nodes of the linked list
    obj.current = obj.first_node;
    if (obj.current.value == 'Null')
        fprintf('The list is empty');
    elseif (obj.current.next == 'Null')
        fprintf('The list has only one element and this is %d',obj.current.value);
    else 
        fprintf('The elements of the list are: [');
        while obj.current.next ~= 'Null'
            fprintf('%d \t',obj.current.value);
            obj.current = obj.current.next;
        end
        fprintf('%d ]\n\n',obj.last_node.value);
    end
end