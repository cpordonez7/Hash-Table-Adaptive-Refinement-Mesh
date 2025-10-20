function [node_counter] = size_list(obj)
    %Return the number of nodes in the linked list
    node_counter = obj.head.counter;
    % node_counter = 0;
    % obj.current = obj.first_node;
    % if (obj.current.value == 'Null')
    %     node_counter = 0;
    % elseif (obj.current.next == 'Null')
    %     node_counter = node_counter + 1;
    % else 
    %     while obj.current.next ~= 'Null'
    %         node_counter = node_counter + 1;
    %         obj.current = obj.current.next;
    %     end
    %     node_counter = node_counter + 1;
    % end
end