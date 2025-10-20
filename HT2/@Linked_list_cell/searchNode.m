function [bool] = searchNode(obj, value)
    %This function search for a node in the linked list, return True or
    %False

    %Check if the linked list is empty
    if obj.first_node.value == 'Null'
        %fprintf('The list is empty');
        bool = 'False';
    else
        obj.current = obj.first_node;
        while (obj.current.value ~= value && obj.current ~= obj.last_node)
            obj.current = obj.current.next;
        end

        if obj.current.value == value
            % bool = Obj.current;
            bool = 'True';
        else
            % fprintf('The node is not in the linked list');
            bool = 'False';
        end
    end
end