function [bool] = isempty(obj)
    %This function checks if the linked list is empty

    %Check if the linked list is empty
    if obj.first_node.value == 'Null'
        bool = 'True';
    else 
        bool = 'False';
    end

 
end