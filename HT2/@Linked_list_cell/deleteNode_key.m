function [] = deleteNode_key(obj,value)
    %Delete a node of the linked list

    %create a helper(auxuliar) node
    helper_node = obj.first_node;

     %Check if the linked list is empty
    if helper_node.value == 'Null'
        fprintf('The list is empty');
    elseif (helper_node.value == value)
        %if the node is at the beginning of the list 
        
        %If there is only one element
        if obj.first_node == obj.last_node
            %Update the first node and end node, the list is empty
            obj.first_node.value = 'Null';
            obj.first_node.cell_struct.i = 0;
            obj.first_node.cell_struct.j = 0;
            obj.first_node.cell_struct.l = 0;
            obj.first_node.cell_struct.key = 'Null';
            obj.first_node.cell_struct.id = 0;
            obj.first_node.next = 'Null';
            obj.first_node.previous = 'Null';
            obj.last_node = obj.first_node;
            obj.head.next = obj.first_node;
        else
            %Update the first node
            obj.first_node = helper_node.next;
            obj.first_node.previous = 'Null';
            obj.head.next = obj.first_node;
        end      
       
        %Delete the helper node
        clear helper_node;
    elseif (obj.last_node.value == value)
        %disp('True')
        %If node is at the end of the list

        %create a helper(auxuliar) node
        helper_node_2 = obj.last_node; 

        if obj.first_node == obj.last_node
            %Update the first node and end node, the list is empty
            obj.first_node.value = 'Null';
            obj.first_node.cell_struct.i = 0;
            obj.first_node.cell_struct.j = 0;
            obj.first_node.cell_struct.l = 0;
            obj.first_node.cell_struct.key = 'Null';
            obj.first_node.cell_struct.id = 0;
            obj.first_node.next = 'Null';
            obj.first_node.previous = 'Null';
            obj.last_node = obj.first_node;
            obj.head.next = obj.first_node;
        else
            %Update the end node
            obj.last_node = obj.last_node.previous;
            obj.last_node.next = 'Null';
        end 
        
        %Delete the helper node
        clear helper_node_2;
        clear helper_node;
    else
        %In other cases

        while helper_node.next ~= 'Null'
            if helper_node.next.value == value
                helper_node_2 = helper_node.next;
                helper_node.next = helper_node_2.next;
                helper_node_2.next.previous = helper_node;
                clear helper_node_2;
            end
            
            helper_node = helper_node.next;

        end
        clear helper_node;
    end

    %counter the number nodes in the linked list
    obj.head.counter =  obj.head.counter - 1;
end