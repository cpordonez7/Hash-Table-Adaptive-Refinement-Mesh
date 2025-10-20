classdef Linked_list_cell < handle
    %   This class, creates a list of the nodes.

%-----------------------------------------------------------------
%               PROPERTIES SECTION
%-----------------------------------------------------------------
    properties (Access = public)
        head           %Head node
        first_node     %Initial node in the linked list 
        last_node      %Last node in the linked list
        current        %current node
    end

%-----------------------------------------------------------------
%               METHODS SECTION
%-----------------------------------------------------------------
    methods
        function [obj] = Linked_list_cell(New_node)
            %UNTITLED2 Construct an instance of this class
            %The linked list is initialized with a value from input.
            obj.first_node = New_node;
            obj.last_node = obj.first_node;
            obj.current = obj.first_node;
            
            %Head node
            obj.head.next = obj.first_node;
            obj.head.counter = 1;
        end
        
%-----------------------------------------------------------------
%               INTERFACES SECTION 
%-----------------------------------------------------------------
    [] = addNode_cell(obj,struct)     %Add a node from the linked list, with input its cell type struct      
    [] = deleteNode_key(obj,value)   %Delete a node from the linked list
    [bool] = searchNode(Obj, value) %Search a node from the linked list
    [Cell_struct] = get_cell_Nodo(Obj, value) %Return the node from the linked list with this key
    [] = print_list(obj)            %Print the linked list
    [bool] = isempty(obj)           %Check if the linked list is empty
    [] = makeEmpty(obj)             %Make empty the linked list
    [] = Update_node_structure_id(obj, value, id)   %Updata id of the cell in the node with key value
    [node_counter] = size_list(obj);   %Return the size of the linked list     


    end %end metodos

end %end class