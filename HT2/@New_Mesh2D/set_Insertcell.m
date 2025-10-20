function [] = set_Insertcell(HTObj,cell_str)
    %set_Insertcell: This function insert a cell in the HT, in the assigned
    %position according to its key. The cell is insert the list end.
    %
    %Input: 
    %HTObj -> Object Hast Table
    %cell_str -> cell structure 
 
    %Indices (i,j) and level.
    % ind_i = cell_str.i;
    % ind_j = cell_str.j;
    % l = cell_str.l;

    % fprintf('%.f \t %.f \n',x1,y1);
    HTObj.Cell_counter_level(cell_str.l + 1) = HTObj.Cell_counter_level(cell_str.l + 1) + 1;

    %Cell Id, enumeration for level, according to the order of insertion (file)
    HTObj.Enumeration_counter = HTObj.Enumeration_counter + 1;
    cell_str.id = HTObj.Enumeration_counter;
    % fprintf('id:%d \n',cell_str.id);

    key = HTObj.key_Calculation(cell_str);     %cell key
    % fprintf('key:%d \n',key);
    cell_str.key = key;                        %Update cell key
    % p = mod(key,HTObj.size_HT);                 
    p = HTObj.hast_function(cell_str.key);     %cell poition in the HT
    % fprintf('p:%d \n',p);
    % fprintf('cell_key:%d \n',cell_str.key);
    % fprintf('cellinsert key:%d \t p:%d \n',key,p);

    % if p == 0 
    %     p = 1;
    % end
    Node_cell = create_node_cell(cell_str);    %Crearte a Node 
    % disp(Node_cell.cell_struct)
    %list size at position p
    size_HT_p = size(HTObj.HT{p+1},2);

    %Check if there is a linked list in position p or it is empty
    if (size_HT_p==0)
        %If size_HT_P = 0 -> the list is empty
        HTObj.HT{p+1} = Linked_list_cell(Node_cell);
    else
        %If size_HT_P ~= 0 -> there is collision, in that position there is at least one cell
        %disp ('There is a collision');
        HTObj.HT{p+1}.addNode_cell(cell_str);
    end
    
end  
