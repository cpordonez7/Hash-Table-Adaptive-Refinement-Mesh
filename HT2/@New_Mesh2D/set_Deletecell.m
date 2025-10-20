function [] = set_Deletecell(HTObj,cell_str)
    %set_Deletecell: This function delete a cell in the HT, in the assigned
    %position according to its key. The list of cells in that position
    %reduces by one unit.

    %Input: 
    %HTObj -> Object Hast Table
    %cell_str -> cell structure 
 
    %Indices (i,j) and level.
    % ind_i = cell_str.i;
    % ind_j = cell_str.j;
    % l = cell_str.l;

    key = HTObj.key_Calculation(cell_str);     %cell key
    % p = mod(key,HTObj.size_HT);              %cell poition in the HT
    p = HTObj.hast_function(key);

    %Process to delete a cell

    %Check if the linked has the cell and delete it. Using the linked list
    %operation
    HTObj.HT{p+1}.deleteNode_key(key);   

    HTObj.Cell_counter_level(cell_str.l + 1) = HTObj.Cell_counter_level(cell_str.l + 1) - 1;

end