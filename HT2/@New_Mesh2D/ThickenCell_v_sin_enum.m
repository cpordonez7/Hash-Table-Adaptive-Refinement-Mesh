function [] = ThickenCell_v_sin_enum(HTObj,cell_struct)
   %This function refine a cell of the mesh, 
    %input: cell struct that you want to refine

    l = cell_struct.l;
    % Id = cell_struct.id;
    % key = cell_struct.key;

    %indixes of the new cell left-below
    ind_i = cell_struct.i;
    ind_j = cell_struct.j;

    %Add a new cell and delete the cell that form the cell to thicken

    if l <= 0
        disp(l)
        fprintf('Error, This cell is in the base level. It is not possible to thicken');
        return;
    end

    if mod(ind_i,2) == 0
        ind_i_new = round(ind_i/2);
    else 
        ind_i_new = round((ind_i-1)/2);
    end

    if mod(ind_j,2) == 0
        ind_j_new = round(ind_j/2);
    else 
        ind_j_new = round((ind_j-1)/2);
    end

    %Index fines cell
    ind_i = 2*ind_i_new;
    ind_j = 2*ind_j_new;
    %first cell
    cell1 = create_cell(ind_i,ind_j,l);

    key_cell1 = HTObj.key_Calculation(cell1);
    p = HTObj.hast_function(key_cell1);
    cell1 =HTObj.HT{p+1}.get_cell_Nodo(key_cell1);
    Id_cell1 = cell1.id;
    set_Deletecell(HTObj,cell1);
    % disp(Id_cell1)

    %Update the HT, add the new node
    cell_str = create_cell(ind_i_new,ind_j_new,l-1);
    set_Insertcell(HTObj,cell_str)

    %Update Id thincked cell
    key_cell_str = HTObj.key_Calculation(cell_str);
    p = HTObj.hast_function(key_cell_str); 

    %Update id of the new cell, leaving the Id of the first fine cell 
    HTObj.HT{p+1}.Update_node_structure_id(key_cell_str,Id_cell1);
    % cellll = HTObj.HT{p+1}.get_cell_Nodo(key_cell_str);
    % disp(cellll)

    %Subtract one because inserting adds
    HTObj.Enumeration_counter = HTObj.Enumeration_counter - 1 ;

    %Delete the four fines cells  
    cell2 = create_cell(ind_i+1,ind_j,l);
    set_Deletecell(HTObj,cell2);
    cell3 = create_cell(ind_i,ind_j+1,l);
    set_Deletecell(HTObj,cell3);
    cell4 = create_cell(ind_i+1,ind_j+1,l);
    set_Deletecell(HTObj,cell4);

        % HTObj.counter_enumeration = HTObj.counter_enumeration - 3;   
        % AQUI CUANDO ENGROSA VER COMO QUEDA LA ENUMAEACION, PORQUE PUEDE
        % ENGROSAR UNA QUE NO SEA LA ULTIMA CELULA REFINADA ENTONCES HAY
        % QUE ACTUALIZAR LA ENUMERACION
end