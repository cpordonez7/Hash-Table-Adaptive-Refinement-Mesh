function [] = ThickenCell(HTObj,cell_struct)
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
    key_cell1 = HTObj.key_Calculation_i_j_l(ind_i,ind_j,l);
    p = HTObj.hast_function(key_cell1);
    cell1 =HTObj.HT{p+1}.get_cell_Nodo(key_cell1);
    Id_cell1 = cell1.id;
    set_Deletecell(HTObj,cell1);

    %Update the HT, add the new node
    cell_str = create_cell(ind_i_new,ind_j_new,l-1);
    set_Insertcell(HTObj,cell_str)


    %Subtract one because inserting adds al final
    HTObj.Enumeration_counter = HTObj.Enumeration_counter - 1 ;

    %%Delete the four fines cells  
    %
    key_cell2 = HTObj.key_Calculation_i_j_l(ind_i+1,ind_j,l);
    p = HTObj.hast_function(key_cell2);
    cell2 =HTObj.HT{p+1}.get_cell_Nodo(key_cell2);
    Id_cell2 = cell2.id;
    
    set_Deletecell(HTObj,cell2);
    %
    key_cell3 = HTObj.key_Calculation_i_j_l(ind_i,ind_j+1,l);
    p = HTObj.hast_function(key_cell3);
    cell3 =HTObj.HT{p+1}.get_cell_Nodo(key_cell3);
    Id_cell3 = cell3.id;

    set_Deletecell(HTObj,cell3);

    %
    key_cell4 = HTObj.key_Calculation_i_j_l(ind_i+1,ind_j+1,l);
    p = HTObj.hast_function(key_cell4);
    cell4 =HTObj.HT{p+1}.get_cell_Nodo(key_cell4);
    Id_cell4 = cell4.id;

    set_Deletecell(HTObj,cell4);

    %Agrega en un vector los Id_i
    vec_Id = [Id_cell1, Id_cell2, Id_cell3, Id_cell4];
    vec_Id = sort(vec_Id);

    %Update Id thincked cell
    key_cell_str = HTObj.key_Calculation(cell_str);
    p = HTObj.hast_function(key_cell_str); 

    %Update id of the new cell, leaving the Id of the first fine cell 
    HTObj.HT{p+1}.Update_node_structure_id(key_cell_str,vec_Id(1));
    % cellll = HTObj.HT{p+1}.get_cell_Nodo(key_cell_str);
    % disp(cellll)

    total_cell= HTObj.Enumeration_counter;  

    if vec_Id(2) < total_cell-2
        new_Id2 = vec_Id(2);
        cell_fin_1 = HTObj.getCell_id(total_cell-2);
        %Update id of the new cell, leaving the Id of the first fine cell 
        key_cell2 = cell_fin_1.key;
        p = HTObj.hast_function(key_cell2);
        HTObj.HT{p+1}.Update_node_structure_id(key_cell2,new_Id2);
        if vec_Id(3) < total_cell-1
            new_Id3 = vec_Id(3);
            cell_fin_2 = HTObj.getCell_id(total_cell-1);
            %Update id of the new cell, leaving the Id of the first fine cell 
            key_cell3 = cell_fin_2.key;
            p = HTObj.hast_function(key_cell3);
            HTObj.HT{p+1}.Update_node_structure_id(key_cell3,new_Id3);
            if vec_Id(4) < total_cell
                new_Id4 = vec_Id(4);
                cell_fin_3 = HTObj.getCell_id(total_cell);
                %Update id of the new cell, leaving the Id of the first fine cell 
                key_cell4 = cell_fin_3.key;
                p = HTObj.hast_function(key_cell4);
                HTObj.HT{p+1}.Update_node_structure_id(key_cell4,new_Id4);
            end
        end    
    end

    %Eliminate three because one remains fixed
    HTObj.Enumeration_counter = HTObj.Enumeration_counter - 3;

    %Updates all propietates of the object
    set_Number_cells_Mesh(HTObj);

        % HTObj.counter_enumeration = HTObj.counter_enumeration - 3;   
        % AQUI CUANDO ENGROSA VER COMO QUEDA LA ENUMAEACION, PORQUE PUEDE
        % ENGROSAR UNA QUE NO SEA LA ULTIMA CELULA REFINADA ENTONCES HAY
        % QUE ACTUALIZAR LA ENUMERACION
end