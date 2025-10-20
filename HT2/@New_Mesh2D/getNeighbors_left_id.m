function [List_Neighbors_left] = getNeighbors_left_id(HTObj,id)
    %Return the ids of the all neighbors cell of the input cell in this
    %direction.

    % fprintf('%d\n',id);
    List_Neighbors_left = [];

    %call the values of (i,j,l) of the cell with id enumeration
    cell_id = HTObj.getCell_id(id);           %input id cell

    %Indixes of the cell 
    ind_i = cell_id.i;
    ind_j = cell_id.j;
    l = cell_id.l;
    N_level = HTObj.N*2^(l);
    M_level = HTObj.M*2^(l);

    %case 1, cell belown-left  
    if (ind_i == 0) && (ind_j == 0)             
        disp('There are not cells in this direction B-L');    
    end

    %Case 2, cell belown
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == 0)                 
        List_Neighbors_left = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_left);      
    end

    %Case 3, cell rigth belown
    if (ind_i == N_level-1) && (ind_j == 0)
        List_Neighbors_left = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_left);
    end

    %Case 4, cell left
    if (ind_i == 0) && (0 < ind_j && ind_j < M_level-1)
        disp('There are not cells in this direction L');
    end

    %Case 5, cell rigth
    if (ind_i == N_level-1) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_left = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_left);
    end

    %Case 6, cell left above  
    if (ind_i == 0) && (ind_j == M_level-1)
         disp('There are not cells in this direction A');          
    end

    %Case 7, cell above
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == M_level-1)
        List_Neighbors_left = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_left); 
    end

    %Case 8, cell rigth above
    if (ind_i == N_level-1) && (ind_j == M_level-1)
        List_Neighbors_left = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_left);
    end

    %Case 9, cell center
    if (0 < ind_i && ind_i < N_level-1) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_left = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_left);   
    end

end