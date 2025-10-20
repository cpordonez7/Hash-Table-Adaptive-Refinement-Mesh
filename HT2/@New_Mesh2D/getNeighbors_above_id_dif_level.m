function [List_Neighbors_above_return] = getNeighbors_above_id_dif_level(HTObj,id)
    %Return the ids of the all neighbors cell and the difference of level
    %of the input cell in this direction.

    % fprintf('%d\n',id);
    List_Neighbors_above = [];

    %call the values of (i,j,l) of the cell with id enumeration
    % p = HTObj.Enumeration_vector{id}(1);
    % key = HTObj.Enumeration_vector{id}(2);
    % 
    % %Return the struct of cell with this key
    % cell_id = HTObj.HT{p}.get_cell_Nodo(key);
    cell_id = HTObj.getCell_id(id);           %input id cell

    %Indixes of the cell 
    ind_i = cell_id.i;
    ind_j = cell_id.j;
    l = cell_id.l;
    N_level = HTObj.N*2^(l);
    M_level = HTObj.M*2^(l);

    %case 1, cell belown-left  
    if (ind_i == 0) && (ind_j == 0)             
        List_Neighbors_above = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_above);

        List_Neighbors_above_return(1,:) = List_Neighbors_above; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_above_return(1,:))
            id_j = List_Neighbors_above_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_above_return(2,:) = List_level_difference;  % only the level difference
    end

    %Case 2, cell belown
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == 0)                 
        List_Neighbors_above = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_above);

        List_Neighbors_above_return(1,:) = List_Neighbors_above; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_above_return(1,:))
            id_j = List_Neighbors_above_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_above_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 3, cell rigth belown
    if (ind_i == N_level-1) && (ind_j == 0)
        List_Neighbors_above = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_above);

        List_Neighbors_above_return(1,:) = List_Neighbors_above; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_above_return(1,:))
            id_j = List_Neighbors_above_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_above_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 4, cell left
    if (ind_i == 0) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_above = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_above);

        List_Neighbors_above_return(1,:) = List_Neighbors_above; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_above_return(1,:))
            id_j = List_Neighbors_above_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_above_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 5, cell rigth
    if (ind_i == N_level-1) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_above = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_above);

        List_Neighbors_above_return(1,:) = List_Neighbors_above; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_above_return(1,:))
            id_j = List_Neighbors_above_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_above_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 6, cell left above  
    if (ind_i == 0) && (ind_j == M_level-1)
         % disp('There are not cells in this direction L-A');   

         %CELL NULL
        List_Neighbors_above_return(1,:) = 0; % only the id cells
        List_Neighbors_above_return(2,:) = 0; % only the level difference
    end

    %Caso 7, cell above
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == M_level-1)
        % disp('There are not cells in this direction A'); 

        %CELL NULL
        List_Neighbors_above_return(1,:) = 0; % only the id cells
        List_Neighbors_above_return(2,:) = 0; % only the level difference
    end

    %Caso 8, cell rigth above
    if (ind_i == N_level-1) && (ind_j == M_level-1)
        % disp('There are not cells in this direction R-A');

        %CELL NULL
        List_Neighbors_above_return(1,:) = 0; % only the id cells
        List_Neighbors_above_return(2,:) = 0; % only the level difference
    end

    %Caso 9, cell center
    if (0 < ind_i && ind_i < N_level-1) && (0 < ind_j && ind_j < M_level-1)
       List_Neighbors_above = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_above);

       List_Neighbors_above_return(1,:) = List_Neighbors_above; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_above_return(1,:))
            id_j = List_Neighbors_above_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_above_return(2,:) = List_level_difference;  % only the level difference
    end

end