function [List_Neighbors_below_return] = getNeighbors_below_id_dif_level(HTObj,id)
    %Return the ids of the all neighbors cell and the difference of level
    %of the input cell in this direction.

    % fprintf('%d\n',id);
    List_Neighbor_below = [];


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
        % disp('There are not cells in this direction B-L'); 

        %CELL NULL
        List_Neighbors_below_return(1,:) = 0; % only the id cells
        List_Neighbors_below_return(2,:) = 0; % only the level difference
    end

    %Case 2, cell belown
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == 0)                 
        % disp('There are not cells in this direction B');

        %CELL NULL
        List_Neighbors_below_return(1,:) = 0; % only the id cells
        List_Neighbors_below_return(2,:) = 0; % only the level difference
    end

    %Caso 3, cell rigth belown
    if (ind_i == N_level-1) && (ind_j == 0)
        % disp('There are not cells in this direction R-B');

        %CELL NULL
        List_Neighbors_below_return(1,:) = 0; % only the id cells
        List_Neighbors_below_return(2,:) = 0; % only the level difference
    end

    %Caso 4, cell left
    if (ind_i == 0) && (0 < ind_j && ind_j < M_level-1)
        %Celulas abajo
        List_Neighbor_below = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbor_below); 

        List_Neighbors_below_return(1,:) = List_Neighbor_below; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_below_return(1,:))
            id_j = List_Neighbors_below_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_below_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 5, cell rigth
    if (ind_i == N_level-1) && (0 < ind_j && ind_j < N_level-1)
        %Celulas abajo
        List_Neighbor_below = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbor_below);

        List_Neighbors_below_return(1,:) = List_Neighbor_below; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_below_return(1,:))
            id_j = List_Neighbors_below_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_below_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 6, cell left above  
    if (ind_i == 0) && (ind_j == M_level-1)
         %Celulas abajo
         List_Neighbor_below = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbor_below);

         List_Neighbors_below_return(1,:) = List_Neighbor_below; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_below_return(1,:))
            id_j = List_Neighbors_below_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_below_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 7, cell above
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == M_level-1)
        %Celulas abajo
        List_Neighbor_below = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbor_below);

        List_Neighbors_below_return(1,:) = List_Neighbor_below; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_below_return(1,:))
            id_j = List_Neighbors_below_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_below_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 8, cell rigth above
    if (ind_i == N_level-1) && (ind_j == M_level-1)
        %Celulas abajo
        List_Neighbor_below = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbor_below);

        List_Neighbors_below_return(1,:) = List_Neighbor_below; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_below_return(1,:))
            id_j = List_Neighbors_below_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_below_return(2,:) = List_level_difference;  % only the level difference
    end

    %Caso 9, cell center
    if (0 < ind_i && ind_i < N_level-1) && (0 < ind_j && ind_j < M_level-1)
        %Celulas abajo
        List_Neighbor_below = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbor_below); 

        List_Neighbors_below_return(1,:) = List_Neighbor_below; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_below_return(1,:))
            id_j = List_Neighbors_below_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_below_return(2,:) = List_level_difference;  % only the level difference
    end

end