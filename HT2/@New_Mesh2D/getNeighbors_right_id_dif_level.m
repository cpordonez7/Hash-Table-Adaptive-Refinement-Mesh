function [List_Neighbors_right_return] = getNeighbors_right_id_dif_level(HTObj,id)
    %Return the ids of the all neighbors cell and the difference of level
    %of the input cell in this direction.

    % fprintf('%d\n',id);
    List_Neighbors_right = [];

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
        List_Neighbors_right = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_right);

        List_Neighbors_right_return(1,:) = List_Neighbors_right; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_right_return(1,:))
            id_j = List_Neighbors_right_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_right_return(2,:) = List_level_difference;  % only the level difference
    end

    %Case 2, cell belown
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == 0)                 
        List_Neighbors_right = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_right);

        List_Neighbors_right_return(1,:) = List_Neighbors_right; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_right_return(1,:))
            id_j = List_Neighbors_right_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_right_return(2,:) = List_level_difference;  % only the level difference
    end

    %Case 3, cell rigth belown
    if (ind_i == N_level-1) && (ind_j == 0)
        % disp('There are not cells in this direction R-B'); 

        %CELL NULL
        List_Neighbors_right_return(1,:) = 0; % only the id cells
        List_Neighbors_right_return(2,:) = 0; % only the level difference
    end

    %Case 4, cell left
    if (ind_i == 0) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_right = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_right); 

        List_Neighbors_right_return(1,:) = List_Neighbors_right; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_right_return(1,:))
            id_j = List_Neighbors_right_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_right_return(2,:) = List_level_difference;  % only the level difference
    end

    %Case 5, cell rigth
    if (ind_i == N_level-1) && (0 < ind_j && ind_j < M_level-1)
        % disp('There are not cells in this direction R');

        %CELL NULL
        List_Neighbors_right_return(1,:) = 0; % only the id cells
        List_Neighbors_right_return(2,:) = 0; % only the level difference
    end

    %Case 6, cell left above  
    if (ind_i == 0) && (ind_j == M_level-1)
        List_Neighbors_right = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_right); 

        List_Neighbors_right_return(1,:) = List_Neighbors_right; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_right_return(1,:))
            id_j = List_Neighbors_right_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_right_return(2,:) = List_level_difference;  % only the level difference
    end

    %Case 7, cell above
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == M_level-1)
        List_Neighbors_right = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_right);

        List_Neighbors_right_return(1,:) = List_Neighbors_right; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_right_return(1,:))
            id_j = List_Neighbors_right_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_right_return(2,:) = List_level_difference;  % only the level difference
    end

    %Case 8, cell rigth above
    if (ind_i == N_level-1) && (ind_j == M_level-1)
        % disp('There are not cells in this direction R-A');

        %CELL NULL
        List_Neighbors_right_return(1,:) = 0; % only the id cells
        List_Neighbors_right_return(2,:) = 0; % only the level difference
    end

    %Case 9, cell center
    if (0 < ind_i && ind_i < N_level-1) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_right = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_right); 

        List_Neighbors_right_return(1,:) = List_Neighbors_right; % only the id cells
        List_level_difference = [];
        for j=1:length(List_Neighbors_right_return(1,:))
            id_j = List_Neighbors_right_return(1,j);

            cell_id_j = HTObj.getCell_id(id_j);         
            List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference 
        end
        List_Neighbors_right_return(2,:) = List_level_difference;  % only the level difference
    end

end