function [List_Neighbors_stencil_RETURN] = getNeighbors_stencil(HTObj, id)
    %Return the id and level diference between the referent cell and her neighbors
    %Contour cell have id zero, level diference zero, is a cell NULL

    List_Neighbors_stencil = [];

    %call the values of (i,j,l) of the cell with id enumeration
    cell_id = HTObj.getCell_id(id);           %input id cell

    %Indixes of the cell 
    ind_i = cell_id.i;
    ind_j = cell_id.j;
    l = cell_id.l;
    N_level = HTObj.N*2^(l);
    M_level = HTObj.M*2^(l);

    % fprintf('i: %d \t j: %d \t l: %d\n',ind_i,ind_j,l);
    %Analyze all cases according to the position of the cell in the mesh

    %case 1 cell belown-left
    if (ind_i == 0) && (ind_j == 0)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_stencil);
        
        List_Neighbors_stencil_RETURN(1,:) = [[0 0] List_Neighbors_stencil]; % only the id cells, incluid the border cells (B, L)
        
        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end    
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference
    end

    %Case 2 cell belown
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == 0)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_stencil);

        List_Neighbors_stencil_RETURN(1,:) = [0 List_Neighbors_stencil]; % only the id cells
        
        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference
    end

    %Case 3 cell rigth belown
    if (ind_i == N_level-1) && (ind_j == 0)
        %Here below cell=0
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_stencil);
        List_Neighbors_left = List_Neighbors_stencil;  %copy
        %Here Rigth cell=0
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_stencil);


        n_v = size(List_Neighbors_left,2);
        List_Neighbors_stencil(:,1:n_v) = [];  %delete the left cell

        List_Neighbors_stencil_RETURN(1,:) = [0 List_Neighbors_left 0  List_Neighbors_stencil]; % only the id cells

        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference

    end

    %Case 4 cell left
    if (ind_i == 0) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_stencil);
        Lista_Neighbors_below = List_Neighbors_stencil;  %copy
        %Here Left cell=0
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_stencil);

        n_v = size(Lista_Neighbors_below,2);
        List_Neighbors_stencil(:,1:n_v) = [];  %delete the belows cell
        
        List_Neighbors_stencil_RETURN(1,:) = [Lista_Neighbors_below 0  List_Neighbors_stencil];   % only the id cells
        
        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference

    end

    %Case 5 cell rigth
    if (ind_i == N_level-1) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_stencil);
        List_Neighbors_below_left = List_Neighbors_stencil;  %copy
        %Here rigth cell=0
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_stencil);

        n_v = size(List_Neighbors_below_left,2);
        List_Neighbors_stencil(:,1:n_v) = [];  %delete the right cell
        
        List_Neighbors_stencil_RETURN(1,:) = [List_Neighbors_below_left 0  List_Neighbors_stencil];   % only the id cells

        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference
        
    end

    %Case 6, cell left above
    if (ind_i == 0) && (ind_j == M_level-1)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_stencil);
        Lista_Neighbors_below_2 = List_Neighbors_stencil;  %copy
        %Here left cell=0
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_stencil);
        %Here above cell=0

        n_v = size(Lista_Neighbors_below_2,2);
        List_Neighbors_stencil(:,1:n_v) = [];  %delete the above cell
        
        List_Neighbors_stencil_RETURN(1,:) = [Lista_Neighbors_below_2 0  List_Neighbors_stencil 0];   % only the id cells
        
        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference

    end

    %Case 7, cell above
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == M_level-1)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_stencil);
        %Here above cell=0 
        
        List_Neighbors_stencil_RETURN(1,:) = [List_Neighbors_stencil 0]; % only the id cells

        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference
    end

    %Case 8, cell rigth above
    if (ind_i == N_level-1) && (ind_j == M_level-1)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_stencil);
        %Here rigth cell=0
        %Here above cell=0
        
        List_Neighbors_stencil_RETURN(1,:) = [List_Neighbors_stencil [0 0]]; % only the id cells

        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference
    end

    %Case 9, cell center
    if (0 < ind_i && ind_i < N_level-1) && (0 < ind_j && ind_j < M_level-1)
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_stencil);
        List_Neighbors_stencil = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_stencil);

        List_Neighbors_stencil_RETURN(1,:) = List_Neighbors_stencil; % only the id cells

        List_level_difference = [];
        for j=1:length(List_Neighbors_stencil_RETURN(1,:))
            id_j = List_Neighbors_stencil_RETURN(1,j);

            if id_j == 0
                List_level_difference = [List_level_difference 0];
            else 
                cell_id_j = HTObj.getCell_id(id_j);         
                List_level_difference = [List_level_difference cell_id_j.l-l]; % only the level difference
            end
        end
        List_Neighbors_stencil_RETURN(2,:) = List_level_difference; % only the level difference
    end

end