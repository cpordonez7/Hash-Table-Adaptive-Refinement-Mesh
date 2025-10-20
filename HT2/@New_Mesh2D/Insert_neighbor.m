function [List_Neighbors_id] = Insert_neighbor(HTObj,ind_i,ind_j,l,Dir,List_Neighbors_id)
    %This function adds the cells that are in the HT and that are neighbors
    %to the list of neighbors.
    %fprintf('i: %d \t j: %d \t l: %d\n',ind_i,ind_j,l);
    % disp(List_Neighbors_id)
    %Key calulation
    cell_neighbor1 = create_cell(ind_i,ind_j,l);
    % disp(cell_neighbor1);
    key = HTObj.key_Calculation(cell_neighbor1);
    p = HTObj.hast_function(key);

    % fprintf('key: %d \t p: %d \n',key,p);

    %If there is a lisked list in this position
    if size(HTObj.HT{p+1},2) ~= 0 && HTObj.HT{p+1}.head.counter ~= 0
        %True or False
        flag = HTObj.HT{p+1}.searchNode(key);

        %Check if the cell is in the HT
        if strcmp(flag, 'True')
            cell_neighbor2 = HTObj.HT{p+1}.get_cell_Nodo(key); 
            l2 = cell_neighbor2.l;
        end

        %Check if the cell is in the HT, then adds id. Or else keep looking
        if strcmp(flag, 'True') && l == l2
            %Add the cell to the list neighbors
            List_Neighbors_id = [List_Neighbors_id [cell_neighbor2.id]];
        else
            List_Neighbors_id = Neighbor_Next_level(HTObj,cell_neighbor1,List_Neighbors_id,Dir);
            List_Neighbors_id = Neighbor_previous_level(HTObj,cell_neighbor1,List_Neighbors_id);
        end 
    else
        % List_Neighbors_id = Neighbor_Next_level(HTObj,cell_neighbor1,List_Neighbors_id,Dir);
        List_Neighbors_id = Neighbor_previous_level(HTObj,cell_neighbor1,List_Neighbors_id);
    end
end