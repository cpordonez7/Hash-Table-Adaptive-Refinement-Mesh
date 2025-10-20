function [List_Neighbors_id] = Neighbor_previous_level(HTObj,cell_neighbor,List_Neighbors_id)
    %This function return the new neighbors cells of the previous level 

    l = cell_neighbor.l;

    if l > 0
        New_cell = calculation_neighbor_previous_level(HTObj,cell_neighbor);
        %Data of the thick cell
        ind_i = New_cell(1,1);
        ind_j = New_cell(2,1);
        l = New_cell(3,1);

        %Key calulation
        cell_neighbor = create_cell(ind_i,ind_j,l);
        key = HTObj.key_Calculation(cell_neighbor);
        p = HTObj.hast_function(key);

        %Check if this cell is in the list, where k is a boolean value.
        %0----> NOT
        %1----> YES
        k1 =isempty(List_Neighbors_id);  %logic value, is emppty the cell neighbors list
        k2=1; %logic value, is emppty the cell repetid list
        
        
        if size(HTObj.HT{p+1},2) ~= 0 && HTObj.HT{p+1}.head.counter ~= 0
            %True or False
            flag = HTObj.HT{p+1}.searchNode(key);

            %Check if the cell is in the HT, then adds i. Or else keep looking
            if strcmp(flag, 'True')
                cell_neighbor2 = HTObj.HT{p+1}.get_cell_Nodo(key);
                l2 = cell_neighbor2.l;

                %if the cell are same, so check if this cell is in the list
                if k1 ~= 1 && l==l2
                    position_cell_repeated = find(List_Neighbors_id == cell_neighbor2.id); %it´s a list
                    k2 = isempty(position_cell_repeated);
                end
             
                %If the list does not have the cell, so insert. Elseif
                %continue the search
                if k2 == 1 && l==l2
                    List_Neighbors_id = [List_Neighbors_id [cell_neighbor2.id]];
                elseif l > 0  && l~=l2
                    List_Neighbors_id = Neighbor_previous_level(HTObj,cell_neighbor,List_Neighbors_id);
                end

            elseif l > 0
                List_Neighbors_id = Neighbor_previous_level(HTObj,cell_neighbor,List_Neighbors_id);
            end

        else
        List_Neighbors_id = Neighbor_previous_level(HTObj,cell_neighbor,List_Neighbors_id);
        end
    end
end