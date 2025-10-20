function [List_Neighbors_id] = Neighbor_Next_level(HTObj,cell_neighbor,List_Neighbors_id,Dir)
    %Esta funcion busca vecinos en un nivel mas fino

    % fprintf('i: %d \t j: %d \t l: %d\n',cell_neighbor.i,cell_neighbor.j,cell_neighbor.l);

    level_max_mesh = HTObj.Max_level;    
    % fprintf('%f \t %f \t %d',celula.x,celula.y,celula.nivel);

    l = cell_neighbor.l;
    % disp(l)
    
    if l < level_max_mesh-1
        New_cell = calculation_neighbor_next_level(HTObj,cell_neighbor,Dir);   %Return cells structs

        num_new_cell= size(New_cell,2);

        for j=1:num_new_cell
            % for j=1:2

            ind_i = New_cell(1,j);
            ind_j = New_cell(2,j);
            l = New_cell(3,j);
            
             %Key calulation
            cell_neighbor = create_cell(ind_i,ind_j,l);
            key = HTObj.key_Calculation(cell_neighbor);
            p = HTObj.hast_function(key);

            if size(HTObj.HT{p+1},2) ~= 0 && HTObj.HT{p+1}.head.counter ~= 0
                %True or False
                flag = HTObj.HT{p+1}.searchNode(key);

                if strcmp(flag, 'True')
                    cell_neighbor2 = HTObj.HT{p+1}.get_cell_Nodo(key); 
                    l2 = cell_neighbor2.l;
                end

                %Check if the cell is in the HT, then adds i. Or else keep looking
                if strcmp(flag, 'True') && l == l2
                    %Add the cell to the list neighbors
                    List_Neighbors_id = [List_Neighbors_id [cell_neighbor2.id]];
                else
                    List_Neighbors_id = Neighbor_Next_level(HTObj,cell_neighbor,List_Neighbors_id,Dir);
                end  
            end
        end
    end
 end 






