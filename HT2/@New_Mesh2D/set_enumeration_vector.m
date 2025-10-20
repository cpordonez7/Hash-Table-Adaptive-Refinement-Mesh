function [] = set_enumeration_vector(HTObj)
    %This function updata the enumeration vector than contain the p and key of the cell
    %according to its enumeration
    
    HTObj.set_Number_cells_Mesh();
    HTObj.Enumeration_vector = cell(HTObj.Number_cells_Mesh,1); %Cells Vector V{k} = [p key].

    %loop for HT
    for p=1:HTObj.size_HT
        
        if size(HTObj.HT{p},2) ~= 0 && HTObj.HT{p}.head.counter ~= 0
            HTObj.HT{p}.current = HTObj.HT{p}.first_node;
            while (HTObj.HT{p}.current ~= HTObj.HT{p}.last_node)
                %date fron cell in the p linked list 
                key = HTObj.HT{p}.current.cell_struct.key;
                id = HTObj.HT{p}.current.cell_struct.id;
                % disp(id)
                HTObj.Enumeration_vector{id} = [p key]; 

                HTObj.HT{p}.current = HTObj.HT{p}.current.next;
            end

            if (HTObj.HT{p}.current == HTObj.HT{p}.last_node)
                key = HTObj.HT{p}.current.cell_struct.key;
                id = HTObj.HT{p}.current.cell_struct.id;
                % disp(id)
                HTObj.Enumeration_vector{id} = [p key];
            end
        end
    end
end