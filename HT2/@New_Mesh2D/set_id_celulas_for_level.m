function [] = set_id_celulas_for_level(HTObj)
    %This function updata the id of the cells for level 

    %vector with the cell number for level
    cell_for_level = HTObj.Cell_counter_level;

    %Counter for id
    temporary_cell_counter = 1;

    %iterator to fill array at each level
    iter_level_matriz = ones(1,HTObj.Max_level+1);

    %Create the lists for level type cells
     list_cell_level = cell(HTObj.Max_level+1,1);

     %Created in each posotion a vector
     for j=1:HTObj.Max_level+1
         list_cell_level{j} = zeros(2,cell_for_level(j));  
     end
   %Loop HT
   for p=1:HTObj.size_HT
       if size(HTObj.HT{p},2) ~= 0 && HTObj.HT{p}.head.counter ~= 0
           HTObj.HT{p}.current = HTObj.HT{p}.first_node;
           while (HTObj.HT{p}.current ~= HTObj.HT{p}.last_node)
               %date fron cell in the p linked list 
               i = HTObj.HT{p}.current.cell_struct.i;
               j = HTObj.HT{p}.current.cell_struct.j;
               l = HTObj.HT{p}.current.cell_struct.l;
               list_cell_level{l+1}(:,iter_level_matriz(l+1)) = [i; j];
               HTObj.HT{p}.current = HTObj.HT{p}.current.next;
               iter_level_matriz(l+1) = iter_level_matriz(l+1) +1;
           end

           if (HTObj.HT{p}.current == HTObj.HT{p}.last_node)
               i = HTObj.HT{p}.current.cell_struct.i;
               j = HTObj.HT{p}.current.cell_struct.j;
               l = HTObj.HT{p}.current.cell_struct.l;
               list_cell_level{l+1}(:,iter_level_matriz(l+1)) = [i; j];
               iter_level_matriz(l+1) = iter_level_matriz(l+1) +1;
           end
        end
   end

    %sort the cells 
    for p=1:HTObj.Max_level+1
        num_cell_for_level = cell_for_level(p);
        % disp(list_cell_level{p})

        arrtr = list_cell_level{p}';
        %Transpone y ordena para i
        [arr_1, indices1] = sortrows(arrtr, 2);
        %Ordena para j
        [arr_2, indices2] = sortrows(arr_1, 1);

        %Ordena para i
        [arr_3, indices3] = sortrows(arr_2, 2);
        % disp(arr_3)
        arr_3 =arr_3';
        

        for iter=1:num_cell_for_level
            i = arr_3(1,iter);
            j = arr_3(2,iter);
            cell_ = create_cell(i,j,p-1);
            key_cell = HTObj.key_Calculation(cell_);
            p_ = HTObj.hast_function(key_cell); 

            %Update id of the cell 
            HTObj.HT{p_+1}.Update_node_structure_id(key_cell,temporary_cell_counter);   
            temporary_cell_counter = temporary_cell_counter + 1;

        end
    end

end