function [] = Refine_region(HTObj,cell_struct1,cell_struct2)
    %This function refine a rectangular region of the mesh. 
    %Input: cell strcuctures of the lower left and upper right ends 
    % of the region.
    %This region has cells all of the same level.

    %cell left-below
    i1 = cell_struct1.i;
    j_1 = cell_struct1.j;
    l1 = cell_struct1.l;
    
    % key = cell_struct.key;

    %cell Right-above
    i_2 = cell_struct2.i;
    j_2 = cell_struct2.j;
    % key = cell_struct.key;

    %Proceso de refinamiento, fija y y se mueve en x
    while j_1 <= j_2
        i_1 = i1;
        while i_1 <= i_2
            cell_to_ref = create_cell(i_1,j_1,l1);    %temporal cell
            key_i = HTObj.key_Calculation(cell_to_ref);
            p = HTObj.hast_function(key_i);

            cell_to_ref = HTObj.HT{p+1}.get_cell_Nodo(key_i);
            HTObj.RefineCell(cell_to_ref);
            
            i_1 = i_1 + 1;
        end
        j_1 = j_1 + 1;
    end

end     
    
      