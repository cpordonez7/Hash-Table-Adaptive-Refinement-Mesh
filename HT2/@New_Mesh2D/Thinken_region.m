function [] = Thinken_region(HTObj,cell_struct1,cell_struct2)
    %This function Thinken a rectangular region of the mesh. 
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

    %In this process it fixes j and moves in i
    while j_1 <= j_2
        i_1 = i1;
        while i_1 <= i_2
            cell_to_thicken = create_cell(i_1,j_1,l1);    %temporal cell
            key_i = HTObj.key_Calculation(cell_to_thicken);
            p = HTObj.hast_function(key_i);

            cell_to_thicken = HTObj.HT{p+1}.get_cell_Nodo(key_i);
            HTObj.ThickenCell(cell_to_thicken);
            
            i_1 = i_1 + 2;
        end
        j_1 = j_1 + 2;
    end

end     
    
      