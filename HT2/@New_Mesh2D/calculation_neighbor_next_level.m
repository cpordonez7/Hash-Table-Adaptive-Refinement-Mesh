function  [List_New_cells] = calculation_neighbor_next_level(HTObj,cell_neighbor,Dir)
    %Find the cells of the next level, according to the indicated direction
    %considering the cases of diagonal cells

    l = cell_neighbor.l;

    %indixes cell left-below
    ind_i = 2*cell_neighbor.i;
    ind_j = 2*cell_neighbor.j;


    %Find the neighbors in this direction (Dir)
    switch Dir
        case 'B'      %Cells below
            List_New_cells = [ind_i ind_i+1; ind_j+1 ind_j+1; l+1 l+1];                     

        case 'L'      %Cells left
            List_New_cells = [ind_i+1 ind_i+1; ind_j ind_j+1; l+1 l+1];

        case 'R'      %Cells Righ
            List_New_cells = [ind_i ind_i; ind_j ind_j+1; l+1 l+1];

        case 'A'      %Cells Above
            List_New_cells = [ind_i ind_i+1; ind_j ind_j; l+1 l+1]; 

        case 'D1'      %Cell left-below
            List_New_cells = [ind_i+1; ind_j+1; l+1]; 

        case 'D2'      %Cell righ-below
            List_New_cells = [ind_i; ind_j+1; l+1]; 

        case 'D3'      %Cell left-above
            List_New_cells = [ind_i+1; ind_j; l+1]; 

        case 'D4'      %Cell righ-above
            List_New_cells = [ind_i; ind_j; l+1];
        otherwise
            disp('Error, ingresó una dirección inválida.')  
    end 

end
