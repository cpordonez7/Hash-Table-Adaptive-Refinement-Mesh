function [List_New_cell] = calculation_neighbor_previous_level(HTObj,cell_neighbor)
    %Retorna la nueva celula en el nivel mas grueso 

    %Calculo de los indices de la celula
    ind_i = cell_neighbor.i;
    ind_j = cell_neighbor.j;
    l = cell_neighbor.l;

    %according the values of the i and j find the new cell
    if mod(ind_i,2) == 0
        ind_i_new = ind_i/2;
    else
        ind_i_new = (ind_i-1)/2;
    end

    if mod(ind_j,2) == 0
        ind_j_new = ind_j/2;
    else
        ind_j_new = (ind_j-1)/2;
    end

    %Celula gruesa
    List_New_cell = [ind_i_new;ind_j_new;l-1];

end
