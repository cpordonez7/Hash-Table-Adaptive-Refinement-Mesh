function [Estruct_Cell_Neihbors] = get_Cell_Neighbors(HTObj, id)
    %Return a structure with the vertices of neighbors cells and a list
    %with the quantity of neightbors by direction.
    %In the list of neighbors only are the neighbors actives, no including
    %the borders cell.

    % %Get the vector with the cell's directions
    % Vec_enum = HTObj.get_enumeration_vector();

    %Difference domial extrem
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;

    %getVecinas_estencil return (id, Dif_level)
    Lista_Id_cell_stencil = HTObj.getNeighbors_stencil(id);
    
    size_cell_stencil = nnz(Lista_Id_cell_stencil(1,:));  %nnz: Calculate the number of elements other than zero

    Lista_Cell_Neighbors = zeros(3,size_cell_stencil);

    %Counter the neighbors for direction
    Num_Neighbors = zeros(4,1);
    %Counter loop While
    iter = 1;
    iter_j = 1;

    %Counter for filling the list of neighbors, only the active ones enter
    iter_Neighbor = 1;

    %Counter of cell for direction
    Conter_cell_direction = 1;   

    %Number of the neighbors in the stencil including borders
    Number_cell_neighbors_stencil = size(Lista_Id_cell_stencil,2);


    %While for filling the list of neighbors. This contains the coordinates
    %of the cell
    while iter_j <= Number_cell_neighbors_stencil
        Id_neighboar = Lista_Id_cell_stencil(1,iter_j);
        if Id_neighboar == 0
            iter_j = iter_j + 1;
        else            
            %call the values of (i,j,l) of the cell with Id_neighboar enumeration
            % p = Vec_enum{Id_neighboar}(1);
            % key = Vec_enum{Id_neighboar}(2);
            % 
            % %Return the struct of cell with this key
            % cell_id = HTObj.HT{p}.get_cell_Nodo(key);
            cell_id = HTObj.getCell_id(Id_neighboar);              %input id cell

            %Indixes of the cell 
            ind_i = cell_id.i;
            ind_j = cell_id.j;
            l = cell_id.l;
            N_level = HTObj.N*2^(l);
            M_level = HTObj.M*2^(l);

            %step size
            dx = (Dif_dom_x)/(N_level);
            dy = (Dif_dom_y)/(M_level);

            x = HTObj.a1 + ind_i*dx;
            y = HTObj.a2 + ind_j*dy; 

            Lista_Cell_Neighbors(1:3,iter_Neighbor) = [x;y;l];
            iter_Neighbor = iter_Neighbor + 1;
            iter_j = iter_j + 1;
        end

    end

    %While to count neighboring cells by direction 
    while iter <= Number_cell_neighbors_stencil
        Id_neighboar = Lista_Id_cell_stencil(1,iter);
        dif_level = Lista_Id_cell_stencil(2,iter); 

        if dif_level > 0 && Id_neighboar ~= 0
            Num_Neighbors(Conter_cell_direction) = 2^dif_level;
            iter = iter + 2^dif_level;

        elseif dif_level <= 0 && Id_neighboar ~= 0
            Num_Neighbors(Conter_cell_direction) = 1;
            iter =iter +1;
        else
            Num_Neighbors(Conter_cell_direction) = 0;
            iter =iter +1;
        end
        Conter_cell_direction = Conter_cell_direction + 1;
    end

    Estruct_Cell_Neihbors.Cell_Neighbors = Lista_Cell_Neighbors;
    Estruct_Cell_Neihbors.Counter_Neighbors_direction = Num_Neighbors;

end