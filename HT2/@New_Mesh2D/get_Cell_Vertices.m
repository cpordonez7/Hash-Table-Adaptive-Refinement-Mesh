function [Lista_Cell_Vertices] = get_Cell_Vertices(HTObj, id)
    %This function return the four vertices of the cell with referent to id
    
    %Difference domial extrem
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;

    %call the values of (i,j,l) of the cell with id enumeration
    cell_id = HTObj.getCell_id(id);              %input id cell

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

    %Cell Vertices
    x1 = x + dx;
    y1 = y + dy;

    Lista_Cell_Vertices = [x x1 x1 x; y y y1 y1]; 

end