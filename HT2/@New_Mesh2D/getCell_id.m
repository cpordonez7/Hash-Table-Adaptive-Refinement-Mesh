function [Cell_struct] = getCell_id(HTObj, id)
    %Return the struct of cell with this id

    %Get the vector with the cell's directions
    % ID_vector = HTObj.get_enumeration_vector;    
    % 
    % p = ID_vector{id}(1);
    % key = ID_vector{id}(2);

    p = HTObj.Enumeration_vector{id}(1);
    key = HTObj.Enumeration_vector{id}(2);

    %Return the struct of cell with this key
    % cell_id = HTObj.HT{p}.get_cell_Nodo(key);

    Cell_struct = HTObj.HT{p}.get_cell_Nodo(key); 
end