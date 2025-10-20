function [] = Refine_Mesh(HTObj)
    %his function refine all mesh.

    % Enume_Vect = get_enumeration_vector(HTObj);
    numcell = length(HTObj.Enumeration_vector);

    %Loop for Enumeration vector
    % for j=1:length(Enume_Vect)
    for j=1:numcell
        cell_to_ref = HTObj.getCell_id(j); 
        HTObj.RefineCell(cell_to_ref);
    end
end