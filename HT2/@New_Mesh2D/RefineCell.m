function [] = RefineCell(HTObj,cell_struct)
    %This function refine a cell of the mesh, 
    %input: cell struct that you want to refine

    l = cell_struct.l;
    Id = cell_struct.id;
    % key = cell_struct.key;
    

    %indixes of the new cell left-below
    ind_i = 2*cell_struct.i;
    ind_j = 2*cell_struct.j;

    set_Deletecell(HTObj,cell_struct);

    %Update the Mesh, add the new nodes

    %FIRT CELL
    cell_1 = create_cell(ind_i,ind_j,l+1);
    set_Insertcell(HTObj,cell_1);     %Insert the cell
    key_cel1 = HTObj.key_Calculation(cell_1);
    p = HTObj.hast_function(key_cel1); 
    % fprintf('cel1 key:%d \t p:%d \n',key_cel1,p);

    %Update id of the new cell, leaving the Id of the thick cell 
    HTObj.HT{p+1}.Update_node_structure_id(key_cel1,Id); 

    %Subtract one because inserting adds (elimina una porque inserta una nueva, entonces el contador se aumenta en una)
    HTObj.Enumeration_counter = HTObj.Enumeration_counter - 1 ;  

    %SECOND CELL
    cell_2 = create_cell(ind_i+1,ind_j,l+1);
    set_Insertcell(HTObj,cell_2);     %Insert the cell
    % key_cel2 = HTObj.key_Calculation(cell_2);
    % p = HTObj.hast_function(key_cel2); 
    % % fprintf('cel2 key:%d \t p:%d \n',key_cel2,p);
    % disp(HTObj.HT{199})



    %tHIRTD CELL
    cell_3 = create_cell(ind_i,ind_j+1,l+1);
    set_Insertcell(HTObj,cell_3);     %Insert the cell
    % key_cel3 = HTObj.key_Calculation(cell_3);
    % p = HTObj.hast_function(key_cel3); 
    % fprintf('cel3 key:%d \t p:%d \n',key_cel3,p);

    %FOURD CELL
    cell_4 = create_cell(ind_i+1,ind_j+1,l+1);
    set_Insertcell(HTObj,cell_4);     %Insert the cell
    % key_cel4 = HTObj.key_Calculation(cell_4);
    % p = HTObj.hast_function(key_cel4); 
    % fprintf('cel4 key:%d \t p:%d \n',key_cel4,p);


    %Updates all propietates of the object
    set_Number_cells_Mesh(HTObj);

end
      
