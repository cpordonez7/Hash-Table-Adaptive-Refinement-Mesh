function [] = Generate_current_mesh_file_creation(HTObj)
    %This function create a file with the mesh that contains the object
    %This file is generates acording to the creation of the mesh 
    
    HTObj.set_Number_cells_Mesh;   %Updata the cells number
    m = HTObj.Number_cells_Mesh;   %Cells number in the mesh

    HTObj.set_enumeration_vector;
    Vec_enum = HTObj.get_enumeration_vector;   %This vector has the p and key de each cell. 
    % Each position of the vector has (p,key) that are the data of a cell in the HT
   
    Num_incog = length(Vec_enum);           %Number of cell
    
    fileID = fopen('File_current_mesh_creation.txt','w');
    fprintf(fileID,'%d \n',m);   % First data: Number of point from grid

    %Difference domial extrem
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;

    %Loop walks on the all cell in the vector
    for j=1:(Num_incog)
    % for j=11:11
        %Data of the cell
        Id = j;           %Id is the vector's position
        % disp(j)
        %calculate index (i,j) the Id cell (F2 cell)
        cell_id = HTObj.getCell_id(Id);
        index_i= cell_id.i;
        index_j= cell_id.j;
        l = cell_id.l;
        
        %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
        h_c = (Dif_dom_x)/(HTObj.N*2^l);
        k_c = (Dif_dom_y)/(HTObj.M*2^l);

        %coordinates
        x = HTObj.a1 + index_i*h_c;
        y = HTObj.a2 + index_j*k_c;

        %Write the coordinates of each cell in the file
        fprintf(fileID,'%.12f \t %.12f \t %d\n',x,y,l);

    end

end 
    
      