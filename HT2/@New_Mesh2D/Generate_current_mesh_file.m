function [] = Generate_current_mesh_file(HTObj)
    %This function create a file with the mesh that contains the object
    %This file is generates acording to the iteration in each linken list 
    
    HTObj.set_Number_cells_Mesh;   %Updata the cells number
    m = HTObj.Number_cells_Mesh;   %Cells number in the mesh
    
    fileID = fopen('File_current_mesh.txt','w');
    fprintf(fileID,'%d \n',m);   % First data: Number of point from grid

    %Difference domial extrem
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;


    %loop for HT
    for p=1:HTObj.size_HT
        
        if size(HTObj.HT{p},2) ~= 0 && HTObj.HT{p}.head.counter ~= 0
            HTObj.HT{p}.current = HTObj.HT{p}.first_node;
            while (HTObj.HT{p}.current ~= HTObj.HT{p}.last_node)
                %date from cell in the p linked list 
                ind_i = HTObj.HT{p}.current.cell_struct.i;
                ind_j = HTObj.HT{p}.current.cell_struct.j;
                l = HTObj.HT{p}.current.cell_struct.l;
                
                %step size
                h = (Dif_dom_x)/(HTObj.N*2^(l));
                k = (Dif_dom_y)/(HTObj.M*2^(l));

                % fprintf('Dif_dom_X: %f\t Dif_dom_Y: %f\t h: %.5f \t k: %.5f \n',Dif_dom_x, Dif_dom_y, h, k);

                x = HTObj.a1 + ind_i*h;
                y = HTObj.a2 + ind_j*k;    

                %Write the coordinates of each cell in the file
                fprintf(fileID,'%.12f \t %.12f \t %d\n',x,y,l);

                HTObj.HT{p}.current = HTObj.HT{p}.current.next;
            end

            if (HTObj.HT{p}.current == HTObj.HT{p}.last_node)
                %date fron last cell in the p linked list 
                ind_i = HTObj.HT{p}.current.cell_struct.i;
                ind_j = HTObj.HT{p}.current.cell_struct.j;
                l = HTObj.HT{p}.current.cell_struct.l;
                
                %step size
                h = (Dif_dom_x)/(HTObj.N*2^(l));
                k = (Dif_dom_y)/(HTObj.M*2^(l));

                x = HTObj.a1 + ind_i*h;
                y = HTObj.a2 + ind_j*k;

                %Write the coordinates of each cell in the file
                fprintf(fileID,'%.12f \t %.12f \t %d\n',x,y,l);
            end
        end
    end

end 
    
      