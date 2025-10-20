function [List_Neighbors_id] = getNeighbors_direction_id(HTObj, id, Dir)
    %Return the id of the neighbors cells to the referent cell to id

    List_Neighbors_id = [];

    %call the values of (i,j,l) of the cell with id enumeration
    cell_id = HTObj.getCell_id(id);           %input id cell

    %Indixes of the cell 
    ind_i = cell_id.i;
    ind_j = cell_id.j;
    l = cell_id.l;
    N_level = HTObj.N*2^(l);
    M_level = HTObj.M*2^(l);

    % fprintf('i: %d \t j: %d \t l: %d\n',ind_i,ind_j,l);
    %Analyze all cases according to the position of the cell in the mesh

    %case 1 cell belown-left
    if (ind_i == 0) && (ind_j == 0)
        switch Dir
            case 'B'      %Cells below
                disp('There are not cell in this direction');                   

            case 'L'      %Cells left
                disp('There are not cell in this direction'); 

            case 'R'      %Cells right
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_id); %return only the id cells

            case 'A'      %Celulas above
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_id);

            otherwise
                disp('Error, You entered an incorrect address.')  
        end
        
        
    end

    %Case 2 cell belown
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == 0)
        switch Dir
            case 'B'      %Cells below
                disp('There are not cell in this direction');                   

            case 'L'      %Cell left
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_id);

            case 'R'      %Cells right
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_id);

            case 'A'      %Cells above
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_id);

            otherwise
                disp('Error, You entered an incorrect address.')  
        end
    end

    %Case 3 cell rigth belown
    if (ind_i == N_level-1) && (ind_j == 0)
        switch Dir
            case 'B'      %Cells below
                disp('There are not cell in this direction');                   

            case 'L'      %Cells left
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_id);

            case 'R'      %Cells right
                disp('There are not cell in this direction');  

            case 'A'      %Cells above
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_id);

            otherwise
                disp('Error, You entered an incorrect address.')  
        end    
        
    end

    %Case 4 cell left
    if (ind_i == 0) && (0 < ind_j && ind_j < M_level-1)
        switch Dir
            case 'B'      %Cells below
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_id);

            case 'L'      %Cells left
                disp('There are not cell in this direction');  

            case 'R'      %Cells right
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_id);

            case 'A'      %Cells above
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_id);

            otherwise
                disp('Error, You entered an incorrect address.')   
        end    
        
    end

    %Case 5 cell rigth
    if (ind_i == N_level-1) && (0 < ind_j && ind_j < M_level-1)
        switch Dir
            case 'B'      %Cells below
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_id);

            case 'L'      %Cells left
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_id);

            case 'R'      %Cells right
                disp('There are not cell in this direction');  

            case 'A'      %Cells above
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_id);

            otherwise
                disp('Error, You entered an incorrect address.')   
        end
     
    end

    %Case 6, cell left above
    if (ind_i == 0) && (ind_j == M_level-1)
        switch Dir
            case 'B'      %Cells below
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_id);

            case 'L'      %Cells left
                disp('There are not cell in this direction'); 

            case 'R'      %Cells right
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_id);

            case 'A'      %Cells above
                disp('There are not cell in this direction');

            otherwise
                disp('Error, You entered an incorrect address.') 
        end
   
    end

    %Case 7, cell above
    if (0 < ind_i && ind_i < N_level-1) && (ind_j == M_level-1)
        switch Dir
            case 'B'      %Cells below
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_id);

            case 'L'      %Cells left
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_id);

            case 'R'      %Cells right
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_id);

            case 'A'      %Cells above
                disp('There are not cell in this direction');

            otherwise
                disp('Error, You entered an incorrect address.')  
        end        
        
    end

    %Case 8, cell rigth above
    if (ind_i == N_level-1) && (ind_j == M_level-1)
        switch Dir
            case 'B'      %Cells below
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_id);

            case 'L'      %Cells left
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_id);

            case 'R'      %Cells right
                disp('There are not cell in this direction');

            case 'A'      %Cells above
                disp('There are not cell in this direction');

            otherwise
                disp('Error, You entered an incorrect address.')   
        end    
        
    end

    %Case 9, cell center
    if (0 < ind_i && ind_i < N_level-1) && (0 < ind_j && ind_j < M_level-1)
        switch Dir
            case 'B'      %Cells below
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j-1,l,'B',List_Neighbors_id);

            case 'L'      %Cells left
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i-1,ind_j,l,'L',List_Neighbors_id);

            case 'R'      %Cells right
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i+1,ind_j,l,'R',List_Neighbors_id);

            case 'A'      %Cells above
                List_Neighbors_id = Insert_neighbor(HTObj,ind_i,ind_j+1,l,'A',List_Neighbors_id);

            otherwise
                disp('Error, You entered an incorrect address.') 
        end 
        
    end

end