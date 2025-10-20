classdef New_Mesh2D < handle
%This class creates an object Hash Table (HT).
%This class is a data structure that allows us to save and search elements
%by middle of a key. The data stored here is of type structure.


%-----------------------------------------------------------------
%               PROPERTIES SECTION
%-----------------------------------------------------------------
    properties (Access = private)
        data_file;              %Data file. This file saves the mesh nodes 
        size_HT;                %Hast Table size
        HT;                     %Auxiliary property to call the object
        Max_level;              %Hast Table max level     
        Mesh_nodes;             %Array that saves the file data. This contain the mesh nodes 

        Number_cells_Mesh;        %Number of cells in the Mesh (Hast Table)        
        Cell_counter_level;    %Enumeration, cell counter vector for level
        Enumeration_counter;   %Enumeration counter of the mesh
        Enumeration_vector;     %Enumeration vecctor. In each positions has [p,key]

        %Domain data 
        N;                      %Number of splits on the axe x. Number of base mesh cells
        M;                      %Number of splits on the axe y. Number of base mesh cells
        a1;                     %Domain lower end, exe x 
        b1;                     %Domain top end, exe x
        a2;                     %Domain lower end, exe y 
        b2;                     %Domain top end, exe y


    end

%-----------------------------------------------------------------
%               METHODS SECTION
%-----------------------------------------------------------------
    methods
        function [HTObj] = New_Mesh2D(mesh_data_file, max_level,grid)
            %Constructor: Initializes the values of the properties with the
            %input parameters and creates the HT.
            %
            %Input parameters:
            %   mesh_data_file - Defines the data file that contains the
            %   nodes of the mesh max_level - Defines the max mesh level(L).
            %   Base level equal to zero 
            %   grid - Structur with number of
            %   partitions on the axes and ends of the domain.
            %   grid.N-> partitions axe x
            %   grid.M-> partitions axe y
            %   grid.a1-> x-axis bottom end
            %   grid.b1-> x-axis top end
            %   grid.a2-> y-axis bottom end
            %   grid.b2-> y-axis top end
            %
            %Data file details (input parameter):
            %   mesh_data_file - This file is an ascci type file whit the 
            %   of the mesh information. This file has the following
            %   structure:
            %    
            %   number of nodes     First data of the file
            %   x1 y1 l             Node coordinates and mesh level
            %
            %   ...
            %--------------------------------------------------------------

             HTObj.data_file = mesh_data_file;
             file= fopen(mesh_data_file,'r');
            %-------------READ THE FILE------------------------------------
             numElem = fscanf(file,'%d \n',[1,1]);   %Number of nodes, first data of the file 
             % disp(numElem)
             HTObj.Mesh_nodes = fscanf(file,'%e %e %d\n',[3 numElem]);  %Array 3xNumElem, save (x,y,l) 
             % disp(HTObj.Mesh_nodes(:,:))
            %Domain data, N,a y b 
             HTObj.N = grid.N;
             HTObj.M = grid.M;
             HTObj.a1 = grid.a1;
             HTObj.b1 = grid.b1;
             HTObj.a2 = grid.a2;
             HTObj.b2 = grid.b2;

             HTObj.Max_level = max_level;    %Maximum level (L: 0, ... L-1)

            %Enumeration counter vector for level
            HTObj.Cell_counter_level = zeros(max_level+1,1);

            %Enumeration counter of the mesh
            HTObj.Enumeration_counter = 0;

            %Hast Table size
             number_cells_max_level = HTObj.N*2^(HTObj.Max_level)*HTObj.M*2^(HTObj.Max_level);  %Number of cells at max level
             HTObj.size_HT = round(0.10*number_cells_max_level);  

            %Initialize the HT
             HTObj.HT = cell(HTObj.size_HT,1);          %Create the cell vector, with empty cells in each position 

            % Filling of the HT: Each point of the file corresponds to the coordinates of the node of a cell,
            % we insert that cell in the HT
             for j=1:numElem
                x = HTObj.Mesh_nodes(1,j);          %x coordinate of the node         
                y = HTObj.Mesh_nodes(2,j);          %y coordinate of the node
                l = HTObj.Mesh_nodes(3,j);          %l level of the node

                %Index mapping (i,j).
                ind_i = HTObj.N*2^(l)*(x-HTObj.a1)/(HTObj.b1-HTObj.a1);
                ind_j = HTObj.M*2^(l)*(y-HTObj.a2)/(HTObj.b2-HTObj.a2);

                %Round
                ind_i = round(ind_i);
                ind_j = round(ind_j);

                %fprintf('i: %d\t j: %d\t l:%d \n',ind_j,ind_j,l);

                cell_str = create_cell(ind_i,ind_j,l);
                %disp(cell_str)

                set_Insertcell(HTObj,cell_str);     %Insert the cell
             end     

            %Count the cell number
             set_Number_cells_Mesh(HTObj);

            %Enumeration counter
             set_enumeration_vector(HTObj)

            % Update id for each cell for level
             % set_id_celulas_for_level(HTObj);

            % %Update id for each cell, according to the reading of the file
            % % set_id_celulas_file(HTObj,numElem); 

        end

%-----------------------------------------------------------------
%               PUBLIC INTERFACES SECTION
%-----------------------------------------------------------------
        [HT] = get_HT(HTObj);                           %Return the HT
        [Enumeration_counter] = get_Enumeration_counter(HTObj);   %This function return number of cell that the HT have in this moment
        [key] = key_Calculation(HTObj,cell_str)             %Cell key calculation 
        [key] = key_Calculation_i_j_l(HTObj,ind_i,ind_j,l)  %Cell key calculation with index and level
        [p] = hast_function(HTObj,key)                      %position of the cell in the HT.
        [] = set_Insertcell(HTObj,cell_str);                %Insert a cell in the HT
        [] = set_Deletecell(HTObj,cell_str);                %Delete a cell in the HT
        [Cell_struct] = getCell_id(HTObj, id)               %Return the struct of the cell with this id
        [] = set_enumeration_vector(HTObj);                 %Updata the enumeration vector
        [Enumeration_counter] = get_enumeration_vector(HTObj); %Return the enumeration vector
        [Size_HT] = get_size_HT(HTObj);                     %Return the size of the HT


        [] = graph_Mesh2D(HTObj);                           %Graph the mesh
        [] = graph_enumeration(HTObj);                      %Graph the enumartion of the cell 
        [] = Graph_Mesh2D_file(HTObj,file_graph)            %Graph the mesh of a file

        [Cell_counter_level] = getCell_counter_level(HTObj); %Enumeration counter vector for level
        [] = set_Number_cells_Mesh(HTObj)                   %keep the Number cell in the Mesh
        [Number_cells_Mesh] = get_Number_cells_Mesh(HTObj)                   %return the Number cell in the Mesh

        [List_Neighbors_id] = getNeighbors_id(HTObj, id)     %Id Neighbors, with input id
        %helper functions of the getVecinas function
        [List_Neighbors_id] = Insert_neighbor(HTObj,ind_i,ind_j,l,Dir,List_Neighbors_id);  
        [List_Neighbors_id] = Neighbor_Next_level(HTObj,cell_neighbor,List_Neighbors_id,Dir);
        [List_New_cells] = calculation_neighbor_next_level(HTObj,cell_neighbor,Dir);
        [List_Neighbors_id] = Neighbor_previous_level(HTObj,cell_neighbor,List_Neighbors_id);
        [List_New_cell] = calculation_neighbor_previous_level(HTObj,cell_neighbor);
        [List_Neighbor_below] = getNeighbors_below_id(HTObj,id);
        [List_Neighbors_above] = getNeighbors_above_id(HTObj,id);
        [List_Neighbors_left] = getNeighbors_left_id(HTObj,id);
        [List_Neighbors_right] = getNeighbors_right_id(HTObj,id);
        [List_Neighbors_stencil_RETURN] = getNeighbors_stencil(HTObj, id);
        [List_Neighbors_left_return] = getNeighbors_left_id_dif_level(HTObj,id);
        [List_Neighbors_below_return] = getNeighbors_below_id_dif_level(HTObj,id);
        [List_Neighbors_right_return] = getNeighbors_right_id_dif_level(HTObj,id);
        [List_Neighbors_above_return] = getNeighbors_above_id_dif_level(HTObj,id);
        [List_Neighbors_id] = getNeighbors_direction_id(HTObj, id, Dir);
        [List_Neighbors_id] = getNeighbors_edges_vertices_id(HTObj, id);
        [List_Neighbors_P] = getNeighbors_Punto(HTObj, P, lc);
        [List_Neighbors_quadr_interp, flag_direction, flag_aux] = getNeighbors_quadratic_interpolation(HTObj,Id,F1,F3,Dir,m);
        


        [] = RefineCell(HTObj,cell_struct);                   %Refine one cell of the mesh 
        [] = Refinar_region(HTObj,cell_struct1,cell_struct2);  %Refine a region of the mesh
        [] = Refine_Mesh(HTObj);                               %Refine all cell
        [] = ThickenCell(HTObj,cell_struct);                   %Thicken a cell
        [] = ThickenCell_v_sin_enum(HTObj,cell_struct);        %Thicken a cell sin sin actualizar la enumracion
        [] = Thinken_region(HTObj,cell_struct1,cell_struct2);  %Thicken a region of the mesh

        [collisions_vector] = get_number_collisions(HTObj);    %Vector of the collisions
        [] = Generate_current_mesh_file(HTObj);                %Create a file with the mesh in the object
        [] = Generate_current_mesh_file_creation(HTObj);       %Create a file with the mesh in the object for creation
        [] = set_initialize_Cell_counter_level(HTObj);         %Initialize the vector to zeros
        [] = set_id_celulas_for_level(HTObj);                  %Update the id for each cell for level
        
       

      
    end         %END METHODS SECTION



end         %END DEFINICION DE CLASE

                %EDN FILE