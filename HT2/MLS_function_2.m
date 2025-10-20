function[A,b] = MLS_function_2(T,ObjHT,A,b,ID_CELL_CENT,Id_cel_vecina,pesoh,pesok,m,dl)
    %This function do the calculate the coefficientes the Matrix the system
    %for aplication of the DF
    %MLS_function_2 is the function with other case for the choose points
    %Calculo de distancias
    %Difference domial extrem
    Dif_dom_x = T.b1 - T.a1;
    Dif_dom_y = T.b2 - T.a2;

    ID_CELL_CENTRAL = ID_CELL_CENT;                    %ID central cell, Stencil DF
    % disp(ID_CELL_CENTRAL);
    %calculate index (i,j) the Id cell (F2 cell)
    cell_id = ObjHT.getCell_id(ID_CELL_CENTRAL);
    index_i = cell_id.i;
    index_j = cell_id.j;
    l_CELL_CEN = cell_id.l;
    % fprintf('id: %d \t i: %d \t j: %d \t l: %d \n',ID_CELL_CENTRAL,index_i,index_j,l_CELL_CEN);

    % %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
    h_cel_cen = (Dif_dom_x)/(T.N*2^l_CELL_CEN);
    k_cel_cen = (Dif_dom_y)/(T.M*2^l_CELL_CEN);

    x_CELL_CEN = T.a1 + index_i*h_cel_cen;
    y_CELL_CEN = T.a2 + index_j*k_cel_cen;
        
    %Coordinates of the center of the cell
    x_CELL_CEN = x_CELL_CEN + h_cel_cen*0.5;   % x value of the curreent cell (center)
    y_CELL_CEN = y_CELL_CEN + k_cel_cen*0.5;   % y value of the curreent cell (center)

    %Auxiliary lists
    P1 = zeros(2,1);

    %coordinates of non-existent cell in stencil
    if dl ~= 0    %if the level defference is no zero, the non-existen cell
        m_dir = m;  %cell direction

        %coordenates of the non-existen cell 
        if m_dir == 1
            P1 = [x_CELL_CEN;y_CELL_CEN-k_cel_cen];
        elseif m_dir == 2
            P1 = [x_CELL_CEN-h_cel_cen;y_CELL_CEN];
        elseif m_dir == 3
            P1 = [x_CELL_CEN+h_cel_cen;y_CELL_CEN];
        else
            P1 = [x_CELL_CEN;y_CELL_CEN+k_cel_cen];
            %P1 = [P1 [x_CELL_CEN;y_CELL_CEN+k_cel_cen]];
        end
    end

    %Calculation Weights 
    %Data of the non-existent cell
    Id_3 = Id_cel_vecina;  %Id existen cell in the stencil in level different  
    Vecinos_Id3 = ObjHT.getNeighbors_edges_vertices_id(Id_3);  %Neighbors of cells of different levels, only Id
    % disp(Id_3)
    %%Se actializa Vecinos_Id3 porque se agrego una celula en la 1 posicion
    Vecinos_Id3 = [Id_3 Vecinos_Id3];
    % disp(Vecinos_Id3)

    %This is the list with neighbors of neighbors not repeated
    Vecinos_Id3 = Add_cell_neighbors(ObjHT,Vecinos_Id3);
    % disp(Vecinos_Id3)

    n1 = size(Vecinos_Id3,2);
    % P_list = [];
    P_list = zeros(2,n1);   %P_list has Id3 cell and her neighbor

    %Coordenates of the neighbors cell 
    %n1+1 because the list has a Id_3 cell in the first positions
    for iter_k=1:n1
        Id_ = Vecinos_Id3(1,iter_k);
        %Id_ = Vecinos_Id3(iter_k);
        %calculate index (i,j) the Id 
        cell_id_ = ObjHT.getCell_id(Id_);
        ind_i_= cell_id_.i;
        ind_j_= cell_id_.j;
        l = cell_id_.l;
        
        %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
        h_c = (Dif_dom_x)/(T.N*2^l);
        k_c = (Dif_dom_y)/(T.M*2^l);

        coordx = T.a1 + ind_i_*h_c;
        coordy = T.a2 + ind_j_*k_c;
        
        P_list(:,iter_k) = [coordx + h_c*0.5;coordy + k_c*0.5];  %Points list
    end 


    [Point_list , Vecinos_Id3] = choose_Neighbors(Vecinos_Id3,P_list,P1);
    % disp(Vecinos_Id3);
    
    % disp(P1)
    % disp(length(Point_list))
    PESOS_W = Weights_wi(P1,Point_list);
    % disp(PESOS_W)

    %Here the matrix is created 
    %in each position of the present cell, write its coefficient, that is
    %the multiplication of the weights Pesos * w_i
    
    n2 = length(PESOS_W);
  
    if (m==1 || m==4)   % B and A cell
        for iter_j=1:n2
            A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j)) = A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j))  + pesok*PESOS_W(iter_j);
        end
    else    %L and R cell
        for iter_j=1:n2
            A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j)) = A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j))  + pesoh*PESOS_W(iter_j);
        end
    end
    % fprintf('%.5f \t  %.5f \t \t %.5f\n',P1(1,iter_j), P1(2,iter_j), Aprox(iter_j));
end

function [Vecinos_Id3] = Add_cell_neighbors(ObjHT,Vecinos_Id3)
    Add_cell_1 = [];
    Add_cell = [];

    %calculation of auxiliary neighbors (Neighbors of neighbors)
    for iter=1:size(Vecinos_Id3,2)
        Add_cell_1 = [Add_cell_1 ObjHT.getNeighbors_id(Vecinos_Id3(:,iter))];
    end

    %Delete repite cell
    n_aux_1 = length(Add_cell_1);
    for iter_i = 1:n_aux_1
        posiciones_cell_repetidas = find(Vecinos_Id3(1,:) == Add_cell_1(iter_i));
        vlr_log = isempty(posiciones_cell_repetidas);

        %If the cell is not in the neighbors add the list  (0:no, 1:si)
        if vlr_log==1
            Vecinos_Id3 = [Vecinos_Id3 Add_cell_1(iter_i)];
            % Add_cell = [Add_cell Add_cell_1(iter_i)];
            % Vecinos_Id3 = [Vecinos_Id3 Vecinos_Id3(1,iter_i)];
            % disp(Vecinos_Id3)
        end
    end

    %Coordenates of auxiliary neighbors cell
    % % % n_aux = size(Add_cell,2);
    % % % Vecinos_Id3_auxiliar = zeros(1,n_aux);
    % % % for iter_k=1:n_aux
    % % %     Id_cel_aux = Add_cell(iter_k);
    % % % 
    % % %     %calculate index (i,j) the Id 
    % % %     % cell_id_aux = ObjHT.getCell_id(Id_cel_aux);
    % % %     % l = cell_id_aux.l;
    % % % 
    % % %     % Vecinos_Id3_auxiliar(:,iter_k) = [Id_cel_aux;l];  %Coloco l para llenar porque es dif_nivel
    % % %     Vecinos_Id3_auxiliar(1,iter_k) = Id_cel_aux;  %Coloco l para llenar porque es dif_nivel
    % % % end 
    % % % 
    % % % Vecinos_Id3 = [Vecinos_Id3 Vecinos_Id3_auxiliar];
    % debe adicionar a la lista de vecinas de Id_3 la nueva celula

end

function [Point_list, Vecinos_Id3_update] = choose_Neighbors(Vecinos_Id3,P_list,P1)
    N_point_choose = 16;
    n_point_Neighbors = size(Vecinos_Id3,2);
    % disp(n_point_Neighbors);

    if n_point_Neighbors < N_point_choose
        N_point_choose = n_point_Neighbors;
        % disp(n_point_Neighbors);
    end
    
    Point_list = zeros(2,N_point_choose);
    Vecinos_Id3_update = zeros(1,N_point_choose);

     %Distance calculation 
    n_aux_2= size(P_list,2);
    Distance_list = zeros(1,n_aux_2);
    for iter_j=1:n_aux_2
        Distance_list(iter_j) = norm(P_list(:,iter_j)-P1,2);
    end
    % disp(Distance_list)
    
    % Ordenar la lista y obtener los índices originales
    [Ordered_distances, index_aux] = sort(Distance_list);
    % disp(Ordered_distances)
    % disp(index_aux)

    for iter_i=1:N_point_choose
        Vecinos_Id3_update(1,iter_i) = Vecinos_Id3(1,index_aux(iter_i)); 
        Point_list(:,iter_i) = P_list(:,index_aux(iter_i));
    end
end 



