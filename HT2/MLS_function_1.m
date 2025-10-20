function[A,b, num_max_points] = MLS_function_1(T,ObjHT,A,b,ID_CELL_CENT,Id_cel_vecina,pesoh,pesok,m,dl)
    %This function do the calculate the coefficientes the Matrix the system
    %for aplication of the DF

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

    n1 = size(Vecinos_Id3,2);
    % P_list = zeros(2,n1+1);   %P_list has Id3 cell and her neighbor, therefor is n1+1 elementns
    P_list = [];

    %Coordenates of the neighbors cell 
    %n1+1 because the list has a Id_3 cell in the first positions
    for iter_k=1:n1+1
        if iter_k == 1
        Id_ = Id_3;
        else
        Id_ = Vecinos_Id3(1,iter_k-1);
        end
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

    %Conditions if P_list has less than 6 neighbors
    % disp(Vecinos_Id3)
    % disp(P_list)
    if size(P_list,2) < 6
        % disp('hay seis puntos minimo')
        [P_list , Vecinos_Id3] = Add_cell_neighbors(ObjHT,T,P_list,P1,Vecinos_Id3);
        % disp(Vecinos_Id3)
    end

    
    % disp(ID_CELL_CENTRAL)
    % disp(Vecinos_Id3)
    % disp(P1)
    % disp(P_list)
   
   
    num_max_points = size(P_list,2);
    % % disp(num_max_points)


    PESOS_W = Weights_wi(P1,P_list);
    % disp(PESOS_W)

    %%Se actializa Vecinos_Id3 porque se agrego una celula 

    %Here the matrix is created 
    %in each position of the present cell, write its coefficient, that is
    %the multiplication of the weights Pesos * w_i
    
    n2 = length(PESOS_W);
    % disp(n2)
  
    if (m==1 || m==4)   % B and A cell
        A(ID_CELL_CENTRAL,Id_3) = A(ID_CELL_CENTRAL,Id_3)  + pesok*PESOS_W(1);
        for iter_j=2:n2
            A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j-1)) = A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j-1))  + pesok*PESOS_W(iter_j);
        end
    else    %L and R cell
        A(ID_CELL_CENTRAL,Id_3) = A(ID_CELL_CENTRAL,Id_3)  + pesoh*PESOS_W(1);
        for iter_j=2:n2
            A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j-1)) = A(ID_CELL_CENTRAL,Vecinos_Id3(1,iter_j-1))  + pesoh*PESOS_W(iter_j);
        end
    end
    % fprintf('%.5f \t  %.5f \t \t %.5f\n',P1(1,iter_j), P1(2,iter_j), Aprox(iter_j));
end

function [P_list,Vecinos_Id3_update] = Add_cell_neighbors(ObjHT,T,P_list,P1,Vecinos_Id3)
    Add_cell_1 = [];
    Add_cell = [];

    %Difference domial extrem
    Dif_dom_x = T.b1 - T.a1;
    Dif_dom_y = T.b2 - T.a2;

    %calculation of auxiliary neighbors (Neighbors of neighbors)
    n_initial = size(Vecinos_Id3,2);
    for iter=1:n_initial
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
            Add_cell = [Add_cell Add_cell_1(iter_i)];
        end
    end

    %Coordenates of auxiliary neighbors cell
    n_aux = size(Vecinos_Id3,2)-n_initial;
    P_list_aux_cell = zeros(2,n_aux);
    Vecinos_Id3_auxiliar = zeros(1,n_aux);
    for iter_k=1:n_aux
        % Id_cel_aux = Add_cell(iter_k);
        Id_cel_aux = Vecinos_Id3(n_initial+iter_k);

        %calculate index (i,j) the Id 
        cell_id_aux = ObjHT.getCell_id(Id_cel_aux);
        ind_i_cel_aux= cell_id_aux.i;
        ind_j_cel_aux= cell_id_aux.j;
        l = cell_id_aux.l;
        
        %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
        h_c = (Dif_dom_x)/(T.N*2^l);
        k_c = (Dif_dom_y)/(T.M*2^l);

        coordx = T.a1 + ind_i_cel_aux*h_c;
        coordy = T.a2 + ind_j_cel_aux*k_c;

        P_list_aux_cell(:,iter_k)  = [coordx + h_c*0.5;coordy + k_c*0.5];  %Points list the auxiliary cell
        Vecinos_Id3_auxiliar(1,iter_k) = Id_cel_aux; 
    end 

    %Distance calculation 
    % % n_aux_2= size(P_list_aux_cell,2);
    % % Distance_list = zeros(1,n_aux_2);
    % % for iter_j=1:n_aux_2
    % %     Distance_list(iter_j) = norm(P_list_aux_cell(:,iter_j)-P1,2);
    % % end

    %Sort_list_distances = zeros(2,n_aux_2);
    %FALTA ORDENAR LA LISTA 

    %Extra cell
    Cell_faltantes = 6 - size(P_list,2);
    % disp(P_list_aux_cell)
    % LL = P_list_aux_cell(:,1:Cell_faltantes);
    % disp(LL)
    
    %Add auxiliary cell, complete the six points
    P_list = [P_list P_list_aux_cell(:,1:Cell_faltantes)];
    % Vecinos_Id3_aumentada = [Vecinos_Id3 Vecinos_Id3_auxiliar(1,1:Cell_faltantes)];
    % Vecinos_Id3_update = zeros(1,N_point_choose);
    Vecinos_Id3_update = [Vecinos_Id3 Vecinos_Id3_auxiliar(1,1:Cell_faltantes)];
    %debe adicionar a la lista de vecinas de Id_3 la nueva celula

end



