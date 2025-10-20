function [A,b] = discretizeDFPoisson_Inter_orden_uno(T,pde,ObjHT)
    %This function creates the matrix A and vector b (rhs), from the
    %discretization of FD, interpolación Lagrange orden uno Dirichlet
    % %input grid, pde and HT
    %Use the vector to calculate the id, then calculate neighbors

    Vec_enum = ObjHT.get_enumeration_vector();   %This vector has the direction de each cell. 
    % Each position of the vector has (i,j) that is the direction of a cell in the HT
    %HT = ObjHT.get_HT;                        %HT

    Num_incog = length(Vec_enum);           %Number of cell

    %Sparse matrix
    A = sparse(Num_incog,Num_incog);

    %vectores
    b = zeros(Num_incog,1); %N*M incógnitas, vector columna
    %Matrix and rhs creration 

    %Difference domial extrem
    Dif_dom_x = T.b1-T.a1;
    Dif_dom_y = T.b2-T.a2;
  

    %Loop walks on the all cell in the vector
    for j=1:(Num_incog)
    % for j=1:14
        %Data of the cell
        Id = j;           %Id is the vector's position
        % disp(j)
        %calculate index (i,j) the Id cell (F2 cell)
        cell_id = ObjHT.getCell_id(Id);
        index_i= cell_id.i;
        index_j= cell_id.j;
        l = cell_id.l;
        
        %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
        h_c = (Dif_dom_x)/(T.N*2^l);
        k_c = (Dif_dom_y)/(T.M*2^l);

        %h and k are worked separately, general case not always the same
        hh = h_c^2;
        kk = k_c^2;

        %Weight of the DF scheme of the five points
        pesoh = 1.0/hh;
        pesok = 1.0/kk;

        coordx = T.a1 + index_i*h_c;
        coordy = T.a2 + index_j*k_c;
        
        %Coordinates of the center of the cell
        x_c = coordx + h_c*0.5; 
        y_c = coordy + k_c*0.5;

        %Current cell
        A(Id,Id) = -2*(pesoh+pesok);
        % fprintf('Id: %d \n',Id);
        
        
        %neighbors Id, current cell.
        vecinas_cell = ObjHT.getNeighbors_stencil(Id);
        % disp(vecinas_cell)
        %getVecinas_estencil return Id and level difference of your neighbors
        % fprintf('id: %d \t coordx: %f \t coory: %f \n',Id,coordx,coordy);

        n = size(vecinas_cell,2); %number of neighbors

        m = 0; %counter for direction cell (1,2,3,4)---> (B,L,R,A).

        iter =1; %counter than does not allow repeating the neighboring cell in the same direction to another 
        while iter<=n  %iter is current neighbor 
        % for iter=1:n
            m = m + 1;  %neighbor's direction (1-4).
                
            %If Id = 0, the cell is NULL, so the cell is in counturn
            if vecinas_cell(1,iter) == 0    %Null Cell (0,0) = (Id, level difference).
                peso1 = h_c* pesoh;
                peso2 = k_c* pesok;

                [A,b] = Condicion_Frontera(A,b,Id,Id,pde,T,x_c,y_c,pesoh,pesok,peso1,peso2,m);  %Pesoh = pesok

            %If Id different zero and level difference = 0
            elseif (vecinas_cell(1,iter) ~= 0 && vecinas_cell(2,iter) == 0)

                %cases are generated according to the address because the weights are not necessarily equal
                if (m==1 || m==4)   % B and A cell
                    % fprintf('%d \t %d \t %d \n',m,vecinas_cell(1,iter),vecinas_cell(2,iter));
                    A(Id,vecinas_cell(1,iter)) = A(Id,vecinas_cell(1,iter)) + 1.0/kk;
                else    %L and R cell
                    % fprintf('%d \t %d \t %d \n',m,vecinas_cell(1,iter),vecinas_cell(2,iter));
                    A(Id,vecinas_cell(1,iter)) = A(Id,vecinas_cell(1,iter)) + 1.0/hh;
                end
                
            %Condition Media interpolation, difference level > 0  and Id different zero    
            elseif vecinas_cell(2,iter) > 0
                dl = vecinas_cell(2,iter); % level difference value 
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,pesoh,pesok,vecinas_cell,iter,m,dl);
                
                iter = iter + 2^dl - 1;   %salta a la vecina en otra direccion

            %Condition Fina-gruesa
            else 
                [A,b] = Interpolacion_Fina_Gruesa(ObjHT,T,pde,A,b,Id,index_i,index_j,pesoh,pesok,vecinas_cell,iter,m);
     
            end
            %Counter cell
            iter =iter+1;
        end
        %Rhs 
        b(Id) = b(Id) + pde.fFun(x_c,y_c);
    end
end 


%INTERPOLATIONS

%Boundary condition
function [A,b]= Condicion_Frontera(A,b,Id1,Id2,pde,T,x_c,y_c,pesoh,pesok,peso1,peso2,m)
%This functions calculate the boundary conditios and add the matrix and rhs
    if m == 1   %Below Cell
        A(Id1,Id2) = A(Id1,Id2) - pesok;
        b(Id1) = b(Id1) -2*pesok*pde.g1Fun(x_c,T.a2);

        % A(Id1,Id2) = A(Id1,Id2) + pesok;
        % b(Id1) = b(Id1) + peso2*pde.g1_DerFun(x_c,T.a2);
    elseif m ==2 %Left Cell
        A(Id1,Id2) = A(Id1,Id2) -pesoh;
        b(Id1) = b(Id1) -2*pesoh*pde.g2Fun(T.a1,y_c);

        % A(Id1,Id2) = A(Id1,Id2) + pesoh;
        % b(Id1) = b(Id1) + peso1*pde.g2_DerFun(T.a1,y_c);
    elseif m ==3 %Right Cell
        A(Id1,Id2) = A(Id1,Id2) -pesoh;
        b(Id1) = b(Id1) -2*pesoh*pde.g4Fun(T.b1,y_c);

        % A(Id1,Id2) = A(Id1,Id2) -pesoh;
        % b(Id1) = b(Id1) -2*pesoh*pde.g4Fun(T.b1,y_c);
    else   %Above
        A(Id1,Id2) = A(Id1,Id2) -pesok;
        b(Id1) = b(Id1) -2*pesok*pde.g3Fun(x_c,T.b2);

        % A(Id1,Id2) = A(Id1,Id2) -pesok;
        % b(Id1) = b(Id1) -2*pesok*pde.g3Fun(x_c,T.b2);
    end
end

%Coated cell interpolation  (Coated:recubierta)             
function [A,b] = Interpolacion_Media(ObjHT,A,b,Id,pesoh,pesok,vecinas_cell,iter,m,dl)
    %iter: desde donde inicia la celula fina, m la direccion vecina 1, 2, 3, 4 (B,L;R;A), dl diferencia
    %de nivel  

    if m == 1
        peso = (1.0/(4^dl))*pesok;
        
        %This loop go through the neighboring sisters in the same direction
        %(below)
        for z=iter:iter+(2^dl)-1
            vecina_local = vecinas_cell(1,z);
            
            A(Id,vecina_local) = peso;

            %Recorre hacia abajo todas las vecinas hermanas de la vecina local
            num_vecina_local = 1;
            while num_vecina_local < (2^dl)
                vecina_local = ObjHT.getNeighbors_below_id(vecina_local);
                A(Id,vecina_local) = peso; 
                num_vecina_local = num_vecina_local + 1;
            end
        end

    elseif m == 2
        peso = (1.0/(4^dl))*pesoh;
        %This loop go through the neighboring sisters in the same direction
        %(left)
        for z=iter:iter+2^dl-1
            vecina_local = vecinas_cell(1,z);
            A(Id,vecina_local) = peso;

            %Recorre hacia izquierda todas las vecinas hermanas de la vecina local
            num_vecina_local = 1;
            while num_vecina_local < 2^dl
                vecina_local = ObjHT.getNeighbors_left_id(vecina_local);
                A(Id,vecina_local) = peso; 
                num_vecina_local = num_vecina_local + 1;
            end
        end

    elseif m == 3
        peso = (1.0/(4^dl))*pesoh;
        %This loop go through the neighboring sisters in the same direction
        %(right)
        for z=iter:iter+2^dl-1
            vecina_local = vecinas_cell(1,z);
            A(Id,vecina_local) = peso;
            
            %Recorre hacia derecha todas las vecinas hermanas de la vecina local
            num_vecina_local = 1;
            while num_vecina_local < 2^dl
                vecina_local = ObjHT.getNeighbors_right_id(vecina_local);
                A(Id,vecina_local) = peso; 
                num_vecina_local = num_vecina_local + 1;
            end
        end

    else
        peso = (1.0/(4^dl))*pesok;
     
        %This loop go through the neighboring sisters in the same direction
        %(Above)
        iterparada =iter+2^dl-1;
        
        for z=iter:iterparada
            vecina_local = vecinas_cell(1,z);
            A(Id,vecina_local) = peso;
            % fprintf('%d\t %d\n',Id,vecina_local);

            %Recorre hacia arriba todas las vecinas hermanas de la vecina local
            num_vecina_local = 1;
            while num_vecina_local < 2^dl
                vecina_local = ObjHT.getNeighbors_above_id(vecina_local);
                
                A(Id,vecina_local) = peso; 
                num_vecina_local = num_vecina_local + 1;
            end
        end
    end
end

%Fina-Gruesa Interpolation
function [A,b] = Interpolacion_Fina_Gruesa(ObjHT,T,pde,A,b,Id,index_i,index_j,pesoh,pesok,vecinas_cell,iter,m)
    %iter: vecina actual, m la direccion
    %Pesok y peso h: Peso de la célula central F2 segun la direccion

     vecina_actual = vecinas_cell(1,iter);  %id cell current
     %coordenadas de G1, celula actual
     cell_current = ObjHT.getCell_id(vecina_actual);
     ind_i = cell_current.i;
     ind_j = cell_current.j;
     l = cell_current.l;
      
     %current cell
     h_c = (T.b1 - T.a1)/(T.N*2^l);
     k_c = (T.b2 - T.a2)/(T.M*2^l);

     coordx = T.a1 + ind_i*h_c;   % x value of the curreent cell
     coordy = T.a2 + ind_j*k_c;   % y value of the curreent cell
        
     %Coordinates of the center of the cell
     x_c = coordx + h_c*0.5; 
     y_c = coordy + k_c*0.5;
    
    if m == 1   %neighbor below, G1, m =1 below
        % peso = 1/kk;
        
        A(Id,Id) = A(Id,Id) + (2.0/3.0)*pesok;   %celula donde está parado, la central F2
        % vecina_actual = vecinas_cell(1,iter);  %id cell current

        A(Id,vecina_actual) = A(Id,vecina_actual) + 0.5*pesok;  %G1

        %neighbor left of the G1 
        G2_list = ObjHT.getNeighbors_left_id_dif_level(vecina_actual);   %G2
        %G2_list has the id and the difference in level of the neighbors

        if G2_list(1,1) == 0.0
            %haga frontera;
            %Pesos
            PesoC1 = pesok*(1.0/12.0);  %Peso below-left
            PesoC2 = pesok*(-1.0/20.0);  %Peso below_Right

            PesoC1_rh = k_c*pesok*(1.0/12.0);
            PesoC2_rh = k_c*pesok*(-1.0/20.0);

            m=2;   %diretion the neighbor (left)

            if mod(index_i,2) == 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m);  %case 1, pesos normales (C1,c2)
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m);  %case 2, pesos invertidos (C2,C1)
            end 

        elseif G2_list(2,1) > 0
            %Media
            m=2;
            iter = 1;
            dl = G2_list(2,1); 
            PesoC1 = pesok*(1.0/12.0);  %Peso below-left
            PesoC2 = pesok*(-1.0/20.0);  %Peso below_Right

            if mod(index_i,2) == 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G2_list,iter,m,dl);
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G2_list,iter,m,dl);
            end
        else
            G2 = G2_list(1,1);

            if mod(index_i,2) == 0
                A(Id,G2) = A(Id,G2) + (1.0/12.0)*pesok;
            else
                A(Id,G2) = A(Id,G2) - (1.0/20.0)*pesok;
            end
        end
        
        %neighbor right of the G1 
        G3_list = ObjHT.getNeighbors_right_id_dif_level(vecina_actual); %G3
        %G3_list has the id and the difference in level of the neighbors

        if G3_list(1,1) == 0.0
            %haga frontera;
            %Pesos
            PesoC1 = pesok*(-1.0/20.0);  %Peso below-Right
            PesoC2 = pesok*(1.0/12.0);  %Peso below_left

            PesoC1_rh = k_c*pesok*(-1.0/20.0);
            PesoC2_rh = k_c*pesok*(1.0/12.0);

            m=3;   %diretion the neighbor (right)

            if mod(index_i,2) == 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m);   %orden de los pesos está acorde, cambian segun el caso
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m);  %case 2, pesos invertidos (C2,C1)
            end 
            
        elseif G3_list(2,1) > 0
            %Media
            m=3;
            iter = 1;
            dl = G3_list(2,1); 
            PesoC1 = pesok*(-1.0/20.0);  %Peso below_Right
            PesoC2 = pesok*(1.0/12.0);  %Peso below-left

            if mod(index_i,2) == 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G3_list,iter,m,dl);
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G3_list,iter,m,dl);
            end

        else
            G3 = G3_list(1,1);
            if mod(index_i,2) == 0
                A(Id,G3) = A(Id,G3) - (1.0/20.0)*pesok;
            else
                A(Id,G3) = A(Id,G3) + (1.0/12.0)*pesok;
            end
        end

        %neighbor above of F2
        F1 = ObjHT.getNeighbors_above_id(Id);
        A(Id,F1) = A(Id,F1) + (-1.0/5.0)*pesok;

        % fprintf('F2:%.lf \t G1:%.lf \t G2:%.lf \t G3: %.lf\t F1: %.lf',A(Id,Id),A(Id,vecina_actual),A(Id,G2),A(Id,G3),A(Id,F1));
        

    elseif m == 2   %neighbor left , G1, m =2 left
        % peso = 1/hh;
        
        A(Id,Id) = A(Id,Id) + (2.0/3.0)*pesoh;   %celula donde está parado, la central 
        % vecina_actual = vecinas_cell(1,iter);  %id cell current
        A(Id,vecina_actual) = A(Id,vecina_actual) + 0.5*pesoh;
        
         %neighbor below of the G1 
         G3_list = ObjHT.getNeighbors_below_id_dif_level(vecina_actual);   %G3
         %G3_list has the id and the difference in level of the neighbors

        if G3_list(1,1) == 0.0
            %haga frontera;
            %Pesos    
            PesoC1 = pesoh*(1.0/12.0);  %Peso Right-below
            PesoC2 = pesoh*(-1.0/20.0);  %Peso left-below

            PesoC1_rh = h_c*pesoh*(1.0/12.0);
            PesoC2_rh = h_c*pesoh*(-1.0/20.0);

            m=1;   %diretion the neighbor ()

            if mod(index_j,2) ~= 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m); 
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m); 
            end    

        elseif G3_list(2,1) > 0
            %Media
            m=1;
            iter = 1;
            dl = G3_list(2,1); 
            PesoC1 = pesoh*(1.0/12.0);  %Peso Right-below
            PesoC2 = pesoh*(-1.0/20.0);  %Peso Left-below

            if mod(index_j,2) ~= 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G3_list,iter,m,dl); 
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G3_list,iter,m,dl); 
            end 

        else
            G3 = G3_list(1,1);
            if mod(index_j,2) ~= 0
                A(Id,G3) = A(Id,G3) - (1.0/20.0)*pesoh;
            else
                A(Id,G3) = A(Id,G3) + (1.0/12.0)*pesoh; 
            end 
  
        end

        %neighbor above of the G1 
        G2_list = ObjHT.getNeighbors_above_id_dif_level(vecina_actual);
        %G2_list has the id and the difference in level of the neighbors

        if G2_list(1,1) == 0.0
            %haga frontera;
            %Pesos
            PesoC1 = pesoh*(-1.0/20);  %Peso left-below
            PesoC2 = pesoh*(1.0/12.0);  %Peso Right-below

            PesoC1_rh = h_c*pesoh*(-1.0/20.0);
            PesoC2_rh = h_c*pesoh*(1.0/12.0);

            m=4;   %diretion the neighbor (above)

            if mod(index_j,2) ~= 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m); 
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m);  
            end 

        elseif G2_list(2,1) > 0
            %Media
            m=4;
            iter = 1;
            dl = G2_list(2,1); 

            PesoC1 = pesoh*(-1.0/20.0);  %Peso above_Right
            PesoC2 = pesoh*(1.0/12.0);  %Peso above-left
            
            if mod(index_j,2) ~= 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G2_list,iter,m,dl); 
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G2_list,iter,m,dl); 
            end 
            
        else
            G2 = G2_list(1,1);
            if mod(index_j,2) ~= 0
                A(Id,G2) = A(Id,G2) + (1.0/12.0)*pesoh;
            else
                A(Id,G2) = A(Id,G2) - (1.0/20.0)*pesoh;
            end
        end

        F1 = ObjHT.getNeighbors_right_id(Id);
        A(Id,F1) = A(Id,F1) + (-1.0/5.0)*pesoh;

    elseif m == 3  %neigbor right
        % peso = 1/hh;
        A(Id,Id) = A(Id,Id) + (2.0/3.0)*pesoh;   %celula donde está parado, la central 
        % disp(A(Id,Id))
        vecina_actual = vecinas_cell(1,iter);  %id cell current
        A(Id,vecina_actual) = A(Id,vecina_actual) + 0.5*pesoh;

        %neighbor below of the G1 
        G2_list = ObjHT.getNeighbors_below_id_dif_level(vecina_actual);   %G2
        %G2_list has the id and the difference in level of the neighbors    
        if G2_list(1,1) == 0.0
            %haga frontera;
            %Pesos    
            PesoC1 = pesoh*(-1.0/20.0);  %Peso Right-below
            PesoC2 = pesoh*(1.0/12.0);  %Peso left-below

            PesoC1_rh = h_c*pesoh*(-1.0/20.0);
            PesoC2_rh = h_c*pesoh*(1.0/12.0);
            
            m=1;   %diretion the neighbor () input peso2

            if mod(index_j,2) == 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m); 
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m); 
            end

        elseif G2_list(2,1) > 0
            %Media
            m=1;
            iter = 1;
            dl = G2_list(2,1); 
            PesoC1 = pesoh*(-1.0/20.0);  %Peso Right-below
            PesoC2 = pesoh*(1.0/12.0);  %Peso Left-below

            if mod(index_j,2) == 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G2_list,iter,m,dl);
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G2_list,iter,m,dl); 
            end
        else
            G2 = G2_list(1,1);
            if mod(index_j,2) == 0
                A(Id,G2) = A(Id,G2) + (1.0/12.0)*pesoh;
            else
                A(Id,G2) = A(Id,G2) - (1.0/20.0)*pesoh;
            end
        end

        %neighbor above of the G1
        G3_list = ObjHT.getNeighbors_above_id_dif_level(vecina_actual);  %G3
        %G3_list has the id and the difference in level of the neighbors

        if G3_list(1,1) == 0.0
            %haga frontera;
            %Pesos    
            PesoC1 = pesoh*(1.0/12.0);  %Peso Right-above
            PesoC2 = pesoh*(-1.0/20.0);  %Peso left-above

            PesoC1_rh = h_c*pesoh*(1.0/12.0);
            PesoC2_rh = h_c*pesoh*(-1.0/20.0);
            
            m=4;   %diretion the neighbor ()

            if mod(index_j,2) == 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m); 
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m); 
            end

        elseif G3_list(2,1) > 0
            %Media
            m=4;
            iter = 1;
            dl = G3_list(2,1); 
            PesoC1 = pesoh*(1.0/12.0);  %Peso Right-above
            PesoC2 = pesoh*(-1.0/20.0);  %Peso left-above

            if mod(index_j,2) == 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G3_list,iter,m,dl);
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G3_list,iter,m,dl); 
            end
        else
            G3 = G3_list(1,1);
            if mod(index_j,2) == 0
                A(Id,G3) = A(Id,G3) - (1.0/20.0)*pesoh;
            else
                A(Id,G3) = A(Id,G3) + (1.0/12.0)*pesoh;
            end
        end

        F1 = ObjHT.getNeighbors_left_id(Id);
        A(Id,F1) = A(Id,F1) + (-1.0/5.0)*pesoh;

    else  %(m=4)
        % peso = 1/kk; %neigbor above

        A(Id,Id) = A(Id,Id) + (2.0/3.0)*pesok;   %celula donde está parado, la central 
        vecina_actual = vecinas_cell(1,iter);  %id cell current
        A(Id,vecina_actual) = A(Id,vecina_actual) + 0.5*pesok;

        %neighbor left of the G1 
        G3_list = ObjHT.getNeighbors_left_id_dif_level(vecina_actual);   %G3
            
        if G3_list(1,1) == 0.0
            %haga frontera;
            %Pesos    
            PesoC1 = pesok*(-1.0/20.0);  %Peso left-above  
            PesoC2 = pesok*(1.0/12.0);  %Peso Right-above  

            PesoC1_rh = k_c*pesok*(-1.0/20.0);
            PesoC2_rh = k_c*pesok*(1.0/12.0);

            m=2;   %diretion the neighbor ()

            if mod(index_i,2) ~= 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m);
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m);
            end

        elseif G3_list(2,1) > 0
            %Media
            m=2;
            iter = 1;
            dl = G3_list(2,1); 
            PesoC1 = pesok*(-1.0/20.0);  %Peso left-above 
            PesoC2 = pesok*(1.0/12.0);  %Peso Right-above

            if mod(index_i,2) ~= 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G3_list,iter,m,dl);
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G3_list,iter,m,dl);
            end
            
        else
            G3 = G3_list(1,1);
            if mod(index_i,2) ~= 0
                A(Id,G3) = A(Id,G3) - (1.0/20.0)*pesok;
            else
                A(Id,G3) = A(Id,G3) + (1.0/12.0)*pesok;
            end
            
        end

        %neighbor right of the G1 
        G2_list = ObjHT.getNeighbors_right_id_dif_level(vecina_actual);   %G2
        %G2_list has the id and the difference in level of the neighbors
        
        if G2_list(1,1) == 0.0
            %haga frontera;
            %Pesos    
            PesoC1 = pesok*(1.0/12.0);  %Peso Right-above
            PesoC2 = pesok*(-1.0/20.0);  %Peso left-above

            PesoC1_rh = k_c*pesok*(1.0/12.0);
            PesoC2_rh = k_c*pesok*(-1.0/20.0);

            m=3;   %diretion the neighbor ()
                
            if mod(index_i,2) ~= 0
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC1,PesoC2,PesoC1_rh,PesoC2_rh,m);
            else
                [A,b]=Condicion_Frontera(A,b,Id,vecina_actual,pde,T,x_c,y_c,PesoC2,PesoC1,PesoC2_rh,PesoC1_rh,m);
            end   

        elseif G2_list(2,1) > 0
            %Media
            m=3;
            iter = 1;
            dl = G2_list(2,1); 
            PesoC1 = pesok*(1.0/12.0);  %Peso Right-above
            PesoC2 = pesok*(-1.0/20.0);  %Peso left-above
            
            if mod(index_i,2) ~= 0
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC1,PesoC2,G2_list,iter,m,dl);
            else
                [A,b] = Interpolacion_Media(ObjHT,A,b,Id,PesoC2,PesoC1,G2_list,iter,m,dl);
            end
            
        else
            G2 = G2_list(1,1);
            if mod(index_i,2) ~= 0
                A(Id,G2) = A(Id,G2) + (1.0/12.0)*pesok;
            else
                A(Id,G2) = A(Id,G2) - (1.0/20.0)*pesok;
            end
        end

        F1 = ObjHT.getNeighbors_below_id(Id);
        A(Id,F1) = A(Id,F1) + (-1.0/5.0)*pesok;
    end
end

