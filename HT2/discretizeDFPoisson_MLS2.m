function [A,b] = discretizeDFPoisson_MLS2(T,pde,ObjHT)
    %This function creates the matrix A and vector b (rhs), from the
    %discretization of FD
    % %input grid, pde and HT
    %Use the vector to calculate the id, then calculate neighbors
    %Tiene MLS2, que hace referencia a ML3 con la forma de MLS1

    Vec_enum = ObjHT.get_enumeration_vector();   %This vector has the p and key de each cell. 
    % Each position of the vector has (p,key) that are the data of a cell in the HT
    HT = ObjHT.get_HT; 

    Num_incog = length(Vec_enum);           %Number of cell

    %Sparse matrix
    A = sparse(Num_incog,Num_incog);

    %vectors
    b = zeros(Num_incog,1);                 %N*M unknows, column vector
   
    %Matrix and rhs creration 

    %Difference domial extrem
    Dif_dom_x = T.b1-T.a1;
    Dif_dom_y = T.b2-T.a2;
   
    %Matrix and rhs creration 

    num_max_point1 = 0;
    num_min_point = 100;
    vec_cont_num_punt = zeros(25,1);
    %Loop walks on the all cell in the vector
    fileID3 = fopen('File_cell_with_i_points_MLS2.txt','w');

    for j=1:(Num_incog)
    % for j=205:205
        %Data of the cell
        %Data of the cell
        Id = j;           %Id is the vector's position
        % disp('center_cell')
        % disp(Id)
        %calculate index (i,j) the Id cell (F2 cell)
        p = Vec_enum{Id}(1);
        key = Vec_enum{Id}(2);
        % disp(p)
        % disp(key)

        %Return the struct of cell with this key
        cell_id = HT{p}.get_cell_Nodo(key);
        % cell_id = ObjHT.getCell_id(Id);

        index_i= cell_id.i;
        index_j= cell_id.j;
        l = cell_id.l;

        % fprintf('id: %d \t i: %d \t j: %d \t l: %d \n',Id,index_i,index_j,l);
        
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

        %calculate index (i,j) the Id cell (F2 cell)
        % index_i= T.N*2^(l)*coordx;
        % index_j= T.M*2^(l)*coordy;
   
        %Current cell
        A(Id,Id) = -2*(pesoh+pesok);
        
        %neighbors Id(Cell_C), current cell.
        vecinas_cell = ObjHT.getNeighbors_stencil(Id); 
        %getVecinas_estencil return Id and level difference of your neighbors
        % fprintf('id: %d \t coordx: %f \t coory: %f \n',Id,coordx,coordy);

        n = size(vecinas_cell,2); %number of neighbors 

        m = 0; %counter for direction cell (1,2,3,4)---> (B,L,R,A).

        iter =1; %counter than does not allow repeating the neighbor cell in the same direction to another 
        while iter<=n  %iter is current neighbor 
        % for iter=1:n
            m = m + 1;  %neighbor's direction (1-4).
                
            %If Id = 0, the cell is NULL, so the cell is in counturn
            if vecinas_cell(1,iter) == 0    %Null Cell (0,0) = (Id, level difference).

                peso1 = h_c* pesoh;
                peso2 = k_c* pesok;

                [A,b] = Condicion_Frontera(A,b,Id,Id,pde,T,x_c,y_c,pesoh,pesok,peso1,peso2,m);  %Pesoh = pesok

            %If Id different zero and level difference = 0, the cell is in
            %the same level
            elseif (vecinas_cell(1,iter) ~= 0 && vecinas_cell(2,iter) == 0)
                %cases are generated according to the address because the weights are not necessarily equal
                if (m==1 || m==4)   % B and A cell
                    A(Id,vecinas_cell(1,iter)) = A(Id,vecinas_cell(1,iter)) + pesok;
                else    %L and R cell
                    A(Id,vecinas_cell(1,iter)) = A(Id,vecinas_cell(1,iter)) + pesoh;
                end
                
            %Condition MLS aproximation, difference level > 0  and Id different zero    
            elseif vecinas_cell(2,iter) ~= 0
                ID_CELL_CENT = Id;
                Id_cel_vecina = vecinas_cell(1,iter);  %central cell's neighbor in a level different, firts neighbor
                dl = vecinas_cell(2,iter); % level difference value 

                % if dl >= 1 
                %     Id_cel_vecina = vecinas_cell(1,iter+1);  %central cell's neighbor in a level different, second neighbor
                %     dl = vecinas_cell(2,iter+1); % level difference value 
                % end


                %tercer y cuarto caso de vecinas
                % if m == 1
                %     DIR = 'B';
                % elseif m ==2
                %     DIR = 'L';
                % elseif m ==3
                %     DIR = 'R';
                % else
                %     DIR = 'A';
                % end
                % vecinos_dir_1 = ObjHT.getNeighbors_direction_id(Id_cel_vecina,DIR);
                % Id_cel_vecina = vecinos_dir_1(1,1);
                % cell_id_ = ObjHT.getCell_id(Id_cel_vecina);
                % dl = cell_id_.l;
     
                [A,b, num_max_points] = MLS_function_3(T,ObjHT,A,b,ID_CELL_CENT,Id_cel_vecina,pesoh,pesok,m,dl);

                if dl >= 1   %The neighbor is in the finer level
                    if m == 1
                        DIR = 'B';
                    elseif m ==2
                        DIR = 'L';
                    elseif m ==3
                        DIR = 'R';
                    else
                        DIR = 'A';
                    end
                    vecinos_dir = ObjHT.getNeighbors_direction_id(Id,DIR);
                    num_vecinas_igual_direccion = length(vecinos_dir); %salta a la vecina en otra direccion
                    % disp('aqui')
                    % disp(num_vecinas_igual_direccion)
                    % num_vecinas_igual_direccion = 2^(dl)-1;  %salta a la vecina en otra direccion
                    iter = iter + num_vecinas_igual_direccion -1;
                end 

                %Cuenta cuantos puntos tienen i puntos
                vec_cont_num_punt(num_max_points) = vec_cont_num_punt(num_max_points) + 1;
                
                %Cuenta los puntos máximos y minimos en la lista de pesos
                %en MLS2
     
                if num_max_points >= num_max_point1
                    num_max_point1 = num_max_points;
                end

                if num_max_points <= num_min_point
                    num_min_point = num_max_points;
                end


            end
            %Counter cell
            iter =iter+1;
        end
        %Rhs 
        b(Id) = b(Id) + pde.fFun(x_c,y_c);
    end
    % % fprintf(fileID3,'num_max_points: %d \t num_min_points: %d\n', num_max_point1,num_min_point);
    % % fprintf('num_max_points: %d \n num_min_points: %d\n', num_max_point1,num_min_point);
    % % fprintf(fileID3,'\n');
    % % fprintf(fileID3,'%d \t',vec_cont_num_punt);
    % % disp(vec_cont_num_punt);
    % % fclose(fileID3);

    fprintf(fileID3,'%d\n', num_max_point1); %Point max number
    fprintf(fileID3,'%d\n', num_min_point); %Point min number
    fprintf(fileID3,'%d \t',vec_cont_num_punt); %Points vector
    % disp(vec_cont_num_punt);
    fclose(fileID3);
end 

%INTERPOLATION

%Boundary condition
function [A,b]= Condicion_Frontera(A,b,Id1,Id2,pde,T,x_c,y_c,pesoh,pesok,peso1,peso2,m)
%This functions calculate the boundary conditios and add the matrix and rhs

%Firts Dirishlet conditions and seconds Neumann contitions
    if m == 1   %Below Cell
        A(Id1,Id2) = A(Id1,Id2) - pesok;
        b(Id1) = b(Id1) -2*pesok*pde.g1Fun(x_c,T.a2);

        % A(Id1,Id2) = A(Id1,Id2) + pesok;
        % b(Id1) = b(Id1) + peso2*pde.g1_DerFun(x_c,T.a2);
    elseif m ==2 %Left Cell
        A(Id1,Id2) = A(Id1,Id2) - pesoh;
        b(Id1) = b(Id1) -2*pesoh*pde.g2Fun(T.a1,y_c);

        % A(Id1,Id2) = A(Id1,Id2) + pesoh;
        % b(Id1) = b(Id1) + peso1*pde.g2_DerFun(T.a1,y_c);
    elseif m ==3 %Right Cell
        A(Id1,Id2) = A(Id1,Id2) - pesoh;
        b(Id1) = b(Id1) -2*pesoh*pde.g4Fun(T.b1,y_c);

        % A(Id1,Id2) = A(Id1,Id2) + pesoh;
        % b(Id1) = b(Id1) + peso1*pde.g3_DerFun(T.b1,y_c);
    else   %Above
        A(Id1,Id2) = A(Id1,Id2) -pesok;
        b(Id1) = b(Id1) -2*pesok*pde.g3Fun(x_c,T.b2);

        % A(Id1,Id2) = A(Id1,Id2) + pesok;
        % b(Id1) = b(Id1) + peso2*pde.g4_DerFun(x_c,T.b2);
    end
end

