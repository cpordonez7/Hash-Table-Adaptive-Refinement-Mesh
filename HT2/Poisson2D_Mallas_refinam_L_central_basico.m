%grid structure , Domine, (a1,b1)x(a2,b2)
clear
 grid.a1 = 0.0;
 grid.a2 = 0.0;
 grid.b1 = 1.0;
 grid.b2 = 1.0;

 %Partitions number, the firts N and M
 N = 4;
 M = 4;
 Rep_Final = 1;  %Doubling of N and M.

 %Temporal vector
 Error_Nor2 = zeros(Rep_Final,1);
 Error_Norh = zeros(Rep_Final,1);
 Error_NormInf = zeros(Rep_Final,1);

  %Difference domial extrem
  Dif_dom_x = grid.b1 - grid.a1;
  Dif_dom_y = grid.b2 - grid.a2;

 
 %Estructur PDE, input number exemplo
 %1 ejm Tesis Rua
 %2 ejm tarea
 %3 ejm arti Andrea  [0,2]x[0,1]
 %4 ejm test gota tanh, Morar la aplitud
 %5 ejm propio manufacturado
 pde = getPDE(1);

 %Create file to write result 
fileID1 = fopen('Result_numer_Poisson_2D_Refinamineto_basico_central.txt','w');
fprintf(fileID1,'\t\t RESULTADOS NUMERICOS POISSON 2D \n\n');

fileID2 = fopen('Result_numer_Poisson_2D_bicg.txt','w');
fprintf(fileID2,'\t\t RESULTADOS NUMERICOS POISSON 2D \n\n');

% fileensayo = fopen('Result_ensayo_Poisson_2D.txt','w');
% fprintf(fileensayo,' x \t y \t\t Aprox \t\t  Exacta \t\t |Aprox-exac| \t\t h_i\n');


fprintf(fileID1,' N \t Norma_2 \t\t Norma_h \t\t  Norma_Inf \t\t Num Cond \t\t Vlr propio max \t Vlr propio min  \t\t Simétrica \t\t max colisiones\n');
fprintf(fileID2,' N \t Norma_2 \t\t Norma_h \t\t  Norma_Inf \t\t Num Cond \t\t Vlr propio max \t Vlr propio min  \t\t Simétrica \t\t max colisiones\n');

%loop for doubling of N and M
 for RepIni=1:Rep_Final
     grid.N = 2^(RepIni-1)*N;
     grid.M = 2^(RepIni-1)*M;
     fprintf(fileID1,'%d x %d &\t',grid.N, grid.M);
     fprintf(fileID2,'%d &\t',grid.N);

     
     profile on    
     tic;
     if RepIni==1
         %Malla con refinamineto de L en el centro básico
         ObjHT = New_Mesh2D('Malla_refinada_L_centro_N_4.txt',2,grid); 
        
     else
         % ObjHT = New_Mesh2D('File_current_mesh.txt',2,grid);
         ObjHT = New_Mesh2D('File_current_mesh.txt',2,grid);
     end
     
     time_create_object = toc;
     disp(['time_create_object: ' num2str(time_create_object) ' segundos']);

    %ENUMERACIÓN POR NIVEL
     % ObjHT.set_id_celulas_for_level
     % ObjHT.set_enumeration_vector

     %PLOT grid
     % ObjHT.graficar();
     % ObjHT.graficar_enumeracion();

     %Calculate values for step size
     h = (Dif_dom_x)/grid.N;
     k = (Dif_dom_y)/grid.M;

     % %mesh ends in center variable
     x1 = grid.a1 + h*0.5;   %First element of the mesh
     y1 = grid.a2 + k*0.5;
     xn = grid.b1 - h*0.5;   %Last element of the mesh
     yn = grid.b2 - k*0.5;
     h2 = h*h;

     x = linspace(x1,xn,grid.N); %Vector de N puntos equiespaciado entre x1 y xn
     y = linspace(y1,yn,grid.M); %Vector de M puntos equiespaciado entre x1 y xn

     %Versiones finales
     tic;
     % [A,B] = discretizeDFPoisson_Inter_orden_uno(grid,pde,ObjHT);
     % [A,B] = discretizeDFPoisson_Inter_cuadratica(grid,pde,ObjHT); 
     [A,B] = discretizeDFPoisson_MLS1(grid,pde,ObjHT);
     % [A,B] = discretizeDFPoisson_MLS2(grid,pde,ObjHT);
     %Input grid, pde and HT.
     time_create_matrix_vector = toc;
     disp(['time_create_matrix_vector: ' num2str(time_create_matrix_vector) ' segundos']);
     
     % Aprox = pcg(A,B,10^-10);   %Gradiente conjugado, para A simétrica
     tic;
      Aprox1 = A\B; %system solution for gauss
     % Aprox1 = bicg(A,B,1e-10,900);   %Gradiente Bi-Conjugado
     time_solver_system = toc;
     disp(['time_solver_system Met Dir: ' num2str(time_solver_system) ' segundos']);

     % % tic;
     % % % Aprox2 = bicg(A,B,1e-10,900);   %gradiente biconjugado
     % % Aprox2 = bicg(A,B,1e-8,1300);  %Caso N=128
     % % time_solver_system = toc;
     % % disp(['time_solver_system Met Bicg: ' num2str(time_solver_system) ' segundos']);

    profile off

    %Plot of sparse matriz
    % spy(A,'r*',7);  


    %calculate enumeration vector and HT
    Vec_enumeation = ObjHT.get_enumeration_vector();   %Vector whit the cells
    % HT = ObjHT.get_HT;

    %Cell number
    num_incog = length(Vec_enumeation);

    %Vector solucion exacta
    Sol_ex = zeros(num_incog,1);

    %Vector de hs, en cada cell
    Vector_h = zeros(num_incog,1);
    x_c = zeros(num_incog,1);
    y_c = zeros(num_incog,1);

    %SOL EXACTA, 
    for j=1:(num_incog)
        cell_j = ObjHT.getCell_id(j); 
        ind_i = cell_j.i;
        ind_j = cell_j.j;
        l = cell_j.l;
        Id = cell_j.id;

        %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
        h_c = (Dif_dom_x)/(grid.N*2^l);
        k_c = (Dif_dom_y)/(grid.M*2^l);
        Vector_h(j) = h_c^(0.5)*k_c^(0.5);
        
        coordx = grid.a1 + ind_i*h_c;    % x value of the curreent cell
        coordy = grid.a2 + ind_j*k_c;    % y value of the curreent cell
    
        %Coordinates of the center of the cell
        x_cen = coordx + h_c*0.5;
        y_cen = coordy + k_c*0.5;
        x_c(j) = x_cen;
        y_c(j) = y_cen;
        Sol_ex(Id) = pde.solexFun(x_c(j),y_c(j));

        % fprintf(fileensayo,'%d \t %.5f \t %.5f \t %.10e \t %.10e \t %.5e \t %.5e \n',Id,x_c(j),y_c(j),Aprox2(j), Sol_ex(j),abs(Aprox2(j)-Sol_ex(j)),Vector_h(j));
        % fprintf(fileensayo,'\n');

    end

    %Calculo de errores METODO DIRECTO
    % Dif = Sol_ex - Aprox;
    Dif_norm_h_1 = (Sol_ex - Aprox1).*Vector_h;   
    Error_Nor2(RepIni) = norm(Sol_ex-Aprox1,2);
    Error_Norh(RepIni) = sqrt(sum(Dif_norm_h_1.^2));
    % Error_Norh(RepIni) = norm((Sol_ex - Aprox).*Vector_h,2);
    % Error_Norh(RepIni) = h*Error_Nor2(RepIni);
    Error_NormInf(RepIni) = norm(Sol_ex-Aprox1,Inf);
    fprintf(fileID1,'%.4e &\t%.4e &\t%.4e &\t',Error_Nor2(RepIni), Error_Norh(RepIni), Error_NormInf(RepIni));
 

    %calculo valores propios y numero de condicionamiento
    % % % valores_propios = eigs(A);
    % % % vlr_prop_max = max(valores_propios);
    % % % vlr_prop_min = min(valores_propios);
    % logial_value =  issymmetric(A);
    % % % fprintf(fileID1,'%.4e &\t%.4e &\t%.4e\n',condest(A), vlr_prop_max, vlr_prop_min);
    fprintf(fileID1,'%.4e \n',condest(A));
    % fprintf(fileID1,'%.10f\t%.10f\t%.10f\n',Error_Nor2(RepIni), Error_Norh(RepIni), Error_NormInf(RepIni));
    
    % % % %Calculo de errores METODO BICG
    % % % % Dif = Sol_ex - Aprox;
    % % % Dif_norm_h_2 = (Sol_ex - Aprox2).*Vector_h;   
    % % % Error_Nor2(RepIni) = norm(Sol_ex-Aprox2,2);
    % % % Error_Norh(RepIni) = sqrt(sum(Dif_norm_h_2.^2));
    % % % % Error_Norh(RepIni) = norm((Sol_ex - Aprox).*Vector_h,2);
    % % % % Error_Norh(RepIni) = h*Error_Nor2(RepIni);
    % % % Error_NormInf(RepIni) = norm(Sol_ex-Aprox2,Inf);
    % % % fprintf(fileID2,'%.4e &\t%.4e &\t%.4e &\t',Error_Nor2(RepIni), Error_Norh(RepIni), Error_NormInf(RepIni));
    % % % 
    % % % % fprintf(fileID,'\n');
    % % % 
    % % % fprintf(fileID2,'%.4e &\t%.4e &\t%.4e\n',condest(A), vlr_prop_max, vlr_prop_min);
    % % % 
    % % % 
    % % % fprintf(fileID2,'%.10f\t%.10f\t%.10f\n',Error_Nor2(RepIni), Error_Norh(RepIni), Error_NormInf(RepIni));

    % % % vec_collisions = ObjHT.get_number_collisions();
    % % % max_collisions = max(vec_collisions);
    % % % fprintf(fileID1,'%.d\n',max_collisions);
    % % % fprintf(fileID2,'%.d\n',max_collisions);

    %calculo error iterativo
    % error_ensayo = 0.0;
    % for j=1:(num_incog)
    %     error_ensayo = error_ensayo + (abs((Sol_ex(j)-Aprox(j))*Vector_h(j)))^2;
    % end
    % NORMAH_ENSAYO = sqrt(error_ensayo);

    % fprintf(fileID,'Normah_ensayo: %.10e \n',NORMAH_ENSAYO);


     %PRINT EXACT ANS NUMERICAL SOLUTIONS
    % for iter=1:N*M
    %     fprintf('%.7e \t %.7e %.5e \n',Sol_ex(iter),Aprox(iter),Dif(iter));
    % end

    % ObjHT.Refinar_malla();
    % ObjHT.Genera_archivo_malla_actual();

    %Updata merf file
    if RepIni ~= Rep_Final
        ObjHT.Refine_Mesh;
        ObjHT.Generate_current_mesh_file;
        Edite_level_of_file('File_current_mesh.txt');
    end
    
 end

 fprintf(fileID1,'\n');
 fprintf(fileID2,'\n');


 %Calculo de las razones 
 fprintf(fileID1,'Razon Norma_2 \t Razon Norma_h \t Razon Norma_Inf \n');
 % fprintf(fileID2,'Razon Norma_2 \t Razon Norma_h \t Razon Norma_Inf \n');
 for iter=1:Rep_Final-1
     Razon1 = Error_Nor2(iter)/Error_Nor2(iter+1);
     Razon2 = Error_Norh(iter)/Error_Norh(iter+1);
     Razon3 = Error_NormInf(iter)/Error_NormInf(iter+1);
     fprintf(fileID1,'%.4f &\t  %.4f &\t %.4f \n',Razon1, Razon2, Razon3);
     % fprintf(fileID2,'%.4e \t  %.4e \t %.4e \n',Razon1, Razon2, Razon3);
 end

% Grafica
% Z1 = reshape(Sol_ex,grid.N,grid.M); %Solucion exacta
% Z2=reshape(Aprox,grid.N,grid.M);   %pasa de vector a matriz, Aproximacion


%Plot solución   
% [X,Y] = meshgrid(x,y);
% contour(X,Y,Z2);
% surf(X,Y,Z2);
% surf(X,Y,Z2);
% plot3(X,Y,Z2)
% meshc(X,Y,Z2)
% pcolor(X,Y,Z2)


%Usar para grafica de aprox o error
 %plot3(x_c,y_c,Aprox1,'*r','MarkerSize',10)

% fclose(fileensayo);
fclose(fileID1);
fclose(fileID2);


