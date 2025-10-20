%grid structure , Domine, (a1,b1)x(a2,b2)
clear
 grid.a1 = 0.0;
 grid.a2 = 0.0;
 grid.b1 = 1.0;
 grid.b2 = 1.0;

 % % %Domain grid test 2
 % grid.a1 = -0.5;
 % grid.a2 = -0.5;
 % grid.b1 = 0.5;
 % grid.b2 = 0.5;

 %Partitions number, the firts N and M
 N = 16;
 M = 16;
 Rep_Final = 2;  %Doubling of N and M.

 %Temporal vector
 Error_Nor2 = zeros(Rep_Final,1);
 Error_Norh = zeros(Rep_Final,1);
 Error_NormInf = zeros(Rep_Final,1);
 tasa_ocupacion_NO_nula = zeros(Rep_Final,1); 
 tasa_ocupacion_nula = zeros(Rep_Final,1); 

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
fileID1 = fopen('Resultados_datos_matriz_discretizacion.txt','w');
fileID2 = fopen('Celulas_por_nivel_en_la_malla.txt','w');
% fileensayo = fopen('Result_ensayo_Poisson_2D.txt','w');
% fprintf(fileensayo,' x \t y \t\t Aprox \t\t  Exacta \t\t |Aprox-exac| \t\t h_i\n');


fprintf(fileID1,' N \t posi_Matriz \t\t posit_Nulas \t\t  Posit_NO_nulas \n');
fprintf(fileID2,' N \t cell en el nivel  i\n');

%loop for doubling of N and M
 for RepIni=1:Rep_Final
     grid.N = 2^(RepIni-1)*N;
     grid.M = 2^(RepIni-1)*M;
     fprintf(fileID1,'%d x %d & \t',grid.N, grid.M);
     fprintf(fileID2,'%d x %d & \t',grid.N, grid.M);

     % Function to create a file of the grid uniform
     % Create_uniform_mesh(grid.N,grid.M,grid.a1,grid.b1,grid.a2,grid.b2);  %N:particiones y extremos del dominio en un eje , N=M
    
     %HT param (file, max nevel,particiones eje, intervalo)
     % ObjHT = New_Mesh2D('malla_c_4.txt',1,4,0,1);
     % ObjHT = New_Mesh2D('Malla_uniforme_base_32.txt',0,grid.N,grid.a1, grid.b1);
     
     % ObjHT = New_Mesh2D('File_malla_4_niv_refinada_N_8.txt',4,grid.N,grid.a1, grid.b1); %Colocar el doble de niveles porque se va a refinar toda la malla
     % ObjHT = New_Mesh2D('Malla_uniforme_base.txt',2,grid); %Colocar el doble de niveles porque se va a refinar toda la malla
     
     
     profile on    
     tic;
     if RepIni==1
         % ObjHT = New_Mesh2D('File_current_mesh.txt',3,grid);
         % ObjHT = New_Mesh2D('Malla_uniforme_base_32.txt',1,grid);
     % ObjHT = New_Mesh2D('Base_uniform_mesh.txt',1,grid);
     % ObjHT = New_Mesh2D('Base_uniform_mesh_6x4_aleatoria.txt',1,grid);
     % ObjHT = New_Mesh2D('Mesh_base_uniforme_6x4.txt',1,grid);
     
     % ObjHT = New_Mesh2D('File_mesh_exp1_32x32_3_leves.txt',3,grid);
     % ObjHT = New_Mesh2D('File_current_mesh.txt',4,grid);
     % ObjHT = New_Mesh2D('malla_4_refinada.txt',2,grid);
     % ObjHT = New_Mesh2D('File_mesh_general_prueba_3.txt',3,grid);
     % ObjHT = New_Mesh2D('Malla_refinada_L_centro_N_4.txt',2,grid);
     % ObjHT = New_Mesh2D('File_current_mesh.txt',2,grid);
      % ObjHT = New_Mesh2D('Malla_uniforme_base_32.txt',1,grid);

      %Esta es la malla test1
      % ObjHT = New_Mesh2D('File_mesh_general_prueba_3.txt',3,grid);
       % ObjHT = New_Mesh2D('File_current_mesh.txt',3,grid);

      % ObjHT = New_Mesh2D('Malla_test_particulas_8x8_3_level.txt',3,grid);
      % ObjHT = New_Mesh2D('File_current_mesh.txt',4,grid);

      %MALLA PARTICULAS 
      % ObjHT = New_Mesh2D('Malla_test2_16x16_4_level.txt',4,grid);
      % ObjHT = New_Mesh2D('Malla_test2_8x8_4_level.txt',4,grid);
      

      %MALAA TRST3 C.I GOTA
      % ObjHT = New_Mesh2D('Malla_test3_32x32_3_level.txt',3,grid);
      % ObjHT = New_Mesh2D('Malla_test3_32x32_3_level_new.txt',3,grid);
      % ObjHT = New_Mesh2D('File_current_mesh.txt',3,grid);
      % ObjHT = New_Mesh2D('Malla_test3_16x16_3_level.txt',3,grid);
     
      % ObjHT = New_Mesh2D('Base_uniform_mesh.txt',1,grid);

      %MALLA TEST2
       % ObjHT = New_Mesh2D('Malla_test1_3_nivel_8x8.txt',3,grid);

      %MALLA TEST4
       % ObjHT = New_Mesh2D('malla_test4_16x16_4levels.txt',4,grid);
      % ObjHT = New_Mesh2D('File_current_mesh.txt',4,grid);
     
     % ObjHT = New_Mesh2D('Malla_refinada_16_base_3_niveles_sin_rin_frontera.txt',3,grid.N,grid.a1, grid.b1);

     %MALLA LAGRANDE TEST1
         % ObjHT = New_Mesh2D('malla_test_1_lagrange_8x8_3_level.txt',3,grid);
         % ObjHT = New_Mesh2D('malla_test_1_lagrange_8x8_3_level_correcta_anidada.txt',3,grid);
        ObjHT = New_Mesh2D('malla_test_1_lagrange_16x16_3_level.txt',3,grid);
        % ObjHT = New_Mesh2D('malla_test_1_lagrange_32x32_3_level.txt',3,grid);

      %MALLA NXM DISTINTOS
      % ObjHT = New_Mesh2D('Malla_refinada_L_centro_6x4.txt',2,grid);
      % ObjHT = New_Mesh2D('Mesh_base_uniforme_6x4.txt',1,grid);
      
        
     else
         % ObjHT = New_Mesh2D('File_current_mesh.txt',2,grid);
         ObjHT = New_Mesh2D('File_current_mesh.txt',3,grid);
     end
     
     time_create_object = toc;
     disp(['time_create_object: ' num2str(time_create_object) ' segundos']);


     ObjHT.set_id_celulas_for_level
     ObjHT.set_enumeration_vector

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

     %Create Matriz and rhs used enumeration vector
     %[A,B] = discretizeDFPoissonEnumeration3(grid,pde,ObjHT); % Matrix used HT
     % [A,B] = discretizeDFPoissonEnumeration4_vector_vecinas(grid,pde,ObjHT);  %Matriz used vector and calculate neighbor
     % [A,B] = discretizeDFPoisson_Cond_Frontera_2_enayo_busca_de_error(grid,pde,ObjHT);
     % [A,B] = discretizeDFPoisson_busca_de_error(grid,pde,ObjHT);   %discretizacion basica
     

     %Versiones finales
     tic;
     % [A,B] = discretizeDFPoisson_Cond_Frontera(grid,pde,ObjHT);
     % [A,B] = discretizeDFPoisson_Cond_Frontera_Neumann(grid,pde,ObjHT);

     % [A,B] = discretizeDFPoisson_Inter_cuadratica_2(grid,pde,ObjHT); %Esta versión es la correta
     [A,B] = discretizeDFPoisson_MLS(grid,pde,ObjHT);
     %Input grid, pde and HT.
     time_create_matrix_vector = toc;
     disp(['time_create_matrix_vector: ' num2str(time_create_matrix_vector) ' segundos']);
    

    %Plot of sparse matriz
    % spy(A,'r*',7); 

    position_nulas = length(find(A==0));
    num_cell_mesh = ObjHT.get_Number_cells_Mesh;
    Entradas_Matriz = num_cell_mesh*num_cell_mesh;
    position_No_nulas = Entradas_Matriz - position_nulas;

    fprintf(fileID1,'%d &\t %d &\t %d \n',Entradas_Matriz,position_nulas,position_No_nulas);

    tasa_ocupacion_NO_nula(RepIni) = position_No_nulas/Entradas_Matriz;
    disp(tasa_ocupacion_NO_nula(RepIni));
    tasa_ocupacion_nula(RepIni) = position_nulas/Entradas_Matriz;


    counter_cell_level = ObjHT.getCell_counter_level;
    niveles_malla = length(counter_cell_level) - 1;

    for iter = 1: niveles_malla
        level_iter = counter_cell_level(iter);
        fprintf(fileID2,'%d &\t',level_iter);
    end
    fprintf(fileID2,'\n');



    %calculate enumeration vector and HT
    Vec_enumeation = ObjHT.get_enumeration_vector();   %Vector whit the cells
    % HT = ObjHT.get_HT;

    %Cell number
    num_incog = length(Vec_enumeation);

    
    % % % vec_collisions = ObjHT.get_number_collisions();
    % % % max_collisions = max(vec_collisions);
    % % % fprintf(fileID1,'%.d\n',max_collisions);
    % % % fprintf(fileID2,'%.d\n',max_collisions);

    % ObjHT.Refinar_malla();
    % ObjHT.Genera_archivo_malla_actual();

    %Updata merf file
    if RepIni ~= Rep_Final
        ObjHT.Refine_Mesh;
        ObjHT.Generate_current_mesh_file;
        Edite_level_of_file('File_current_mesh.txt');
    end
    
 end
 fprintf(fileID1,'\n\n');

fprintf(fileID1,'Vector tasa de ocupacion nula \n');
 fprintf(fileID1,'%.5f \t',tasa_ocupacion_nula);
 fprintf(fileID1,'\n\n');

 fprintf(fileID1,'Vector tasa de ocupacion NO nula \n');
 fprintf(fileID1,'%.5f \t',tasa_ocupacion_NO_nula);

 fprintf(fileID1,'\n');
 fprintf(fileID2,'\n');
% fclose(fileensayo);
fclose(fileID1);
fclose(fileID2);

