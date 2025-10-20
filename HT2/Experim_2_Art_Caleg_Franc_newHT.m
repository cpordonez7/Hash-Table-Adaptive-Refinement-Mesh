%Problema de las particulas sin funcion RK4 y refinamiento de malla
%Este archivo refina y engrosa
clear

%Mesh data [-0.5,0.5]x[-0.5,0.5]
grid.a1 = -0.5;
grid.b1 = 0.5;
grid.a2 = -0.5;
grid.b2 = 0.5;
grid.N = 16;
grid.M = 16;

%Difference domial extrem
Dif_dom_x = grid.b1 - grid.a1;
Dif_dom_y = grid.b2 - grid.a2;

%number of levels to refine (3--> Refine three times, with four levels in the mesh)
%Numero de veces a refinar 2--> tres niveles
Num_leves_refine =3;            %max level in the mesh (L) (Number of max levels)(l+1)

%time interval
t0 = 0;
tf = 0.6;
time_interval = [t0,tf];

%%Initial particle creation (Initial conditions)
rng(0)                                                  %seed(random)
Num_particles = 1000;
y1 = zeros(Num_particles,1);                           %Coordinates points in Y (all zeros)
x1 = -0.45 + (0.45-(-0.45)).*(rand(Num_particles,1));  %Coordinates points in X (random in (-0.45,0.45))
y0 = zeros(2*Num_particles,1);

%Parameters
beta = 500.0;
mu = 10;

BETAC = beta/(2*pi);
MUC = 4*mu;

%Initial condition vector (x0,y0,x1,y1,x2,y2,x3,y3,...,xn,yn) 
for iter=1:Num_particles
    y0(2*iter-1) = x1(iter);
    y0(2*iter) = y1(iter);
end

%Plot particles
% % % plot(y0(1:2:2*Num_particles),y0(2:2:2*Num_particles),'.b','MarkerSize',7)
% axis([-0.5 0.5 -0.5 0.5])
% pause(0.3)
% hold on

%Create uniform mesh file ('Malla_uniforme_base.txt')
Create_uniform_mesh(grid.N,grid.M,grid.a1,grid.b1,grid.a2,grid.b2);

profile on
tic;
%Mesh creation, max level:4, N:16 M:16, a1:-0.5, b1=0.5 a2:-0.5, b2=0.5
ObjMesh = New_Mesh2D('Base_uniform_mesh.txt',Num_leves_refine,grid);
% ObjHT.graficar();
time_create_object = toc;
 disp(['time_create_object: ' num2str(time_create_object) ' segundos']);

Key_celulas_Refinadas = zeros(Num_leves_refine,2*Num_particles);  %key and level of the refined cell 
iter_key = 1;

% % % %Loop that goes through the particle vector and if it is in the refined mesh
% % % for iter=1:Num_particles
% % %     % disp(iter)
% % %     l = 0;
% % %     x = x1(iter);
% % %     y = y1(iter);
% % % 
% % %     while l < Num_leves_refine 
% % %         % disp(l)
% % %         HT = ObjMesh.get_HT();
% % %         %calculate of cell indices
% % %         ind_i = grid.N*2^(l)*(x-grid.a1)/(Dif_dom_x);
% % %         ind_j = grid.M*2^(l)*(y-grid.a2)/(Dif_dom_y);
% % % 
% % %         ind_i = floor(ind_i);
% % %         ind_j = floor(ind_j);
% % % 
% % %         % %Calculate of the vertices of the cell with indices (ind_i,ind_j) 
% % %         % x = grid.a1 + ind_i*(Dif_dom_x)/(grid.N*2^(l));
% % %         % y = grid.a2 + ind_j*(Dif_dom_y)/(grid.M*2^(l));
% % % 
% % %         %Key calulation
% % %         cell_particle = create_cell(ind_i,ind_j,l);
% % %         key = ObjMesh.key_Calculation(cell_particle);
% % %         p = ObjMesh.hast_function(key);
% % %         % fprintf('i:%d \t j:%d \t l:%d key:%d\n',ind_i,ind_j,l,key);
% % % 
% % %         size_HT_p = size(HT{p+1},2);
% % %         if (size_HT_p ~=0)
% % %             %True or False
% % %             flag = HT{p+1}.searchNode(key);
% % %             % disp(flag)
% % % 
% % %             if strcmp(flag, 'True')
% % %                 cell_particle = HT{p+1}.get_cell_Nodo(key);
% % %                 l2 = cell_particle.l;
% % %             end
% % % 
% % %             %Check if the cell is in the HT. If the cell exists, so refine.
% % %             if strcmp(flag, 'True') && l == l2
% % %                 %Refine and save the key and level of the cell refine
% % %                 Key_celulas_Refinadas(l+1,iter_key) = cell_particle.key; 
% % %                 Key_celulas_Refinadas(l+1,iter_key+1) = l2; 
% % %                 ObjMesh.RefineCell(cell_particle);
% % %             end
% % %         end
% % % 
% % %         l = l+1;
% % % 
% % %     end
% % %     iter_key = iter_key + 2;
% % % 
% % % end
% % % % disp(Key_celulas_Refinadas)
% % % 
% % %  ObjMesh.graph_Mesh2D();
% % %  pause(0.03)
% % %  % 
% % %  hold off

 % % % %calculate enumeration vector and HT
 % % % % ObjMesh.set_enumeration_vector;
 % % % % Vec_enumeation = ObjMesh.get_enumeration_vector;
 % % % 
 % % % %THICKEN PROCESS
 % % % HT = ObjMesh.get_HT();
 % % % 
 % % % l_cel = Num_leves_refine;
 % % % while l_cel > 0
 % % %     it_part = 1;
 % % %     while it_part <= 2*Num_particles
 % % %     % while it_part <= 2
 % % %         %Cell data that was refinate
 % % %         key_cel = Key_celulas_Refinadas(l_cel,it_part);
 % % %         l_part = Key_celulas_Refinadas(l_cel,it_part+1);
 % % %         if key_cel ~= 0.0 || l_part ~= 0
 % % %             p = ObjMesh.hast_function(key_cel);
 % % %             % fprintf('key:%d \t p:%d \t l_cel:%d \n',key_cel,p,l_cel);
 % % %             size_HT_p = size(HT{p+1},2);
 % % %             if (size_HT_p ~=0)
 % % %                 %True or False
 % % %                 flag = HT{p+1}.searchNode(key_cel);
 % % %                 % disp(flag)
 % % % 
 % % %                if strcmp(flag, 'True')
 % % %                    cell_part = HT{p+1}.get_cell_Nodo(key_cel);
 % % %                    l2 = cell_part.l;
 % % %                    % disp(l2)
 % % %                end
 % % % 
 % % %                %Check if the cell is in the HT. If the cell exists, so Thinken.
 % % %                if strcmp(flag, 'True') && l_cel == l2
 % % %                    %Thinken 
 % % %                    ObjMesh.ThickenCell(cell_part);
 % % %                    % disp(cell_part)
 % % %                end
 % % %             end
 % % %         %case then the cell is in the zero level ans is the origen    
 % % %         elseif key_cel == 0.0 && l_part == 0.0 && l_cel == 1 
 % % %            % disp(it_part)
 % % %             p = ObjMesh.hast_function(key_cel);
 % % %             size_HT_p = size(HT{p+1},2);
 % % %             if (size_HT_p ~=0)
 % % %                 %True or False
 % % %                 flag = HT{p+1}.searchNode(key_cel);
 % % %                 % disp(flag)
 % % % 
 % % %                if strcmp(flag, 'True')
 % % %                    % disp(it_part)
 % % %                    cell_part = HT{p+1}.get_cell_Nodo(key_cel);
 % % %                    l2 = cell_part.l;
 % % %                end
 % % % 
 % % %                %Check if the cell is in the HT. If the cell exists, so Thinken.
 % % %                %If it exists, it means it comes from a refinement process; 
 % % %                %otherwise, it's the label for a cell that has already been refined. 
 % % %                if strcmp(flag, 'True') && l_cel == l2
 % % %                    %Thinken 
 % % %                    ObjMesh.ThickenCell(cell_part);
 % % %                end
 % % %             end
 % % %         end 
 % % % 
 % % %         it_part = it_part + 2;
 % % %     end
 % % %     l_cel = l_cel - 1;
 % % % end
 % ObjMesh.graph_Mesh2D(); %NO ME ESTÁ GRAFICANDO EL ENGROSE junto
    % % % pause(0.00001)
    % % % hold off
 % clear Id_celulas_Refinadas;      


%Solver SEDO Ode45
%Returns a array with the solution at each time step, t1 --> [x0 y0 x1 y1
%x2 y2 ... xn yn], t2--->[x0' y0' x1' y1' x2' y2' ... xn' yn']...

%Solver SEDO
tic;
[t, solution] = ode45(@(t,y0) SEDO_particles(t,y0,BETAC,MUC), time_interval, y0);
time_solver_system = toc;
disp(['time_solver_system: ' num2str(time_solver_system) ' segundos']);

tic;
disp(length(t))
iter = 1;
%Do, while the end time is reached or maximun iterations are reached
while (iter <= length(t))
% while (iter == 300)
    %in each time step works all the particles
    %RK4 vector shape

    %Update particles
    x1 = solution(iter,1:2:2*Num_particles);
    y1 = solution(iter,2:2:2*Num_particles);

    %Plot the particles at each time step
    plot(x1,y1,'.r','MarkerSize',5)
    % axis([grid.a1 grid.a2 grid.a1 grid.a2])
    % pause(0.003)
    % hold off
    % pause(0.003)


    %Create base mesh
    % ObjHT = HastTable_V3('Malla_uniforme_base.txt',4,16,-0.5, 0.5);
    % ObjHT.graficar();

    % Key_celulas_Refinadas = zeros(Num_leves_refine,2*Num_particles);  %key and level of the refined cell 
    iter_key = 1;

    %Loop that goes through the particle vector and if it is in the refined mesh
    for iter_j=1:Num_particles
        % disp(iter)
        l = 0;
        x = x1(iter_j);
        y = y1(iter_j);

        if (abs(x1(iter_j)) <= 0.5 && abs(y1(iter_j))<=0.5)

            while l < Num_leves_refine 
                % disp(l)
                HT = ObjMesh.get_HT();
                %calculate of cell indices
                ind_i = grid.N*2^(l)*(x-grid.a1)/(Dif_dom_x);
                ind_j = grid.M*2^(l)*(y-grid.a2)/(Dif_dom_y);

                ind_i = floor(ind_i);
                ind_j = floor(ind_j);

                % %Calculate of the vertices of the cell with indices (ind_i,ind_j) 
                % x = grid.a1 + ind_i*(Dif_dom_x)/(grid.N*2^(l));
                % y = grid.a2 + ind_j*(Dif_dom_y)/(grid.M*2^(l));

                %Key calulation
                cell_particle = create_cell(ind_i,ind_j,l);
                key = ObjMesh.key_Calculation(cell_particle);
                p = ObjMesh.hast_function(key);
                % fprintf('i:%d \t j:%d \t l:%d key:%d\n',ind_i,ind_j,l,key);

                size_HT_p = size(HT{p+1},2);
                if (size_HT_p ~=0)
                    %True or False
                    flag = HT{p+1}.searchNode(key);
                    % disp(flag)

                    if strcmp(flag, 'True')
                        cell_particle = HT{p+1}.get_cell_Nodo(key);
                        l2 = cell_particle.l;
                    end

                    %Check if the cell is in the HT. If the cell exists, so refine.
                    if strcmp(flag, 'True') && l == l2
                        %Refine and save the key and level of the cell refine
                        Key_celulas_Refinadas(l+1,iter_key) = cell_particle.key; 
                        Key_celulas_Refinadas(l+1,iter_key+1) = l2; 
                        ObjMesh.RefineCell(cell_particle);
                    end
                end

                l = l+1;

            end
            iter_key = iter_key + 2;

        end
    end

    ObjMesh.set_enumeration_vector();
    ObjMesh.graph_Mesh2D();
    pause(0.00001)
    hold off

    % if iter==300
    %     ObjMesh.Generate_current_mesh_file;
    %     % num_cell = ObjMesh.get_Number_cells_Mesh;
    %     % disp(num_cell)
    % end

    %THICKEN PROCESS OR EACH TIME STEP
    HT = ObjMesh.get_HT();

    l_cel = Num_leves_refine;
    while l_cel > 0
        it_part = 1;
        while it_part <= 2*Num_particles
        % while it_part <= 2
         %Cell data that was refinate
            key_cel = Key_celulas_Refinadas(l_cel,it_part);
            l_part = Key_celulas_Refinadas(l_cel,it_part+1);
            if key_cel ~= 0.0 || l_part ~= 0
                p = ObjMesh.hast_function(key_cel);
                % fprintf('key:%d \t p:%d \t l_cel:%d \n',key_cel,p,l_cel);
                size_HT_p = size(HT{p+1},2);
                if (size_HT_p ~=0)
                    %True or False
                    flag = HT{p+1}.searchNode(key_cel);
                    % disp(flag)

                    if strcmp(flag, 'True')
                        cell_part = HT{p+1}.get_cell_Nodo(key_cel);
                        l2 = cell_part.l;
                        % disp(l2)
                    end

                    %Check if the cell is in the HT. If the cell exists, so Thinken.
                    if strcmp(flag, 'True') && l_cel == l2
                        %Thinken 
                        %%% ObjMesh.ThickenCell(cell_part);
                        ObjMesh.ThickenCell_v_sin_enum(cell_part);
                        % % % ObjMesh.set_enumeration_vector();
                        % disp(cell_part)
                    end
                end
            elseif key_cel == 0.0 && l_part == 0.0 && l_cel == 1 %This is the case then the cell is un the level zero
                % disp(it_part)
                p = ObjMesh.hast_function(key_cel);
                size_HT_p = size(HT{p+1},2);
                if (size_HT_p ~=0)
                    %True or False
                    flag = HT{p+1}.searchNode(key_cel);
                    % disp(flag)

                    if strcmp(flag, 'True')
                        % disp(it_part)
                        cell_part = HT{p+1}.get_cell_Nodo(key_cel);
                        l2 = cell_part.l;
                    end

                    %Check if the cell is in the HT. If the cell exists, so Thinken.
                    %If it exists, it means it comes from a refinement process; 
                    %otherwise, it's the label for a cell that has already been refined. 
                    if strcmp(flag, 'True') && l_cel == l2
                        %Thinken 
                        % % % ObjMesh.ThickenCell(cell_part);
                        ObjMesh.ThickenCell_v_sin_enum(cell_part);
                        % % % ObjMesh.set_enumeration_vector();
                    end
                end
            end     
            it_part = it_part + 2;
        end
        l_cel = l_cel - 1;
    end

    %iter while, iter that runs through all the times in the interval
    iter = iter + 1;  

    % ObjHT.graficar();
    % pause(0.001)
    % hold off

end
time_refine_thicken = toc;
disp(['time_refine_thicken: ' num2str(time_refine_thicken) ' segundos']);

profile off









