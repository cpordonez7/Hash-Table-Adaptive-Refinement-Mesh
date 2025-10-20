%Problema de las particulas sin funcion RK4 y refinamiento de malla
%Este archivo refina y engrosa tiene error se quiere que engrose todas sino
%las que se mueven 
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
Flag_matrix = zeros(2,Num_particles);   %Aux matrix: This has values [0/1,level] for each particles. 1--> The particles are iquals


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

iter = 1;
%Update particles
 x1 = solution(iter,1:2:2*Num_particles);
 y1 = solution(iter,2:2:2*Num_particles);
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



%%Cycle counter, starts in iter=2 because in t=1 are the initial conditions
iter_t = 2;
%Do, while the end time is reached or maximun iterations are reached
while (iter_t <= length(t))
% while (iter_t <= 2)
    % fprintf('tiiiempo: %d\n',iter_t);
    % disp(Flag_matrix);
    % disp(Id_celulas_Refinadas);
    %in each time step works all the particles
    %Particles at t = iter
    x1 = solution(iter_t,1:2:2*Num_particles);
    y1 = solution(iter_t,2:2:2*Num_particles);

    plot(x1,y1,'.r','MarkerSize',5)

    %PROCESO DE ENGROSE PARA CADA PASO DE TIEMPO
    %calculate enumeration vector and HT
    HT = ObjMesh.get_HT();

    l_cel = Num_leves_refine ;
    while l_cel > 0
        for it_part=1:Num_particles
        % for it_part=1:1
            if Flag_matrix(1,it_part) ~= 1  
                % disp('hola');
                x = x1(it_part);
                y = y1(it_part);
                l = l_cel-1;

                %calculate of cell indices
                ind_i = grid.N*2^(l)*(x-grid.a1)/(grid.b1-grid.a1);
                ind_j = grid.M*2^(l)*(y-grid.a2)/(grid.b2-grid.a2);

                ind_i = floor(ind_i);
                ind_j = floor(ind_j);

                %Look for the cell
                %Key calulation
                cell_particle = create_cell(ind_i,ind_j,l);
                key1 = ObjMesh.key_Calculation(cell_particle);

                cell_particle_fina = create_cell(ind_i,ind_j,l_cel);
                key2 = ObjMesh.key_Calculation(cell_particle_fina);

                p = ObjMesh.hast_function(key2);
         
                % fprintf('i:%d \t j:%d \t l:%d key:%d\n',ind_i,ind_j,l,key);
               
                size_HT_p = size(HT{p+1},2);
                if (size_HT_p ~=0)
                    %True or False
                    flag = HT{p+1}.searchNode(key2);
                    % disp(flag)

                    if strcmp(flag, 'True')
                        cell_particle_fina = HT{p+1}.get_cell_Nodo(key2);
                        l2 = cell_particle_fina.l;
                    end
                end
                
                key_cel = Key_celulas_Refinadas(l_cel,2*it_part-1);
                l_part = Key_celulas_Refinadas(l_cel,2*it_part);

                if key_cel ~= key1 || l ~=l_part
                        p_old = ObjMesh.hast_function(key_cel);
                        size_HT_p_old = size(HT{p_old+1},2);
                        if l_part ~= 0 && l ~= 0 && (size_HT_p_old ~=0)
                            flag_old = HT{p_old+1}.searchNode(key_cel);
                            if strcmp(flag_old, 'True')
                                cell_old = HT{p_old+1}.get_cell_Nodo(key_cel);
                                if cell_old.l ~= 0
                                    ObjMesh.ThickenCell(cell_old);
                                    Key_celulas_Refinadas(l_cel,2*it_part-1) = 0;
                                    Key_celulas_Refinadas(l_cel,2*it_part) = 0;
                                end
                            end
                            % disp('holaaaaa')
                            
                        end
                elseif key_cel == key1 && l == l_part
                    Flag_matrix(1,it_part) = 1;
                    Flag_matrix(2,it_part) = l_cel;
                end

            end  %If flag
        end
        l_cel = l_cel - 1;
    end
    % disp(Id_celulas_Refinadas);

    %%Refine
    %Loop that goes through the particle vector and refine from where the
    %flag indicates
    iter_key = 1;  %Counter for mark the particles 
    for iter=1:Num_particles
        l = Flag_matrix(2,iter);
        % disp(l)

        while l < Num_leves_refine
            x = x1(iter);
            y = y1(iter);

            %calculate of cell indices
            ind_i = grid.N*2^(l)*(x-grid.a1)/(grid.b1-grid.a1);
            ind_j = grid.M*2^(l)*(y-grid.a2)/(grid.b2-grid.a2);

            ind_i = floor(ind_i);
            ind_j = floor(ind_j);

            % fprintf('parti: %d \t x: %.5f \t y: %.5f \t l: %d\n',iter,x,y,l);
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

    % disp(Id_celulas_Refinadas);

    ObjMesh.graph_Mesh2D();
    pause(0.00001)
    hold off

    %Initialize Flag in zeros
    Flag_matrix = zeros(2,Num_particles);

    %iter while, iter that runs through all the times in the interval
    iter_t = iter_t + 1;  

    % ObjHT.graficar();
    % pause(0.001)
    % hold off

end

time_refine_thicken = toc;
disp(['time_refine_thicken: ' num2str(time_refine_thicken) ' segundos']);

profile off









