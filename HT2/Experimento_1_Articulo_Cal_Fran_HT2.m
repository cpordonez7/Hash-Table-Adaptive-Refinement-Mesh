%Experimento 1  HT2
%Estructura malla , Dominio, (a1,b1)x(a2,b2)
 clear
 grid.a1 = 0.0;
 grid.a2 = 0.0;
 grid.b1 = 1.0;
 grid.b2 = 1.0;

 %Difference extremes in the domain
 Dif_dom_x = grid.b1 - grid.a1;
 Dif_dom_y = grid.b2 - grid.a2;

 %Number of partitions of the base mesh
 grid.N = 32;
 grid.M = 32;
 Num_level_refin = 5;  %number of levels to refine (l) (l>=0)
 %%profile on
 %Create the uniform mesh
 Create_uniform_mesh(grid.N,grid.M,grid.a1,grid.b1,grid.a2,grid.b2);  %N:particiones y extremos del dominio en un eje , N=M

 tic;
 % ObjHT = New_Mesh2D('Malla_uniforme_base_32.txt',Num_level_refin,grid); %archivo malla uniforme, numero de niveles, dominio
 ObjHT = New_Mesh2D('Base_uniform_mesh.txt',Num_level_refin,grid);
 time_create_object = toc;
 disp(['time_create_object: ' num2str(time_create_object) ' segundos']);

 %Time
 tic;
 %Loop that goes through the mesh while there is refinement  
 for iter=1:Num_level_refin

    %get the enumeration vector and HT
    Vec_enumeation = ObjHT.get_enumeration_vector();   %Vector whit the cells
    % HT = ObjHT.get_HT;

    %Cell number
    num_incog = length(Vec_enumeation);
    % disp(num_incog)

    %Vectors, coord center cell
    x_c = zeros(num_incog,1);
    y_c = zeros(num_incog,1);

    aprox = zeros(num_incog,1);

    %function validation 
    for j=1:(num_incog)
        cell_j = ObjHT.getCell_id(j); 
        ind_i = cell_j.i;
        ind_j = cell_j.j;
        l = cell_j.l;
        Id = cell_j.id;

        %Domaint data [a1,b1]x[a2,b2] for each cell in her nevel
        h_c = (Dif_dom_x)/(grid.N*2^l);
        k_c = (Dif_dom_y)/(grid.M*2^l);
        
        coordx = grid.a1 + ind_i*h_c;    % x value of the curreent cell
        coordy = grid.a2 + ind_j*k_c;    % y value of the curreent cell
     
        
        %Coordinates of the center of the cell
        x_cen = coordx + h_c*0.5;
        y_cen = coordy + k_c*0.5;

        x_c(j) = x_cen;
        y_c(j) = y_cen;

        % fprintf('%.10f \t %.10f \n',x_c(iter),y_c(iter));

        %Phi function 
        z = phi(x_cen,y_cen);
        % fprintf('%.10f \t %.10f \t %.10f \n',x_c(iter),y_c(iter),z);
        aprox(j) = z;
        % fprintf('%.10f \t %.10f \t %.10f \n',x_c(iter),y_c(iter),aprox(iter));

        %conditions functions refinement
        left_cond = (-0.999+0.05*l);
        right_cond = (0.999-0.05*l);

        %Refine
        if left_cond < z && z < right_cond 
            ObjHT.RefineCell(cell_j);
        end
    end  

    %Update the vector with the new inserted cells
    ObjHT.set_enumeration_vector;
 end
 time_refinement = toc;
 disp(['time_refinement: ' num2str(time_refinement) ' segundos']); 

 %%profile off

 % for m=1:num_incog
 %     fprintf('%.10f \t %.10f \t %.10f \n',x_c(m),y_c(m),aprox(m));
 % end

 % Grafica
% % % scatter3(x_c, y_c, aprox, 60, aprox, 'filled');
% % % colorbar;
% % % xlabel('Eje X');
% % % ylabel('Eje Y');
% % % zlabel('Eje Z');
% title('Mapa de calor en 3D');

% % % ObjHT.graph_Mesh2D();

 % Z1 = reshape(Sol_ex,grid.N,grid.M); %Solucion exacta
 % Z2=reshape(aprox,grid.N,grid.M);   %pasa de vector a matriz, Aproximacion


 %Plot solución   
 % [X,Y] = meshgrid(x,y);
 % contour(X,Y,Z2);
 % surf(X,Y,Z2);
 % surf(X,Y,Z2);
 % plot3(X,Y,Z2)
 % meshc(X,Y,Z2)
 % pcolor(X,Y,Z2)

 % plot3(x,y,Aprox,'MarkerSize',10)
%Con la función scatter se grafica
%Function phi original
function z = phi(x,y)
    xc = 0.5;
    yc = 0.5;
    d = sqrt((x-xc)^2+(y-yc)^2);

    z = tanh(75*(0.25-d));
end

% function z = phi(x,y)
%     xc = 1.0/13.0;
%     yc = 1.0/6.0;
%     d = sqrt((x-xc)^2+(y-yc)^2);
% 
%     z = 2*tanh(45*(0.5-d));
% end





