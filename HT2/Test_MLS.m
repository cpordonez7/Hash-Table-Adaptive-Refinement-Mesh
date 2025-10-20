%%CEL 93
L=[0.0625 -0.0625 0.1875 0.015625 0.046875 0.078125 0.109375 -0.015625 0.140625  ;
    -0.4375 -0.4375  -0.4375 -0.359375  -0.3593750 -0.359375 -0.359375 -0.359375  -0.359375];

P = [ 0.015625000000000;
  -0.390625000000000];

% %%CEL 205
%  L=[-0.0625 -0.109375 -0.078125 -0.046875 -0.0156250 -0.1875 0.0625 -0.140625 0.015625;
%     0.4375 0.359375   0.359375 0.359375 0.359375 0.4375 0.4375 0.359375 0.359375];
% 
%  P =[-0.046875;
%    0.390625];

m = length(L);
W = zeros(m,m);
P_matrix = zeros(m,6);

for iter =1:m
    W(iter,iter)= 1.0/norm(L(:,iter)-P);
    P_matrix(iter,:) = base(L(:,iter));
end

WP = W*P_matrix;

fileID1 = fopen('Test_Matrices_MLS.txt','w');
fprintf(fileID1,'\t PUNTO A APROX \n\n');

fprintf(fileID1,'P = [%f %f]',P);
fprintf(fileID1,'\n');

fprintf(fileID1,'\t LISTA PUNTOS \n\n');
fprintf(fileID1,'L = [');
for iter=1:2
    fprintf(fileID1,'%f  ',L(iter,:));
    fprintf(fileID1,'\n');
end
fprintf(fileID1,']');

fprintf(fileID1,'\n');


fprintf(fileID1,'\t Matriz W \n\n');
for iter=1:m
    fprintf(fileID1,'%f  ',W(iter,:));
    fprintf(fileID1,'\n');
end
fprintf(fileID1,'\n');

fprintf(fileID1,'\t Matriz P \n\n');
for iter=1:m
    fprintf(fileID1,'%f  ',P_matrix(iter,:));
    fprintf(fileID1,'\n');
end
fprintf(fileID1,'\n');

fprintf(fileID1,'\t\t Matriz WP \n\n');
for iter=1:m
    fprintf(fileID1,'%f  ',WP(iter,:));
    fprintf(fileID1,'\n');
end
fprintf(fileID1,'\n');
fprintf(fileID1,'Rango WP: %d  ',rank(WP));

fprintf(fileID1,'\n');

fprintf(fileID1,'\t\t factor \n\n');
for iter=1:m
    fprintf(fileID1,'%f  ',WP(iter,1)/WP(iter,6));
    fprintf(fileID1,'\n');
end

fprintf(fileID1,'\n');
% Aplica la eliminación Gaussiana
fprintf(fileID1,'\t\t Eliminación Gaussiana \n\n');
rref_WP = rref(WP);
for iter=1:m
    fprintf(fileID1,'%f  ',rref_WP(iter,:));
    fprintf(fileID1,'\n');
end

fclose(fileID1);



% % % Crea una matriz con los vectores como columnas
% % A = [WP(:,1), WP(:,6)];
% % 
% % % Aplica la eliminación Gaussiana
% % rref_A = rref(A)
% % 
% % % Verifica si hay alguna fila nula en la parte de coeficientes
% % dependent = any(all(rref_A(:, 1:end-1) == 0, 2));
% % 
% % %1--> si , 0-->no
% % if dependent
% %     disp('Los vectores son linealmente dependientes.');
% % else
% %     disp('Los vectores son linealmente independientes.');
% % end







 %Aplica la descomposición QR
[Q, R] = qr(WP)

% Calcula la matriz de dependencia lineal
dependent_columns = [];
tolerance = 1e-10; % Define una tolerancia para considerar si un valor es cero

for iter = 1:min(size(WP))
    if abs(diag(R(iter,iter))) < tolerance % Si el elemento de la diagonal es aproximadamente cero
        dependent_columns = [dependent_columns, iter];
    end
end

disp("Las columnas dependientes son:");
disp(dependent_columns);

b = base(P);






function B = base(P1)
        %Bases
        x = P1(1);
        y = P1(2);
        B = [1 x y x*y x*x y*y];
end



