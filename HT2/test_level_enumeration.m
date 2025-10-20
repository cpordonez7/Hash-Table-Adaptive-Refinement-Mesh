%Para ennumerar de abajo hacia arriba y de izq a derecha, se varía primero
%i, luego j y nuevamente i

%Para ennumerar de isq a derecha y de abajo hacia arriba , se varía primero
%j, luego i y nuevamente j


% Arreglo bidimensional de ejemplo
arr = [0 2 3 3 1 0 3; 0 2 2 0 0 2 1];
arrtr = arr';
% disp(arr)
% Transpone t ordena para i
[arr_1, indices1] = sortrows(arr', 1);

%Ordena para j
[arr_2, indices2] = sortrows(arr_1, 2);

%Ordena para i
[arr_3, indices3] = sortrows(arr_2, 1)


% % % % lista_enumaration = cell(3,1);
% % % % 
% % % % lista_enumaration{1} = [1 1; 3 2; 0 2];
% % % % lista_enumaration{2} = [0 2 3 3 1 0 3; 0 2 2 0 0 2 1];
% % % % 
% % % % 
% % % % 
% % % % arrtr = lista_enumaration{1}';
% % % % % disp(arr)
% % % % % Transpone t ordena para i
% % % % [arr_1, indices1] = sortrows(arrtr, 1);
% % % % 
% % % % %Ordena para j
% % % % [arr_2, indices2] = sortrows(arr_1, 2);
% % % % 
% % % % %Ordena para i
% % % % [arr_3, indices3] = sortrows(arr_2, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% lista_enumaration = cell(1,1);
% 
% lista_enumaration{1} = zeros(2,2);
% lista_enumaration{2} = zeros(2,7);




