   
%Matriz de ejemplo
matriz = [4, 2, 3; 
          1, 6, 5; 
          7, 8, 9];

% % Fila de referencia (por ejemplo, la primera fila)
% fila_referencia = matriz(1, :);
% 
% % Obtener los índices ordenados de las columnas
% [~, indices_ordenados] = sort(fila_referencia);
% 
% % Reordenar las columnas de la matriz
% matriz_ordenada = matriz(:, indices_ordenados);
% 
% % Mostrar la matriz original y la matriz ordenada
% disp('Matriz original:');
% disp(matriz);
% 
% disp('Matriz ordenada por columnas con referencia a la primera fila:');
% disp(matriz_ordenada);

% Lista de ejemplo
% lista = [3, 1, 4, 1, 5, 9, 2, 6];
lista = [0.0884 0.0884 0.2500 0.2500 0.0884 0.0884 0.3536 0.0884 0.3536 0.3536];

% Ordenar la lista y obtener los índices originales
[lista_ordenada, indices] = sort(lista);

% Mostrar la lista ordenada y los índices originales
disp('Lista ordenada:');
disp(lista_ordenada);

disp('Índices originales de los elementos:');
disp(indices);
