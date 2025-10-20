%%Cells with i points 
%Read two files, MLS1 and ML2 

%Ex: MESH TEST4
 % file1= fopen('Cells_con_i_puntos_MLS1_TEST4_32.txt','r');
 % file2= fopen('Cells_con_i_puntos_MLS2_TEST4_32.txt','r');

 file1= fopen('File_cell_with_i_points_MLS1.txt','r'); %File MLS1
 file2= fopen('File_cell_with_i_points_MLS2.txt','r'); %File MLS2

 %Fist data: Max number
 num_max1 = fscanf(file1,'%d \n',1);
 num_min1 = fscanf(file1,'%d \n',1);

 %Second data: Min number
 num_max2 = fscanf(file2,'%d \n',1);
 num_min2 = fscanf(file2,'%d \n',1);

 %Cells for i-points
 cell_MLS1 = fscanf(file1,'%d');  
 cell_MLS1 = cell_MLS1(:)'; % to row vector

 cell_MLS2 = fscanf(file2,'%d');  
 cell_MLS2 = cell_MLS2(:)'; % to row vector

 if num_min1<num_min2
     dif1= num_min2 - num_min1;
     cell_MLS2 = [zeros(1,dif1), cell_MLS2(num_min2:num_max2)];
     min = num_min1;
 else
     dif1= num_min1 - num_min2;
     cell_MLS1 = [zeros(1,dif1), cell_MLS1(num_min1:num_max1)];
     min = num_min2;
 end

 if num_max1<num_max2
     dif2= num_max2 - num_max1;
     cell_MLS1 = [cell_MLS1(num_min1:num_max1), zeros(1,dif2)];
     max = num_max2;
 else
     dif2= num_max1 - num_max2;
     cell_MLS2 = [cell_MLS2(num_min2:num_max2), zeros(1,dif2)];
     max = num_max1;
 end

i_points = min:max;
bar(i_points, [cell_MLS1',cell_MLS2']);


% Bars graph
legend('MLS1','MLS2');
xlabel('i-points'); ylabel('Number of cells with i points');
% %title('----');

fclose(file1);
fclose(file2);