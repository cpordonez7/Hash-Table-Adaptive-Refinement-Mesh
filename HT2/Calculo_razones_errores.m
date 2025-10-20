
%Calculo de las razones de errores

fileErrores = fopen('Data_errores.txt','r');
                %-------------LEER ARCHIVO------------------
numData = fscanf(fileErrores,'%d\n',[1,1]); 
Errores = fscanf(fileErrores,'%e %e %e\n',[numData 1]);
% disp(Errores)

Datos_ = numData/3;


fclose(fileErrores);

fileID1 = fopen('Razones_errores.txt','w');

 %Calculo de las razones 
 fprintf(fileID1,'Razon Norma_2 \t Razon Norma_h \t Razon Norma_Inf \n');
 % fprintf(fileID2,'Razon Norma_2 \t Razon Norma_h \t Razon Norma_Inf \n');
 Razon1 = zeros(Datos_-1,1);
 Razon2 = zeros(Datos_-1,1);
 Razon3 = zeros(Datos_-1,1);

 for iter=0:Datos_-2
     % disp(iter)
     Razon1(iter+1) = Errores(3*iter+1,1)/Errores(3*(iter+1)+1,1);
     % disp(Errores(iter,1))
     % disp(Errores(iter,2))
     Razon2(iter+1) = Errores(3*iter+2,1)/Errores(3*(iter+1)+2,1);
     Razon3(iter+1) = Errores(3*iter+3,1)/Errores(3*(iter+1)+3,1);
     % fprintf(fileID1,'%.4f \t  %.4f \t %.4f \n',Razon1, Razon2, Razon3);
     % fprintf(fileID2,'%.4e \t  %.4e \t %.4e \n',Razon1, Razon2, Razon3);
 end

 for iter=1:Datos_-1
     fprintf(fileID1,'%.4f \n',Razon1(iter));
     % fprintf(fileID2,'%.4e \t  %.4e \t %.4e \n',Razon1, Razon2, Razon3);
 end
fprintf(fileID1,'\n');
 for iter=1:Datos_-1
     fprintf(fileID1,'%.4f \n',Razon2(iter));
     % fprintf(fileID2,'%.4e \t  %.4e \t %.4e \n',Razon1, Razon2, Razon3);
 end
fprintf(fileID1,'\n');
 for iter=1:Datos_-1
     fprintf(fileID1,'%.4f \n',Razon3(iter));
     % fprintf(fileID2,'%.4e \t  %.4e \t %.4e \n',Razon1, Razon2, Razon3);
 end
 

 fclose(fileID1);