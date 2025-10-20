function [] = Edite_level_of_file(file_input)
    %Edit file: subtract one from the level
    %Input: file

    file = fopen(file_input,'r');
    %-------------READ FILE------------------
    numElem = fscanf(file,'%d\n',[1,1]);  %Number of the points, first data of the file
    Data_mesh = fscanf(file,'%e %e %d\n',[3 numElem]);  %Array 3xNumElem, save (x,y,l)

    fileID = fopen('File_current_mesh.txt','w');

    fprintf(fileID,'%d\n',numElem);   % First data: Number of point from grid
   
    %iterate the data
    for j=1:numElem
        x = Data_mesh(1,j);
        y = Data_mesh(2,j);
        l = Data_mesh(3,j)-1;
        fprintf(fileID,'%.12f\t%.12f\t%d\n',x,y,l);
    end
end