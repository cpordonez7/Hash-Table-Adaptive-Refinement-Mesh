function [] = Create_uniform_mesh(N,M,a1,b1,a2,b2)
    %This function creates a uniform grid, save the points in a file
    %Input: 
    %N, number of partitions on the x-axis.
    %M, number of partitions on the y-axis
    %(a1,b1) endpoints of the domai in x.
    %(a2,b2) endpoints of the domai in y.

    fileID = fopen('Base_uniform_mesh.txt','w');

    fprintf(fileID,'%d \n',N*M);   % First data: Number of point from grid

    %step size
    h = (b1-a1)/N;
    k = (b2-a2)/M;
   
    %fix j and advance in i
    for j=1:M
        y = a2 + (j-1)*k;
        for i=1:N
            x = a1 + (i-1)*h;
            fprintf(fileID,'%.12f\t%.12f\t%d\n',x,y,0);
            % fprintf('%.12f\t%.12f\t%d\n',x,y,0);

        end
    end
end

