function [] = Graph_Mesh2D_file(HTObj,file_graph)
    %This function plots the grid with points in the entered file.

    %Colors vector
    Color = {'k','r','b','c','g'};
     % Color = {'k','k','k','k'};

    %Read the file
    file= fopen(file_graph,'r');
    %-------------READ FILE------------------
    numNodes = fscanf(file,'%d\n',[1,1]);  %Number of points, first data in the file
    Grid_Nodos = fscanf(file,'%e %e %d\n',[3 numNodes]);  %Array 3xNumNodes, save (x,y,l)

    % numNodos = HTObj.getNumNodos();  %Nume de nodos en la malla
    % Grid_Nodos = HTObj.Nodosmalla; 

    %difference between the domain extremes
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;

    hold on

    %Iterate through the points file
    for j=1:numNodes
        x = Grid_Nodos(1,j);  %(x,y,l) of each point
        y = Grid_Nodos(2,j);
        l = Grid_Nodos(3,j);

        %Draw the line according to the level
        h = Dif_dom_x/(HTObj.N*2^l);
        k = Dif_dom_y/(HTObj.M*2^l);

        line([x x+h],[y y],'color',Color{l+1},'LineWidth',1.0) %Linea de (x,y) a (x+h,y)
        line([x x],[y y+k],'color',Color{l+1},'LineWidth',1.0) %Linea de (x,y) a (x,y+k)

    end
    axis('off');  %No Muestra los ejes

end

