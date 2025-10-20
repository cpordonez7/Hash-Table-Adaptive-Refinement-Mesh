function [] = graph_Mesh2D(HTObj)
    %This function graph the mesh that the HT has.

     %Colors vector
    % Color = {'k','r','b','g','c'};
    Color = {'k','k','k','k','k'};
    %Draw the conrourn of the mesh
    line([HTObj.a1 HTObj.b1],[HTObj.a2 HTObj.a2],'color',Color{1},'LineWidth',1.0) %Linea de (0,0) to (1,0)
    line([HTObj.a1 HTObj.a1],[HTObj.a2 HTObj.b2],'color',Color{1},'LineWidth',1.0) %Linea de (0,0) to (0,1)
    line([HTObj.b1 HTObj.b1],[HTObj.a2 HTObj.b2],'color',Color{1},'LineWidth',1.0) %Linea de (1,0) to (1,1)
    line([HTObj.a1 HTObj.b1],[HTObj.b2 HTObj.b2],'color',Color{1},'LineWidth',1.0) %Linea de (0,1) to (1,1)

    %Difference domial extrem
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;

    hold on

    for p=1:HTObj.size_HT
    % for p=1:1
        %Check if the linked list is empty

        % if size(HTObj.HT{p},2) == 0
        %     fprintf('The list in the position %d is empty \n',p);
        if size(HTObj.HT{p},2) ~= 0 && HTObj.HT{p}.head.counter ~= 0
            HTObj.HT{p}.current = HTObj.HT{p}.first_node;
            while (HTObj.HT{p}.current ~= HTObj.HT{p}.last_node)
                %date fron cell in the p linked list 
                ind_i = HTObj.HT{p}.current.cell_struct.i;
                ind_j = HTObj.HT{p}.current.cell_struct.j;
                l = HTObj.HT{p}.current.cell_struct.l;
                % key = HTObj.HT{p}.current.cell_struct.key;

                % fprintf('p: %d \t key: %d \t i: %d \t j: %d \t l: %d\n',p, key, ind_i, ind_j, l);

                h = (Dif_dom_x)/(HTObj.N*2^(l));
                k = (Dif_dom_y)/(HTObj.M*2^(l));

                % fprintf('Dif_dom_X: %f\t Dif_dom_Y: %f\t h: %.5f \t k: %.5f \n',Dif_dom_x, Dif_dom_y, h, k);

                x = HTObj.a1 + ind_i*h;
                % x = HTObj.a1 + ind_i*(HTObj.b1-HTObj.a1)/(HTObj.N*2^(l));
                y = HTObj.a2 + ind_j*k;
                
                % fprintf('x: %f \t y: %f \n',x,y);

                % fprintf('x:%lf \t y:%lf\n ',x,y);

                %only black color
                line([x x+h],[y y],'color',Color{1},'LineWidth',1.0) %Linea de (x,y) a (x+h,y)
                line([x x],[y y+k],'color',Color{1},'LineWidth',1.0) %Linea de (x,y) a (x,y+k)
                % line([x x+h],[y y],'color',Color{l+1},'LineWidth',1.0) %Linea de (x,y) a (x+h,y)
                % line([x x],[y y+k],'color',Color{l+1},'LineWidth',1.0) %Linea de (x,y) a (x,y+k)

                HTObj.HT{p}.current = HTObj.HT{p}.current.next;
            end

            if (HTObj.HT{p}.current == HTObj.HT{p}.last_node)
                %date fron last cell in the p linked list 
                ind_i = HTObj.HT{p}.current.cell_struct.i;
                ind_j = HTObj.HT{p}.current.cell_struct.j;
                l = HTObj.HT{p}.current.cell_struct.l;
                % key = HTObj.HT{p}.current.cell_struct.key;
                % fprintf('p: %d \t key: %d \t i: %d \t j: %d \t l: %d\n',p, key, ind_i, ind_j, l);

                h = (Dif_dom_x)/(HTObj.N*2^(l));
                k = (Dif_dom_y)/(HTObj.M*2^(l));

                % fprintf('Dif_dom_X: %f\t Dif_dom_Y: %f\t h: %.5f \t k: %.5f \n',Dif_dom_x, Dif_dom_y, h, k);

                x = HTObj.a1 + ind_i*h;
                y = HTObj.a2 + ind_j*k;

                % fprintf('x: %f \t y: %f \n',x,y);

                %only black color
                line([x x+h],[y y],'color',Color{1},'LineWidth',1.0) %Linea de (x,y) a (x+h,y)
                line([x x],[y y+k],'color',Color{1},'LineWidth',1.0) %Linea de (x,y) a (x,y+k)

                % line([x x+h],[y y],'color',Color{l+1},'LineWidth',1.0) %Linea de (x,y) a (x+h,y)
                % line([x x],[y y+k],'color',Color{l+1},'LineWidth',1.0) %Linea de (x,y) a (x,y+k)
            end

        end
    end

    axis('off');  %Shows the axes
end
