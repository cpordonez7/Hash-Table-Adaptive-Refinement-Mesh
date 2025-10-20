function [] = graph_enumeration(HTObj)
    %This function graph the enumeration of the cell, according to its
    %identification.

    %Difference domial extrem
    Dif_dom_x = HTObj.b1-HTObj.a1;
    Dif_dom_y = HTObj.b2-HTObj.a2;

    hold on
    % iter = 0;
    for p=1:HTObj.size_HT
    % for p=1:1
        %Check if the linked list is empty

        if size(HTObj.HT{p},2) ~= 0 && HTObj.HT{p}.head.counter ~= 0
            
            HTObj.HT{p}.current = HTObj.HT{p}.first_node;
            while (HTObj.HT{p}.current ~= HTObj.HT{p}.last_node)
                % iter =iter + 1;
                % disp(iter)
                %date fron cell in the p linked list 
                ind_i = HTObj.HT{p}.current.cell_struct.i;
                ind_j = HTObj.HT{p}.current.cell_struct.j;
                l = HTObj.HT{p}.current.cell_struct.l;
                id = HTObj.HT{p}.current.cell_struct.id;

                h = (Dif_dom_x)/(HTObj.N*2^(l));
                k = (Dif_dom_y)/(HTObj.M*2^(l));

                x = HTObj.a1 + ind_i*h;
                % x = HTObj.a1 + ind_i*(HTObj.b1-HTObj.a1)/(HTObj.N*2^(l));
                y = HTObj.a2 + ind_j*k;
                
                t = 16/(l+1);
                text(x+h*0.3,y+k*0.5,num2str(id),'FontSize',t); 

                % disp(id)
                HTObj.HT{p}.current = HTObj.HT{p}.current.next;
            end

            if (HTObj.HT{p}.current == HTObj.HT{p}.last_node)
                % iter =iter + 1;
                %disp(HTObj.HT{p})
                %date fron last cell in the p linked list 
                ind_i = HTObj.HT{p}.current.cell_struct.i;
                ind_j = HTObj.HT{p}.current.cell_struct.j;
                l = HTObj.HT{p}.current.cell_struct.l;
                id = HTObj.HT{p}.current.cell_struct.id;
                % disp(HTObj.HT{p}.current.cell_struct)

                h = (Dif_dom_x)/(HTObj.N*2^(l));
                k = (Dif_dom_y)/(HTObj.M*2^(l));


                x = HTObj.a1 + ind_i*h;
                y = HTObj.a2 + ind_j*k;

                t = 16/(l+1);
                text(x+h*0.3,y+k*0.5,num2str(id),'FontSize',t); 

            end

        end
    end

    axis('off');  %Shows the axes
end
