function [List_Neighbors_P] = getNeighbors_Punto(HTObj, P, lc)
    %Return the id of the neighbors cells to the referent point P(x,y)
    % disp(P)
    %Difference domial extrem
    Dif_dom_x = HTObj.b1 - HTObj.a1;
    Dif_dom_y = HTObj.b2 - HTObj.a2;

    x = P(1,1);         %x coordinate of the node         
    y = P(2,1);         %y coordinate of the node
    l = lc;

    %Se busca en qué nivel está ese punto
    flag_1 = 0;
    iter = 0;
    while flag_1 == 0 && l<=HTObj.Max_level-1
        l = l + iter;
        % disp('level')
        % disp(l)
        HT = HTObj.get_HT();
        ind_i = HTObj.N*2^(l)*(x-HTObj.a1)/(Dif_dom_x);
        ind_j = HTObj.M*2^(l)*(y-HTObj.a2)/(Dif_dom_y);

        ind_i = floor(ind_i);
        ind_j = floor(ind_j);
        % disp(ind_i)
        % disp(ind_j)

        %Key calulation
        cell_punto = create_cell(ind_i,ind_j,l);
        key = HTObj.key_Calculation(cell_punto);
        p = HTObj.hast_function(key);

        size_HT_p = size(HT{p+1},2);
        if (size_HT_p ~=0)
            %True or False
            flag = HT{p+1}.searchNode(key);
            % disp(flag)

            if strcmp(flag, 'True')
                cell_punto = HT{p+1}.get_cell_Nodo(key);
                lc = cell_punto.l;
                id_cel  = cell_punto.id;
                if l == lc
                    flag_1 = 1;
                    l = lc;
                end
            end
        end
        iter = iter + 1;
    end 
    
    % disp(id_cel)

    List_Neighbors_below_P = HTObj.getNeighbors_below_id(id_cel);
    List_Neighbors_left_P = HTObj.getNeighbors_left_id(id_cel);

    Neighbors_below = List_Neighbors_below_P(1,1);
    Neig_Neig_below = HTObj.getNeighbors_left_id(Neighbors_below);

    num_neig_NNb = length(Neig_Neig_below);

    Neighbors_diagonal = Neig_Neig_below(1,num_neig_NNb);

    List_Neighbors_P = [id_cel Neighbors_below List_Neighbors_left_P(1,1) Neighbors_diagonal];
    
end