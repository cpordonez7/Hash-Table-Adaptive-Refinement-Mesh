function [List_Neighbors_quadr_interp, flag_direction, flag_aux] = getNeighbors_quadratic_interpolation(HTObj,Id,F1,F3,Dir,m)
    %This function get the neightbors for makes the cuadratic interpolation
    %Id is a cell center, F1 and F3 are neighbors cells in the m direction
    %in orden (L_R) or (B-A).
    %Dir B,L, R, A m 1, 2, 3, 4 Dir = m

    flag_aux = 0;
    F2 = HTObj.getNeighbors_direction_id(F1,Dir);
    F4 = HTObj.getNeighbors_direction_id(F3,Dir);

    if length(F1) ~= 1 && length(F2) ~= 1%en vez de F1 va F4
        fprintf('Error, cannot apply quadratic interpolation for this mesh');
        return;
    end


    if m == 1 || m == 4 %(Below-Above)
        %Dif_level countain (Id, dif_level in this direction)
        Dif_level_1 = HTObj.getNeighbors_left_id_dif_level(F1); %left
        Dif_level_2 = HTObj.getNeighbors_right_id_dif_level(F3); %Right

        if m == 1
            Dif_level_3 = HTObj.getNeighbors_below_id_dif_level(F2);
        else 
            Dif_level_3 = HTObj.getNeighbors_above_id_dif_level(F2);
        end
  

        if (Dif_level_1(2,1) == 0 && Dif_level_1(1,1) ~= 0)
        % if Dif_level_1(2,1) == 0 && Dif_level_1(2,2) == 0
            flag_direction = 1;
            F5 = Dif_level_1(1,1);
            list_vec = HTObj.getNeighbors_left_id_dif_level(F2);
            F6 = list_vec(1,1);
            F7 = Id;

        elseif (Dif_level_2(2,1) == 0 && Dif_level_2(1,1) ~= 0)
        % elseif Dif_level_2(2,1) == 0 && Dif_level_2(2,2) == 0
            flag_direction = 2;
            F5 = Dif_level_2(1,1);
            list_vec = HTObj.getNeighbors_right_id_dif_level(F4);
            F6 = list_vec(1,1);
            F7 = Id;
        elseif (Dif_level_3(2,1) == 0 && Dif_level_3(1,1) ~= 0)
        % elseif Dif_level_3(2,1) == 0 && Dif_level_3(2,2) == 0
            flag_direction = 3;
            F5 = Dif_level_3(1,1);
            if m ==1
                list_vec = HTObj.getNeighbors_below_id_dif_level(F4);
                F6 = list_vec(1,1);
            else
                F6 = F5;
                list_vec = HTObj.getNeighbors_above_id_dif_level(F4);
                F5 = list_vec(1,1);
            end

            if (Dif_level_1(2,1) < 0 && Dif_level_1(1,1) ~= 0)
                flag_aux = 1;
                F7 = Dif_level_1(1,1);
            elseif (Dif_level_2(2,1) < 0 && Dif_level_2(1,1) ~= 0)
                flag_aux = 2;
                F7 = Dif_level_2(1,1);
            else
                fprintf('Error, cannot apply quadratic interpolation');
                return;
            end
        else

            fprintf('Error, cannot apply quadratic interpolation');
            return;
        end

        if m == 1
            List_Neighbors_quadr_interp = [F1 F2 F3 F4 F5 F6 F7];
        else
            List_Neighbors_quadr_interp = [F3 F4 F1 F2 F5 F6 F7];
        end

    else
        Dif_level_1 = HTObj.getNeighbors_below_id_dif_level(F1);
        Dif_level_2 = HTObj.getNeighbors_above_id_dif_level(F3);

        if m == 2
            Dif_level_3 = HTObj.getNeighbors_left_id_dif_level(F2);
        else 
            Dif_level_3 = HTObj.getNeighbors_right_id_dif_level(F2);
        end

        if (Dif_level_1(2,1) == 0 && Dif_level_1(1,1) ~= 0)
        % if Dif_level_1(2,1) == 0 && Dif_level_1(2,2) == 0
            flag_direction = 1;
            F5 = Dif_level_1(1,1);
            list_vec = HTObj.getNeighbors_below_id_dif_level(F2);
            F6 = list_vec(1,1);
            F7 = Id;

        elseif (Dif_level_2(2,1) == 0 && Dif_level_2(1,1) ~= 0)
        % elseif Dif_level_2(2,1) == 0 && Dif_level_2(2,2) == 0
            flag_direction = 2;
            F5 = Dif_level_2(1,1);
            list_vec = HTObj.getNeighbors_above_id_dif_level(F4);
            F6 = list_vec(1,1);
            F7 = Id;

        elseif (Dif_level_3(2,1) == 0 && Dif_level_3(1,1) ~= 0)
        % elseif Dif_level_3(2,1) == 0 && Dif_level_3(2,2) == 0
            flag_direction = 3;
            F5 = Dif_level_3(1,1);
            if m ==2
                F6 = F5;
                list_vec = HTObj.getNeighbors_left_id_dif_level(F4);
                F5 = list_vec(1,1);
            else
                list_vec = HTObj.getNeighbors_right_id_dif_level(F4);
                F6 = list_vec(1,1);
            end

            if (Dif_level_1(2,1) < 0 && Dif_level_1(1,1) ~= 0)
                flag_aux = 1;
                F7 = Dif_level_1(1,1);
            elseif (Dif_level_2(2,1) < 0 && Dif_level_2(1,1) ~= 0)
                flag_aux = 2;
                F7 = Dif_level_2(1,1);
            else
                fprintf('Error, cannot apply quadratic interpolation');
                return;
            end
        else
            fprintf('Error, cannot apply quadratic interpolation');
            return;
        end

        if m == 2
            List_Neighbors_quadr_interp = [F3 F4 F1 F2 F5 F6 F7];
        else
            List_Neighbors_quadr_interp = [F1 F2 F3 F4 F5 F6 F7];
        end
    end

end