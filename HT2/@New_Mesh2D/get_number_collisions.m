function [collisions_vector] = get_number_collisions(HTObj)
    %This function counts the number of collisiones that there are in each
    %lisked list of the HT.

    collisions_vector = zeros(HTObj.size_HT,1);

    for p=1:HTObj.size_HT

        if size(HTObj.HT{p},2) ~= 0 && HTObj.HT{p}.head.counter ~= 0
            collisions_vector(p) = HTObj.HT{p}.size_list;

        else
            collisions_vector(p) = 0;            

        end
    end
end
