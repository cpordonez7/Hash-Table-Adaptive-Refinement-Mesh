function [] = set_initialize_Cell_counter_level(HTObj)
    %This function_initialize the vector counter level,It is used when making a remapper 

    HTObj.Cell_counter_level = zeros(max_level+1,1);

end