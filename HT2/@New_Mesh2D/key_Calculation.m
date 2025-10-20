function [key] = key_Calculation(HTObj,cell_str)
    %key_Calculation: This function calculate the cell key.
    %
    %Input: 
    %HTObj -> Object Hast Table
    %cell_str -> cell structure 
 
    %Indices (i,j) and level.
    ind_i = cell_str.i;
    ind_j = cell_str.j;
    l = cell_str.l;

    l_max = HTObj.Max_level;            %Max level
    m = HTObj.N*2^(l_max);              %number of horizontal cell in the max level
    key = (ind_i+ind_j*m)*2^(l_max-l);          %key

    % key = (ind_i+ 1) + (HTObj.N*2^(l)*ind_j);  
    
end