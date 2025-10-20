function [key] = key_Calculation_i_j_l(HTObj,ind_i,ind_j,l)
    %key_Calculation: This function calculate the cell key.
    %
    %Input: 
    %HTObj -> Object Hast Table
    %(i,j,l) -> index and level of the cell structure 
 

    l_max = HTObj.Max_level;            %Max level
    m = HTObj.N*2^(l_max);              %number of horizontal cell in the max level
    key = (ind_i+ind_j*m)*2^(l_max-l);          %key

    % key = (ind_i+ 1) + (HTObj.N*2^(l)*ind_j);  
    
end