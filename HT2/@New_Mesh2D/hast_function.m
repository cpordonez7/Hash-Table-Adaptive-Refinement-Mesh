function [p] = hast_function(HTObj,key)
    %hast_function: This function calculates the position of the cell in
    %the HT.

    %Input: 
    %HTObj -> Object Hast Table
    %key -> key cell

    % p = mod(key,HTObj.size_HT);   %Modulo (Met. division)

    ra = (sqrt(5)-1)/2.0;
    kr = key*ra;
    p =  floor(HTObj.size_HT*(kr-floor(kr)));     %Met. Fibonacci
    
end