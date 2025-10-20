function [cell_strc] = create_cell(i,j,l)
%This function returns a structure type cell.
%
%Input: 
% (i,j)-> indices of coordinates.
% l-> cell level.
%
%Output : Structure with (i,j,l,key) values. 
%All cells are initialized with null key.
 
    cell_strc.i = i;
    cell_strc.j = j;
    cell_strc.l = l;
    cell_strc.key = 'Null';
    cell_strc.id = 0;
end