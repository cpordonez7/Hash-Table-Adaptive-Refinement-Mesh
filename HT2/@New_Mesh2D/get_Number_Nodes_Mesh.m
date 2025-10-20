function [Number_Nodes_Mesh] = get_Number_Nodes_Mesh(HTObj)
    %This function return the number of nodes that the nodes array 
    %Calcula el numero de nodos en el arreglo que ingresa, nodos en la malla. 
    
    Number_Nodes_Mesh = length(HTObj.Mesh_nodes);

end