classdef create_node_cell < handle
    %   This class creates a object type node, this object has atributos
    %   key(value), previous (previous node) and next (pointer, next node)

%-----------------------------------------------------------------
%               PROPERTIES SECTION
%-----------------------------------------------------------------
    properties (Access = public)
        cell_struct %cell 
        next        %pointer
        previous    %previous node
        value       %value node - key
    end

%-----------------------------------------------------------------
%               METHODS SECTION
%-----------------------------------------------------------------
    methods
        function [obj] = create_node_cell(cell_struct)
            %UNTITLED2 Construct an instance of this class
            % Initialize the node properties.
            obj.cell_struct = cell_struct;
            obj.value = cell_struct.key;
            obj.previous = 'Null';
            obj.next = 'Null';
        end



%-----------------------------------------------------------------
%               SECCION INTERFACES
%-----------------------------------------------------------------
    %En esta seccion se agregan los métodos cuandolos las propiedades de
    %los objetos son privadas
    end

end