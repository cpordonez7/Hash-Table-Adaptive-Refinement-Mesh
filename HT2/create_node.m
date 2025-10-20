classdef create_node < handle
    %   This class creates a object type node, this object has atributos Id
    %   (value), previous (previous node) and next (pointer, next node)

%-----------------------------------------------------------------
%               PROPERTIES SECTION
%-----------------------------------------------------------------
    properties (Access = public)
        next        %pointer
        previous    %previous node
        value       %value node
    end

%-----------------------------------------------------------------
%               METHODS SECTION
%-----------------------------------------------------------------
    methods
        function [obj] = create_node(value)
            %UNTITLED2 Construct an instance of this class
            % Initialize the node properties.
            obj.value = value;
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