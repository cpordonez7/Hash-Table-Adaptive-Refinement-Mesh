classdef Node < handle
    %   This class creates an object with attributes value and pointer

%-----------------------------------------------------------------
%               PROPERTIES SECTION
%-----------------------------------------------------------------
    properties (Access = public)
        key         %node value
        next        %pointer  
        
    end

%-----------------------------------------------------------------
%               METHODS SECTION
%-----------------------------------------------------------------
    methods
        function [NodeObj] = Node(key)
            %Constructor: Construct an instance of this class. 
            % Initializes the key value with the input parameter and by
            % default it makes the pointer equal to null.

            NodeObj.key = key;
            NodeObj.next = 'null';
        end



%-----------------------------------------------------------------
%               INTERFACES SECTION
%-----------------------------------------------------------------
    %In this section the methods are added when the properties of objects
    %are private
    end

end