%Test HT
%Examples of mesh operations and other methods in the HT

%%Creating a uniform mesh
%
%1. Define the domain
% Domain [a1,b1]x[a2,b2] and partitions in x-axis = N and y-axis = M
grid.N = 4;
grid.M = 4;
grid.a1 = 0.0;
grid.a2 = 0.0;
grid.b1 = 1.0;
grid.b2 = 1.0;

% Creating the mesh file
% Input: Domain information and partitions along the axes
% Output: Mesh data file named 'Base_uniform_mesh.txt'
% This file has a structure where the first data is the cell number in the
% mesh, followed by (x, y, l) coordinates and the cell level
Create_uniform_mesh(grid.N,grid.M,grid.a1,grid.b1,grid.a2,grid.b2);

%Creatind the Object
%Input: Data file, maximum number of levels and grid
%Output: An object with all the grid information
ObjHT = New_Mesh2D('Base_uniform_mesh.txt',2,grid);

%Validation of methods

%Graph: Use the graph method of the object
ObjHT.graph_Mesh2D;

%Enumeration
%Initially, the mesh has a lexicographic enumeration at level zero
ObjHT.graph_enumeration;

%Refine
% 1. Get the structures of the cells to refine with id parameter.
Cell6 = ObjHT.getCell_id(6);
Cell7 = ObjHT.getCell_id(7);

% 2. Refine the cells, calling the RefineCell method
ObjHT.RefineCell(Cell6);
ObjHT.RefineCell(Cell7);

figure
%Graph
ObjHT.graph_Mesh2D;
ObjHT.graph_enumeration;

%NOTE: When making a change to the Refine or Thicken mesh, it is necessary
%to update the enumeration vector, It allows the handling of the mesh
ObjHT.set_enumeration_vector;

%Change the enumeration: Enumeration by level
ObjHT.set_id_celulas_for_level;

figure
%Graph
ObjHT.graph_Mesh2D;
ObjHT.graph_enumeration;

%Thicken
% 1. Get the structures of the cells to refine with id parameter
Cell19 = ObjHT.getCell_id(19);

% 2. Thicken the cell, calling the ThickenCell method
ObjHT.ThickenCell(Cell19);
ObjHT.set_enumeration_vector;

figure
%Graph
ObjHT.graph_Mesh2D;
ObjHT.graph_enumeration;

%Generate file of the current mesh
%This method generartes a file with the data of the current mesh with the
%name 'File_current_mesh.txt'

% ObjHT.Generate_current_mesh_file;




