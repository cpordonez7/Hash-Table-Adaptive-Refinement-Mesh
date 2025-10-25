%Test New Mesh2D

% Dominio
% malla test 1 [0x1]^2
% malla test 2 [-0.5,0.5]^2

%Structur for domain
grid.N = 32;
grid.M = 32;
grid.a1 = 0.0;
grid.a2 = 0.0;
grid.b1 = 1.0;
grid.b2 = 1.0;

% grid.N = 4;
% grid.M = 4;
% grid.a1 = 0.0;
% grid.a2 = 0.0;
% grid.b1 = 1.0;
% grid.b2 = 1.0;

%created a object
% Create_uniform_mesh(grid.N,grid.M,grid.a1,grid.b1,grid.a2,grid.b2);
% ObjHT = New_Mesh2D('Base_uniform_mesh.txt',2,grid);
ObjHT = New_Mesh2D('File_current_mesh_creation.txt',3,grid);
% ObjHT = New_Mesh2D('malla_test_1_lagrange_8x8_3_level.txt',3,grid);
% ObjHT = New_Mesh2D('malla_test_1_lagrange_8x8_3_level_correcta_anidada.txt',3,grid);
% ObjHT = New_Mesh2D('Malla_test2_16x16_4_level.txt',4,grid);
% ObjHT = New_Mesh2D('Malla_test2_8x8_4_level.txt',4,grid);
% ObjHT = New_Mesh2D('Malla_test_particulas_8x8_3_level.txt',3,grid);
% ObjHT = New_Mesh2D('malla_test4_16x16_4levels.txt',4,grid);
% ObjHT = New_Mesh2D('File_current_mesh.txt',3,grid);
% ObjHT = New_Mesh2D('Malla_test3_32x32_3_level_new.txt',3,grid);
 % ObjHT = New_Mesh2D('Mesh_base_uniforme_6x4.txt',1,grid);
% ObjHT = New_Mesh2D('Malla_refinada_L_centro_6x4.txt',2,grid);

% ObjHT.set_id_celulas_for_level
% ObjHT.set_enumeration_vector
% cell_1 = ObjHT.getCell_id(1);
% cell_6 = ObjHT.getCell_id(6);
% cell_11 = ObjHT.getCell_id(11);
% ObjHT.RefineCell(cell_1);
% ObjHT.RefineCell(cell_6);
% ObjHT.RefineCell(cell_11);
% 
% ObjHT.set_enumeration_vector;
% cell_6 = ObjHT.getCell_id(6);
% ObjHT.ThickenCell(cell_6);
% % ObjHT.set_enumeration_vector;
% cell_11 = ObjHT.getCell_id(11);
% ObjHT.ThickenCell(cell_11);
% % ObjHT.set_enumeration_vector;
% cell_1 = ObjHT.getCell_id(1);
% ObjHT.ThickenCell(cell_1);
% ObjHT.set_id_celulas_for_level;
% ObjHT.set_enumeration_vector;

% cell_15 = ObjHT.getCell_id(15);
% ObjHT.ThickenCell(cell_15);
% ObjHT.set_enumeration_vector;
% cell_19 = ObjHT.getCell_id(19);
% ObjHT.ThickenCell(cell_19);
% cell_6 = ObjHT.getCell_id(6);
% cell_11 = ObjHT.getCell_id(11);
% ObjHT.RefineCell(cell_1);
% ObjHT.set_enumeration_vector;
% % cell_16 = ObjHT.getCell_id(16);
% % cell_9 = ObjHT.getCell_id(9);
% % cell_11 = ObjHT.getCell_id(11);
% % 
% % ObjHT.Refine_region(cell_9,cell_16);
% % ObjHT.RefineCell(cell_11);
% % 
% % ObjHT.set_enumeration_vector;
% % ObjHT.set_id_celulas_for_level
% % % % 
% % % % cell_1 = ObjHT.getCell_id(1);
% % % % cell_9 = ObjHT.getCell_id(9);
% % % % cell_3 = ObjHT.getCell_id(3);
% % % % cell_4 = ObjHT.getCell_id(4);
% % % % cell_11 = ObjHT.getCell_id(11);
% % % % cell_22 = ObjHT.getCell_id(22);
% % % % cell_25 = ObjHT.getCell_id(25);
% % % % cell_26 = ObjHT.getCell_id(26);
% % % % cell_33 = ObjHT.getCell_id(33);
% % % % cell_38 = ObjHT.getCell_id(38);
% % % % cell_39 = ObjHT.getCell_id(39);
% % % % cell_49 = ObjHT.getCell_id(49);
% % % % cell_59 = ObjHT.getCell_id(59);
% % % % cell_46 = ObjHT.getCell_id(46);
% % % % cell_62 = ObjHT.getCell_id(62);
% % % % cell_63 = ObjHT.getCell_id(63);
% % % % 
% % % % 
% % % % ObjHT.Refine_region(cell_1,cell_9);
% % % % ObjHT.Refine_region(cell_3,cell_11);
% % % % ObjHT.Refine_region(cell_4,cell_22);
% % % % ObjHT.Refine_region(cell_25,cell_26);
% % % % ObjHT.Refine_region(cell_38,cell_39);
% % % % ObjHT.Refine_region(cell_49,cell_59);
% % % % ObjHT.Refine_region(cell_62,cell_63);
% % % % ObjHT.RefineCell(cell_33);
% % % % ObjHT.RefineCell(cell_46);
% % % % % % % 
% % % % % % % 
% % % % ObjHT.set_enumeration_vector;
% % % % ObjHT.set_id_celulas_for_level
% % % % ObjHT.set_enumeration_vector
% % % % cell_126 = ObjHT.getCell_id(126);
% % % % cell_137 = ObjHT.getCell_id(137);
% % % % cell_52 = ObjHT.getCell_id(52);
% % % % cell_65 = ObjHT.getCell_id(65);
% % % % cell_74 = ObjHT.getCell_id(74);
% % % % cell_75 = ObjHT.getCell_id(75);
% % % % % 
% % % % ObjHT.Refine_region(cell_126,cell_137);
% % % % ObjHT.Refine_region(cell_52,cell_65);
% % % % ObjHT.RefineCell(cell_74);
% % % % ObjHT.RefineCell(cell_75);


% % ObjHT.RefineCell(cell_82);
% ObjHT.Refine_region(cell_90,cell_91);
% 
% ObjHT.ThickenCell(cell_29);
% ObjHT.ThickenCell(cell_125);
% ObjHT.ThickenCell(cell_70);
% ObjHT.ThickenCell(cell_168);

% ObjHT.set_enumeration_vector;
% % cell_109 = ObjHT.getCell_id(109);
% % cell_114 = ObjHT.getCell_id(114);
% % cell_132 = ObjHT.getCell_id(132);
% % cell_143 = ObjHT.getCell_id(143);
% % cell_47 = ObjHT.getCell_id(47);
% % cell_49 = ObjHT.getCell_id(49);
% % cell_59 = ObjHT.getCell_id(59);
% % cell_62 = ObjHT.getCell_id(62);
% % cell_63 = ObjHT.getCell_id(63);
% % % % % 
% % % % % 
% % % % % % 
% % ObjHT.Refine_region(cell_1,cell_9);
% % 
% % ObjHT.Refine_region(cell_50,cell_63);
% % ObjHT.Refine_region(cell_72,cell_73);
% % ObjHT.Refine_region(cell_109,cell_114);
% % ObjHT.Refine_region(cell_132,cell_143);
% ObjHT.Refine_region(cell_62,cell_63);
% % ObjHT.RefineCell(cell_101);
% % ObjHT.RefineCell(cell_22);
% % ObjHT.RefineCell(cell_82);


% % % 
% % % 
% % % ObjHT.set_enumeration_vector;
% % % 
% % % cell_106 = ObjHT.getCell_id(106);
% % % cell_142 = ObjHT.getCell_id(142);
% % % cell_146 = ObjHT.getCell_id(146);
% % % cell_148 = ObjHT.getCell_id(148);
% % % cell_150 = ObjHT.getCell_id(150);
% % % cell_126 = ObjHT.getCell_id(126);
% % % cell_128 = ObjHT.getCell_id(128);
% % % 
% % % ObjHT.ThickenCell(cell_106);
% % % ObjHT.ThickenCell(cell_142);
% % % ObjHT.ThickenCell(cell_146);
% % % ObjHT.ThickenCell(cell_148);
% % % ObjHT.ThickenCell(cell_150);
% % % ObjHT.ThickenCell(cell_126);
% % % ObjHT.ThickenCell(cell_128);

% % ObjHT.set_enumeration_vector;
% % ObjHT.set_id_celulas_for_level
% % 
% % cell_2 = ObjHT.getCell_id(2);
% % cell_11 = ObjHT.getCell_id(11);
% % cell_51 = ObjHT.getCell_id(51);
% % cell_62 = ObjHT.getCell_id(62);
% % cell_69 = ObjHT.getCell_id(69);
% % cell_70 = ObjHT.getCell_id(70);
% % cell_123 = ObjHT.getCell_id(123);
% % cell_134 = ObjHT.getCell_id(134);
% % cell_100 = ObjHT.getCell_id(100);
% % cell_105 = ObjHT.getCell_id(105);

% % % 
% % % 
% % % % % 
% ObjHT.Refine_region(cell_2,cell_11);
% 
% ObjHT.Refine_region(cell_51,cell_62);
% ObjHT.Refine_region(cell_69,cell_70);
% ObjHT.Refine_region(cell_123,cell_134);
% % ObjHT.Refine_region(cell_100,cell_105);
% 
% ObjHT.set_enumeration_vector;
% ObjHT.set_id_celulas_for_level


% % % 
% % % 
% % % 
% % % % 
% cell_ = create_cell(1,1,0);
% cell_.id = 6;
% 
% cell_2 = create_cell(2,2,0);
% cell_2.id = 16;
% Mesh2D_HTObj.Refine_region(cell_,cell_2);
% % % % 
% % % % % 
% % % % cell_4 = create_cell(3,3,1);
% % % % cell_4.id = 18;
% % % % 
% % % % cell_2 = create_cell(2,2,1);
% % % % cell_2.id = 6;
% % % % 
% % % % Mesh2D_HTObj.Refine_region(cell_2,cell_4);
% % % 
% % % % 
% % % % cell_5 = create_cell(0,3,0);
% % % % cell_5.id = 13;
% % % % 
% % % % cell_6 = create_cell(0,6,1);
% % % % cell_6.id = 13;
% % % % 
% % % % Mesh2D_HTObj.RefineCell(cell_);
% % % % Mesh2D_HTObj.RefineCell(cell_2);
% % % % Mesh2D_HTObj.RefineCell(cell_3);
% % % % Mesh2D_HTObj.RefineCell(cell_4);
% % % % Mesh2D_HTObj.RefineCell(cell_5);
% % % % Mesh2D_HTObj.RefineCell(cell_6);
% % % 
% % % 
