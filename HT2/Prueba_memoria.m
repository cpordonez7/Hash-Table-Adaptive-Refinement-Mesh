function Prueba_memoria
ens = cell(2,3);
% ens{1} = {struct('x',7,'y',7), struct('x',8,'y',8),struct('x',9,'y',9)};
% ens{2} = {struct('x',0,'y',0), struct('x',1,'y',1),struct('x',2,'y',2)};
% a = ens{1};
% b = ens{2};
% a1= ens{1}(1);
% a2 = ens{1}(2);
% b1 = ens{2}(1);


ens = cell(2,1);
cel1 = create_cell(0,0,0);
cel1.key = 1;
cel2 = create_cell(1,1,1);
cel2.key = 2;
cel3 = create_cell(2,2,2);
cel3.key = 3;
Node1 = create_node_cell(cel1); 
Node2 = create_node_cell(cel2);

ens{1} = Linked_list_cell(Node1);
ens{2}= Linked_list_cell(Node2);
ens{2}.addNode_cell(cel3);
a = ens{1};
b = ens{2};
% print_list(a)
whos 