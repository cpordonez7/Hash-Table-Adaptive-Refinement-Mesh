function [pde] = getPDE(option)
%Esta funcion retorna pde tipo estructura. Esta hecha para tres ejm

    switch option
        case 1      %Ejemplo 1. Tesis profe Cata
            pde.fFun = @(x,y) funF(x,y);
            pde.g1Fun = @(x,y) funG1(x,y);
            pde.g2Fun = @(x,y) funG2(x,y);
            pde.g3Fun = @(x,y) funG3(x,y);
            pde.g4Fun = @(x,y) funG4(x,y);   
            pde.g1_DerFun = @(x,y) funDerG1(x,y);
            pde.g2_DerFun = @(x,y) funDerG2(x,y);
            pde.g3_DerFun = @(x,y) funDerG3(x,y);
            pde.g4_DerFun = @(x,y) funDerG4(x,y);
            pde.solexFun = @(x,y) funphi(x,y); 

        case 2      %Ejemplo 2. Tarea PDE
            pde.fFun = @(x,y) funF2(x,y);
            pde.g1Fun = @(x,y) funG12(x,y);
            pde.g2Fun = @(x,y) funG22(x,y);
            pde.g3Fun = @(x,y) funG32(x,y);
            pde.g4Fun = @(x,y) funG42(x,y);    
            pde.solexFun = @(x,y) funphi2(x,y);

        case 3      %Ejemplo 3. Art
            pde.fFun = @(x,y) funF3(x,y);
            pde.g1Fun = @(x,y) funG13(x,y);
            pde.g2Fun = @(x,y) funG23(x,y);
            pde.g3Fun = @(x,y) funG33(x,y);
            pde.g4Fun = @(x,y) funG43(x,y);    
            pde.solexFun = @(x,y) funphi3(x,y);

        case 4    %Ejemplo 4. gota
            pde.fFun = @(x,y) funF4(x,y);
            pde.g1Fun = @(x,y) funG14(x,y);
            pde.g2Fun = @(x,y) funG24(x,y);
            pde.g3Fun = @(x,y) funG34(x,y);
            pde.g4Fun = @(x,y) funG44(x,y);    
            pde.solexFun = @(x,y) funphi4(x,y);

        case 5    %Ejemplo 5. Problema manufacturado mio
            pde.fFun = @(x,y) funF5(x,y);
            pde.g1Fun = @(x,y) funG15(x,y);
            pde.g2Fun = @(x,y) funG25(x,y);
            pde.g3Fun = @(x,y) funG35(x,y);
            pde.g4Fun = @(x,y) funG45(x,y);   
            pde.g1_DerFun = @(x,y) funDerG15(x,y);
            pde.g2_DerFun = @(x,y) funDerG25(x,y);
            pde.g3_DerFun = @(x,y) funDerG35(x,y);
            pde.g4_DerFun = @(x,y) funDerG45(x,y);
            pde.solexFun = @(x,y) funphi5(x,y);
        otherwise
            disp('Error, ingrese un valor válido.')  
    end
end

%Función F
function u = funF(x,y)
    u = -8.0*pi*pi*sin(2*pi*x)*sin(2*pi*y);
end

%Condiciones de contorno G1, G2, G3 Y G4
function u = funG1(x,y)
       u = sin(2*pi*x)*sin(2*pi*y); 
end

function u = funG2(x,y)
    u = sin(2*pi*x)*sin(2*pi*y);
end

function u = funG3(x,y)  
    u    = sin(2*pi*x)*sin(2*pi*y);
end
function u = funG4(x,y)
    u = sin(2*pi*x)*sin(2*pi*y);
end

%Derivadas
function u = funDerG1(x,y)
    u = cos(2*pi*y)*2*pi*sin(2*pi*x);
end

function u = funDerG2(x,y)
    u = cos(2*pi*x)*2*pi*sin(2*pi*y);
end

function u = funDerG3(x,y)
    u = cos(2*pi*x)*2*pi*sin(2*pi*y);
end

function u = funDerG4(x,y)
    u = cos(2*pi*y)*2*pi*sin(2*pi*x);
end

%Función Phi
function solex = funphi(x,y) 
    solex = sin(2*pi*x)*sin(2*pi*y);
    %y1 = x*exp(y);
end

%%%%% FUNCIONES SEGUNDO EJEMPLO %%%%%% 

%Función F
function u = funF2(x,y)
    u = -(exp(10*(x^2+y^2))/50.0)*(x*(-x+1)*(-200*y^4+200*y^3-50*y^2+30*y-1)+y*(-y+1)*(-200*x^4+200*x^3-50*x^2+30*x-1));
end

%Condiciones de contorno G12, G22, G32 Y G42
function u = funG12(x,y)
       u = (1.0/100.0)*x*(1-x)*y*(1-y)*exp(10*(x^2+y^2)); 
end

function u = funG22(x,y)
    u = (1.0/100.0)*x*(1-x)*y*(1-y)*exp(10*(x^2+y^2));
end

function u = funG32(x,y)  
    u    = (1.0/100.0)*x*(1-x)*y*(1-y)*exp(10*(x^2+y^2));
end
function u = funG42(x,y)
    u = (1.0/100.0)*x*(1-x)*y*(1-y)*exp(10*(x^2+y^2));
end

%Función Phi
function solex = funphi2(x,y) 
    solex = (1.0/100.0)*x*(1-x)*y*(1-y)*exp(10*(x^2+y^2));
end



% EJEMPLO 3

%Función F
function u = funF3(x,y)
    u = x*exp(y);
end

%Condiciones de contorno G1, G2, G3 Y G4 (B,L;A,R)
function u = funG13(x,y)
       u = x; 
end

function u = funG23(x,y)
    u = 0.0;
end

function u = funG33(x,y)  
    u    = exp(1)*x;
end
function u = funG43(x,y)
    u = 2*exp(y);
end

%Función Phi
function solex = funphi3(x,y) 
    solex = x*exp(y);
    %y1 = x*exp(y);
end

% EJEMPLO 4

%Función F
function u = funF4(x,y)
    xc = 1.0/13.0;
    yc = 1.0/6.0;
    d = (x-xc)^2+(y-yc)^2;
    d1 = sqrt(d);
    u2 = ((3*45)*(1.0/6.0-y)*(1-(tanh(45*d1-22.5))^2)*(y-1.0/6.0))/(d^(3.0/2.0));
    u1 = ((3*45)*(1.0/13.0-x)*(1-(tanh(45*d1-22.5))^2)*(x-1.0/13.0))/(d^(3.0/2.0));
    u3 = ((3*4050)*(1-(tanh(45*d1-22.5))^2)*((x-1.0/13.0)^2)*tanh(45*d1-22.5))/(d);
    u4 = ((3*4050)*(1-(tanh(45*d1-22.5))^2)*((y-1.0/6.0)^2)*tanh(45*d1-22.5))/(d);
    u5 = ((3*90)*(1-(tanh(45*d1-22.5))^2))/(d1);


    % u2 = ((75)*(1.0/6.0-y)*(1-(tanh(75*d1-37.5))^2)*(y-1.0/6.0))/(d^(3.0/2.0));
    % u1 = ((75)*(1.0/13.0-x)*(1-(tanh(75*d1-37.5))^2)*(x-1.0/13.0))/(d^(3.0/2.0));
    % u3 = ((11250)*(1-(tanh(75*d1-37.5))^2)*((x-1.0/13.0)^2)*tanh(75*d1-37.5))/(d);
    % u4 = ((11250)*(1-(tanh(75*d1-37.5))^2)*((y-1.0/6.0)^2)*tanh(75*d1-37.5))/(d);
    % u5 = ((150)*(1-(tanh(75*d1-37.5))^2))/(d1);
    
    u = -u1 - u2 + u3 + u4 - u5;
end

%Condiciones de contorno G1, G2, G3 Y G4 (B,L;A,R)
function u = funG14(x,y)
    xc = 1.0/13.0;
    yc = 1.0/6.0;
    d = sqrt((x-xc)^2+(y-yc)^2);

    u = 3*tanh(45*(0.5-d));
    % u = tanh(75*(0.5-d));
end

function u = funG24(x,y)
    xc = 1.0/13.0;
    yc = 1.0/6.0;
    d = sqrt((x-xc)^2+(y-yc)^2);

    u = 3*tanh(45*(0.5-d));
    % u = tanh(75*(0.5-d));
end

function u = funG34(x,y)  
    xc = 1.0/13.0;
    yc = 1.0/6.0;
    d = sqrt((x-xc)^2+(y-yc)^2);

    u = 3*tanh(45*(0.5-d));
    % u = tanh(75*(0.5-d));
end
function u = funG44(x,y)
    xc = 1.0/13.0;
    yc = 1.0/6.0;
    d = sqrt((x-xc)^2+(y-yc)^2);

    u = 3*tanh(45*(0.5-d));
    % u = tanh(75*(0.5-d));
end

%Función Phi
function solex = funphi4(x,y) 
    xc = 1.0/13.0;
    yc = 1.0/6.0;
    d = sqrt((x-xc)^2+(y-yc)^2);

    solex = 3*tanh(45*(0.5-d));
    % solex = tanh(75*(0.5-d));
    
    %y1 = x*exp(y);
end



%%%%%%%%%%%%%%%%%%%%%%%%5
% EJEMPLO 5

%Función F
function u = funF5(x,y)
    
    u = -4*y^2*exp(0.1*x)*sin(y^2)+2*(2*x+2)*cos(x+y)-2*(x^2+2*x)*sin(x+y)+0.01*exp(0.1*x)*sin(y^2)+2*exp(0.1*x)*cos(y^2)+2*sin(x+y);
   
end

%Condiciones de contorno G1, G2, G3 Y G4 (B,L;A,R)
function u = funG15(x,y)

    u = exp(0.1*x)*sin(y^2)+(x^2+2*x)*sin(x+y);
end

function u = funG25(x,y)

    u = exp(0.1*x)*sin(y^2)+(x^2+2*x)*sin(x+y);
end

function u = funG35(x,y)

    u = exp(0.1*x)*sin(y^2)+(x^2+2*x)*sin(x+y);
end

function u = funG45(x,y)

    u = exp(0.1*x)*sin(y^2)+(x^2+2*x)*sin(x+y);
end

%Derivadas
function u = funDerG15(x,y)
    u = 2*y*exp(0.1*x)*cos(y*y)+(x*x+2*x)*cos(x+y);
end

function u = funDerG25(x,y)
    u = (2*x+2)*sin(x+y)+(x*x+2*x)*cos(x+y)+0.1*exp(0.1*x)*sin(y*y);
end

function u = funDerG35(x,y)
    u = (2*x+2)*sin(x+y)+(x*x+2*x)*cos(x+y)+0.1*exp(0.1*x)*sin(y*y);
end

function u = funDerG45(x,y)
    u = 2*y*exp(0.1*x)*cos(y*y)+(x*x+2*x)*cos(x+y);
end


%Función Phi
function solex = funphi5(x,y) 

    solex = exp(0.1*x)*sin(y^2)+(x^2+2*x)*sin(x+y);

end
