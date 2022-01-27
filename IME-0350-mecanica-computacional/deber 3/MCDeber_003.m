% ********** Deber 03 de Mecánica Computacional - Módulo 1 *********
clear
clc

%Constantes del problema
k = 240; %[W/mK]
w = 0.020; %[m]
W = 0.040; %[m]
T_s = 50; %[°C] Temp de superficie
T_r = 20; %[°C] Temp del refrigerante
dx = 0.005; %[m]
dy = 0.005; %[m]

%Valor ingresado por usuario
h = input('ingrese un vector con valores de h: '); 

i = 1;
T = [];
Q =[];
while i <= length(h)
    
    %Se arma un sistema "Ax = b". El vector "x" contendrá las temperaturas.
    
    A = [  -3-h(i)*dy/k 1 0 0 2 0 0; 
           1 -4-2*dy*h(i)/k 1 0 0 2 0;
           0 1 -2-h(i)*dx/k 0 0 0 1;
           0 0 0 -2 1 0 0;
           1 0 0 1 -4 1 0;
           0 1 0 0 1 -4 1 ;
           0 0 1 0 0 2 -4 ;
        ];
    
    b = [ -T_r*h(i)*dy/k
           -2*T_r*dy*h(i)/k
           -h(i)*T_r*dx/k
           -T_s
           -T_s 
           -T_s
           -T_s
         ];
     
    %Anexar los valores de T hallados para cada valor de h
    T = [T (A\b);];
    q = -8*h(i)*dy *((1/2)*(T_r - T(1,i)) + (T_r - T(2,i)) + (1/2)*(T_r - T(3,i)));
    Q = [Q q;];
    
    %Actualizar el contador
    i = i+1;
end 

R = table(h',T', Q');




