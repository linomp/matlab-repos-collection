clear all
clc

%**** CONSTANTES DEL PROBLEMA ******
L = 0.075; %[m]
K = 75; % [W/m*K]
q = 1.5 * 10 ^6; %[W/m^3]
h_1 = 1000; %[W/m^2*K]
h_2 = 700; %[W/m^2*K] 
T_fluido_1 = 30; %[°C]
T_fluido_2 = 50; %[°C]

%Cálculo de condiciones de borde
Ts_1 = (q*L/h_1) + T_fluido_1; %[°C]
Ts_2 = (q*L/h_2) + T_fluido_2; %[°C]
%***********************************


%******** Parte 1: RESOLUCIÓN POR DIFERENCIAS FINITAS ******
dx = 0.015/1
nx = (2*L/dx)+1
A = zeros(nx,nx); 
B = zeros(nx,1);

A(1,:) = [(-h_1-K/dx) K/dx zeros(1,nx-2)];
A(nx,:) = [zeros(1,nx-2) K/dx (-h_2-K/dx)];

B(1,1) = (-q*dx)-h_1*T_fluido_1;
B(nx,1) = (-q*dx)-h_2*T_fluido_2;

fila_A = [];

 for i=2:nx-1
     fila_A = [zeros(1,i-2) 1 -2 1 zeros(1,nx-1-i)];
     A(i,:) = fila_A;
     B(i,:) = [(-q*dx^2/K)];
 end 

T_1 = A\B;
x_1 = (-L:dx:L);

subplot(1,3,1);
plot(x_1,T_1)
set(gca,'XTick',(-L:0.030:L));
xlim([-L L]);
title('Por medio de diferencias finitas');
grid on
xlabel('x [m]')
ylabel('T [°C]')
%***********************************************************


%************* Parte 2: CON SOLVER DE MATLAB ***************

%Se utiliza el solver de ecuaciones diferenciales con valores 
%de frontera conocido como "bv4pc".

solinit = bvpinit(linspace(-L,L,100),[140 140]);
sol = bvp4c(@bvp4ode,@bvp4bc,solinit);

subplot(1,3,2);
x_2 = linspace(-L, L);
T_2 = deval(sol,x_2);
plot(x_2,T_2(1,:));

set(gca,'XTick',(-L:0.030:L));
xlim([-L L]);
title('Distribución de Temperatura - Solver de Matlab')
grid on
xlabel('x [m]')
ylabel('T [°C]')
%***********************************************************


%***** Parte 3: MÉTODO ANALÍTICO - SOLUCIÓN EXACTA *********

%Se utiliza la solución obtenida en clase

T_3 = @(x) ((q*L^2)/(2*K))*(1-(x.^2)/(L^2)) + ...
    ((Ts_2 - Ts_1)/2) *(x./L)  +  ((Ts_1 + Ts_2)/2); 

x_3 = (-L : L/10^6 : L);

subplot(1,3,3);
plot(x_3,T_3(x_3))

set(gca,'XTick',(-L:0.030:L));
set(gcf,'color','w')
xlim([-L L]);
title('Distribución de Temperatura -  Sol. Exacta')
grid on
xlabel('x [m]')
ylabel('T [°C]')
%***********************************************************
 