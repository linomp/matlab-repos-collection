% PROBLEMA 2D - TRANSITORIO
% JOMINY
%
% Probeta:
%   |-  2.5 cm -|  (1 in.)
    %%%%%%%%%%%%% -
    %           % |
    %           % |
    %           % 10 cm. (4 in.)
    %           % | 
    %           % |
    %           % |
    %%%%%%%%%%%%% -
    
clear 
clc

%Condiciones del Experimento
h_agua = 10000; % [W/K*m^2]
T_agua = 15; %°C
h_aire = 5; %[W/K*m^2]
T_aire = 25; %°C

%Temperatura inicial de la probeta
To = 925; %[°C]  

%Propiedades del material (AISI 1018)
K = 51.9; % [W/m*K]   
densidad = 7.872*10^3; % [kg/m^3]
Cp = 486 ;% [J/(kg·K)]
alpha = K /(densidad*Cp); %[m^2/s]

%Dimensiones de la Probeta
B = 0.025; % base de la geometría [m]
H = 0.1; % altura de la geometría [m]

%Generación de la red nodal 
%(dx = dy = 0.0050 -> tamaño: 21 fil. x  6 col.) -> dt máx: 0.2346 [sg]  
dx = 0.0025*2;  % m
dy = 0.0025*2;  % m
nX = (B / dx)+1; 
nY = (H / dy)+1;

x = 1000.*(0:dx:B); 
y = 1000.*(H:-dy:0);
[X,Y] = meshgrid(x,y);

%Se generan matrices de unos de 40 filas por 10 columnas 
%donde se almacenarán las temperaturas de los tiempos p y p+1
T_p = ones(nY,nX);
T_p1 = ones(nY,nX);

%Declaración de vectores vacíos para almacenar temperaturas en 
%bordes superior e inferior de la geometría.
bot_edge_Temp = [];
top_edge_Temp = [];

%Declaración de vectores vacíos para almacenar temperaturas en 
%puntos específicos de la geometría.
punto_control_1 = [];
punto_control_2 = [];
punto_control_3 = [];
punto_control_4 = [];
punto_control_5 = [];
punto_control_6 = [];
punto_control_7 = [];
punto_control_8 = [];

%Fo, Bi, Time Step

dt = 0.01; %[sg] -- max: 0.2346 [sg] 
tMax = 10^4; %[sg]  
t = []; %Vector vacío para guardar los tiempos 
p = 1;
maxSteps = uint32(tMax/dt);
tol = 1e-7;

Fo = alpha*dt /(dx^2)
Bi_agua = h_agua*dx/K
Bi_aire = h_aire*dx/K

%Comprobación de criterios de convergencia
if (Fo > 1/4 || Fo*(2+Bi_agua) > 1/2 || Fo*(2+Bi_aire) > 1/2 || Fo*(1+Bi_agua) > 1/4 || Fo*(1+Bi_aire) > 1/4)
    display('NO CUMPLE CRITERIO DE ESTABILIDAD NUMÉRICA');
else
    display('SI CUMPLE CRITERIO DE ESTABILIDAD NUMÉRICA');
end

%%

display('Simulación en curso ...');

while p < maxSteps+1      
    
    %En la primera iteración, setear todos los nodos a temperatura To.
    if p == 1
        T_p1 = To.*T_p; 
    else
        %En cada iteración, se almacenan las temperaturas de la anterior iteración
        % en una matriz de temperaturas en tiempo p.
        T_p = T_p1; 
        
        %Para luego calcular las nuevas temperaturas(tiempo p+1).
        for n = 1:nX
            for m = 1:nY

                %Conducción en nodos interiores
                if(n >=2 && n <= nX-1 && m >= 2 && m <= nY-1)
                    T_p1(m,n) = (Fo*(T_p(m-1,n)+T_p(m,n+1)+T_p(m+1,n)+T_p(m,n-1)))+ ... 
                        (1-4*Fo)*T_p(m,n);  
                
                %Convección del aire en la pared izquierda (sin esquinas)
                elseif(n == 1 && m >= 2 && m <= nY-1)           
                    T_p1(m,n) = Fo*((2*(T_p(m,n+1)))+T_p(m+1,n)+T_p(m-1,n)+(2*Bi_aire*T_aire)) + ...
                        (1-4*Fo-2*Bi_aire*Fo)*T_p(m,n);
                   
                %Convección del aire en la pared derecha (sin esquinas)  
                elseif(n == nX && m >= 2 && m <= nY-1)   
                     T_p1(m,n) = Fo*((2*(T_p(m,n-1)))+T_p(m+1,n)+T_p(m-1,n)+(2*Bi_aire*T_aire))+ ...
                         (1-4*Fo-2*Bi_aire*Fo)*T_p(m,n);
                      
                %Convección del aire en la pared superior (sin esquinas)     
                elseif( m == 1 && n >= 2 && n<= nX-1)           
                     T_p1(m,n) = Fo*((2*(T_p(m+1,n)))+T_p(m,n+1)+T_p(m,n-1)+(2*Bi_aire*T_aire))+ ...
                         (1-4*Fo-2*Bi_aire*Fo)*T_p(m,n);
                       
                %Convección del agua en la pared inferior (sin esquinas)
                elseif( m == nY && n >= 2 && n <= nX-1) 
                     T_p1(m,n) = Fo*((2*(T_p(m-1,n)))+T_p(m,n+1)+T_p(m,n-1)+(2*Bi_agua*T_agua))+ ...
                         (1-4*Fo-2*Bi_agua*Fo)*T_p(m,n);
                     
                %Esquina superior derecha
                elseif(m == 1 && n == nX) 
                    T_p1(m,n) = (2*Fo)*((T_p(m,n-1)+T_p(m+1,n))+(2*Bi_aire*T_aire))+ ...
                         (1-4*Fo-4*Bi_aire*Fo)*T_p(m,n);
                
                %Esquina superior izquierda 
                elseif(m == 1 && n == 1)              
                    T_p1(m,n) = (2*Fo)*((T_p(m,n+1)+T_p(m+1,n))+(2*Bi_aire*T_aire))+ ...
                         (1-4*Fo-4*Bi_aire*Fo)*T_p(m,n);
                     
                %Esquina inferior derecha
                elseif(m == nY &&  n == nX)                           
                    T_p1(m,n) =  (2*Fo)*((T_p(m,n-1)+T_p(m-1,n))+(2*Bi_agua*T_agua))+ ...
                         (1-4*Fo-4*Bi_agua*Fo)*T_p(m,n);
                        
                %Esquina inferior izquierda
                elseif(m == nY  &&  n == 1)                    
                    T_p1(m,n) =  (2*Fo)*((T_p(m,n+1)+T_p(m-1,n))+(2*Bi_agua*T_agua))+ ...
                         (1-4*Fo-4*Bi_agua*Fo)*T_p(m,n);   

                end
            end
        end
    end
    
    t = [t (p*dt)]; %anexar cada tiempo al vector de tiempos.
       
    bot_edge_Temp = [bot_edge_Temp;T_p1(nY,3);];  
    punto_control_1 = [punto_control_1;T_p1(nY-1,3)];
    punto_control_2 = [punto_control_2;T_p1(nY-2,3)];
    punto_control_3 = [punto_control_3;T_p1(nY-3,3)];
    punto_control_4 = [punto_control_4;T_p1(nY-4,3)];    
    punto_control_5 = [punto_control_5;T_p1(nY-5,3)];
    punto_control_6 = [punto_control_6;T_p1(nY-7,3)];
    punto_control_7 = [punto_control_7;T_p1(nY-10,3)];
    punto_control_8 = [punto_control_8;T_p1(nY-15,3)];  
    top_edge_Temp = [top_edge_Temp;T_p1(1,3);]; 
       
    %Graficar cada 350 pasos
     if(p == 1 || uint16(p/350) == p/350)
        
        h1 = subplot(1,2,1);  
        ax1 = get(h1,'position');  
        set(h1,'position',ax1);           
        contourf(X,Y,T_p1,60) 
        title('Jominy End Quench');
        xlabel('X [mm]');
        ylabel('Y [mm]');
        colormap(jet(250))
        caxis( [T_agua To] );
        barra = colorbar;
        set(barra,'YTick',[T_agua:35:To]);
        set(get(barra,'title'),'string','T [°C]');          
    
        subplot(1,2,2);        
        semilogx(t,bot_edge_Temp,'b',t,punto_control_1,'--',...
            t,punto_control_2,'-',t,punto_control_3,'-.',t,punto_control_4,'--',...
            t,punto_control_5,'-',t,punto_control_6,'-.',t,punto_control_7,'--',...
            t,punto_control_8,t,top_edge_Temp,'r');    
       %legend('Y = 0','Y = 5','Y = 10','Y = 15','Y = 20','Y = 25','Y = 35','Y = 50','Y = 75','Y = 100')
        grid on
        title('Temperatures at control points [°C]');
        xlabel('t [sg]');      
         
        drawnow
     end 
    
     %Criterio de parada :: ¿Es la diferencia de temperaturas de un mismo
     %punto en dos puntos temporales consecutivos menor a la tolerancia?
     % Si -> break  
     % No -> seguir iterando
    if length(top_edge_Temp)>=2
         if (abs((top_edge_Temp(length(top_edge_Temp))-top_edge_Temp(length(top_edge_Temp)-1)))<=tol)
            break
         end
    end
    
    p = p + 1;
    
end 


display('Simulación Terminada!');
Tiempo_transcurrido = p*dt

h1 = subplot(1,2,1);  
ax1 = get(h1,'position');  
set(h1,'position',ax1);           
contourf(X,Y,T_p1,60) 
title('Jominy End Quench');
xlabel('X [mm]');
ylabel('Y [mm]');
colormap(jet(250))
caxis( [T_agua To] );
barra = colorbar;
set(barra,'YTick',[T_agua:35:To]);
set(get(barra,'title'),'string','T [°C]');    

subplot(1,2,2);
semilogx(t,bot_edge_Temp,'b',t,punto_control_1,'--',...
            t,punto_control_2,'-',t,punto_control_3,'-.',t,punto_control_4,'--',...
            t,punto_control_5,'-',t,punto_control_6,'-.',t,punto_control_7,'--',...
            t,punto_control_8,t,top_edge_Temp,'r'); 
grid on
legend('y = 0 mm','y = 5 mm','y = 10 mm','y = 15 mm','y = 20 mm','y = 25 mm','y = 35 mm','y = 50 mm','y = 75 mm','y = 100 mm')       
title('Temperatures at control points [°C]');
xlabel('t [sg]');     


 