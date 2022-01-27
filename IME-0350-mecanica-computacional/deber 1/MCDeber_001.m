% ********** Deber 01 de Mecánica Computacional - Módulo 1 *********
clear
clc
% --------------- Derivadas exactas por método analítico -------------

%Creación de la función 
syms x;
f(x) = 3*x^9
%Obtención de las derivadas exactas de la función y
%transformación a funciones anónimas.
df_Ex = matlabFunction(diff(f,x)); %primera derivada
d2f_Ex = matlabFunction(diff(diff(f(x),x),x)); %segunda derivada
% -------------------------------------------------------------------

% --------- Derivadas aproximadas por diferencias centradas -----------

fun = matlabFunction(f(x)); %transformación de f(x) a función anónima

% para almacenar las aproximaciones de la 1era derivada
dfcent = zeros(2,1); 
dfsup = zeros (2,1);
dfinf = zeros(2,1);

% para almacenar las aproximaciones de la 1era derivada
d2fcent = zeros(2,1);
d2fsup = zeros (2,1);
d2finf = zeros(2,1);

h = [0.50; 0.0001]; % valores que serán tomados como el paso
X = 3; % punto donde se evaluarán las derivadas para compararlas

%para almacenar los errores porcentuales (1era derivada)
e1_cent = zeros(2,1);
e1_sup = zeros(2,1); 
e1_inf = zeros(2,1); 

%para almacenar los errores porcentuales (2da derivada)
e2_cent = zeros(2,1); 
e2_sup = zeros(2,1); 
e2_inf = zeros(2,1); 

for i =1:2  
  %1era derivada
  dfsup(i,1)= (-fun(X+2*h(i))+4*fun(X+h(i))-3*fun(X)) / (2*h(i));  %forward
  dfinf(i,1)= (3*fun(X)-4*fun(X-h(i))+fun(X-2*h(i))) / (2*h(i));   %backward
  dfcent(i,1) = (fun(X+h(i))-fun(X-h(i))) / (2*h(i));            %centered 
  
  %2da derivada
  d2fsup(i,1)= (-fun(X+3*h(i))+4*fun(X+2*h(i))-5*fun(X+h(i))... %forward          
                    +2*fun(X)) / ((h(i))^2);
  d2finf(i,1)= (2*fun(X)-5*fun(X-h(i))+4*fun(X-2*h(i))... %backward
                    -fun(X-3*h(i))) / ((h(i))^2);
  d2fcent(i,1) = (fun(X+2*h(i,1))-2*fun(X+h(i,1))+fun(X))...    %centered 
                    / ((h(i,1))^2);  
  
  
  %Error porcentual - 1era derivada      
  e1_cent(i,1) = 100*abs(df_Ex(X)-dfcent(i,1)) ...
                          / df_Ex(X);
  e1_sup(i,1) =  100*abs(df_Ex(X)-dfsup(i,1)) ...
                          / df_Ex(X);                  
  e1_inf(i,1) =  100*abs(df_Ex(X)-dfinf(i,1)) ...
                          / df_Ex(X);                 
             
  
  %Error porcentual - 2da derivada      
  e2_cent(i,1) = 100*abs(d2f_Ex(X)-d2fcent(i,1)) ...
                          / d2f_Ex(X); 
  e2_sup(i,1) =  100*abs(d2f_Ex(X)-d2fsup(i,1)) ...
                          / d2f_Ex(X);
  e2_inf(i,1) =  100*abs(d2f_Ex(X)-d2finf(i,1)) ...
                          / d2f_Ex(X);
end

%RESULTADOS
display('Para x = 3');

%Resultado por diferenciación simbólica (Derivada exacta)
Primera_Derivada_Exacta = df_Ex(X)
Segunda_Derivada_Exacta = d2f_Ex(X)

%Resultados de diferencias centradas
Diferencias_Centradas = table(h, dfcent(:,1),e1_cent,d2fcent(:,1),e2_cent);
Diferencias_Centradas.Properties.VariableNames = {'h','Primera_Derivada','Error_Primera_Derivada','Segunda_Derivada','Error_Segunda_Derivada'}

%Resultados de diferencias hacia adelante
Hacia_Adelante = table(h, dfsup(:,1),e1_sup, d2fsup(:,1),e2_sup);
Hacia_Adelante.Properties.VariableNames = {'h','Primera_Derivada','Error_Primera_Derivada','Segunda_Derivada','Error_Segunda_Derivada'}

%Resultados de diferencias hacia atrás
Hacia_Atras = table(h, dfinf(:,1),e1_inf, d2finf(:,1),e2_inf);
Hacia_Atras.Properties.VariableNames = {'h','Primera_Derivada','Error_Primera_Derivada','Segunda_Derivada','Error_Segunda_Derivada'}

% ----------------------------------------------------------------------