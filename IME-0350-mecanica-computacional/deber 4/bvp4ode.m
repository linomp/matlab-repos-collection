function dydx = bvp4ode(x,y)

q = 1.5 * 10 ^6; %[W/m^3]
K = 75; % [W/m*K]

dydx = [ y(2)  -q/K]; 

