
% FEM Final Project
% Prof. Alex Jerves
% 22/12/2017
%
% Contact: lino.mediavilla@estud.usfq.edu.ec

close all; clc; clear all

% Material & section properties
A = ((5.87e-2)*pi^2); %[m^2] cross section of all beams
E = 200e9; % [Pa] young modulus
rho = 7850; % [kg/m^3] density
L = [5 3 3]; % [m] height of each floor
a = 3; % [m] base dimension
T = 10; % [sg] simulation time

% Other problem data
n_dofs = 22;
dt = (1e-3) / 1.5;
omega_1 = (pi*2); % [rad] force 1 freq.
theta_1 = 120 * pi/180; % [rad] force 1 angle
A1 = 1e14; % force 1 amplitude
omega_2 = (pi*2); % [rad] force 2 freq.
theta_2 = 60 * pi/180; % [rad] force 2 angle
A2 = 1e14; % force 2 amplitude

ts = (T/dt) + 1 ; % # de nodos temporales 
t = linspace(0,T,ts);

% Construct elements with:
% node coordinates, name tags, dofs at both ends (the so called "EFT"), 
% elemental K & M in global coordinates, length and angle.
[nodeLists,lst,nSet,coords] = get_the_lists(L, a);
tags = 'A':'V';
truss = struct('elements',{cell(length(tags),1)});
for i = 1:length(tags) 
    % Assign corresponding nodes and dofs to the element
    nodes = cell(2,1);
    EFT = [];
    for j = 1:2
        nodes{j} = nSet.nodes{nodeLists{i}(j)};
        EFT = [EFT nodes{j}.dofs];
    end 
    
    % Get K & M from an element's node coordinates (internally computes
    % angle and element length just based on the coords).
    [k, m, l, ang] = get_Ke_n_Me(nodes);       
    truss.elements{i,1} = struct('tag', tags(i), 'nodes', nodes, ... 
                                'nodeList', nodeLists{i}, 'EFT', EFT, ...
                                'k', k, 'm', m, 'length',l,'angle',ang);
    truss.elements{i,1} = truss.elements{i,1}(1);
end

% Assemble Global K & M
% i.e. sum contributions of each element to each node according to EFT vectors
K = zeros(n_dofs);
M = zeros(n_dofs);
for i = 1:size(truss.elements,1)
    K ([truss.elements{i}.EFT],[truss.elements{i}.EFT]) =  ... 
                   K ([truss.elements{i}.EFT],[truss.elements{i}.EFT]) + ...
                   truss.elements{i}.k;
    M ([truss.elements{i}.EFT],[truss.elements{i}.EFT]) = ... 
                   M ([truss.elements{i}.EFT],[truss.elements{i}.EFT]) + ...
                   truss.elements{i}.m;    
end
K = E * A * K;
M = (rho / E) * M;

%% Solution to Matrix Eqn. system 

u = zeros(n_dofs,ts); % solution vector: as many rows as dofs, as many columns as timesteps

% Initial conditions
u (:,1) = 0; % no initial displacement
u (:,2) = u (:,1); % no initial velocity

% Boundary conditions
u (1:2,:) = zeros(2,ts); % no displacements allowed by bot left support
u (15:16,:) = zeros(2,ts); % no vertical displacement allowed by bot right support

M = M([3:15 , 17:end],[3:15 , 17:end]);
K = K([3:15 , 17:end],[3:15 , 17:end]);

I = -((1/dt^2)*( M + K )); 

for s = 3:ts    
   % Evaluate forces at current time
   f1 = A1*cos(omega_1*t(s));
   f2 = A2*cos(omega_2*t(s));
   f = zeros(22,1);
   f(7) = f1*cos(theta_1);
   f(8) = f1*sin(theta_1);
   f(21) = f2*cos(theta_2); 
   f(22) = f2*cos(theta_2); 
   
   % Solve for the displacements at current timestep
   u([3:15 , 17:end],s) =  I \ ( (M * (1/(dt^2)) * (u([3:15 , 17:end],s-2) ... 
                        - 2 * u([3:15 , 17:end],s-1))) + f([3:15 , 17:end]) ); 
end
 
%% Post process

F = animateTruss(lst,coords,nSet,ts,theta_1,theta_2,u);
 