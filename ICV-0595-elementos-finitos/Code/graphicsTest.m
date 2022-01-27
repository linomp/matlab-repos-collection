close all; clc; clear

% Truss geometry setup for animation
% Truss treated as a directed graph: adjacency list, edges and shit

L = [10 10 10];
a = 10;
nSet = struct('nodes', {cell(11,1)}); 

% Assign node coords
coords = [ [zeros(4,1); (a/2)*ones(3,1); a*ones(4,1)], ... 
           [0; L(1); sum(L(1:2)); sum(L(1:3)); ... 
            sum(L(1:2)) + L(3)/2 ; L(1) + L(2)/2; L(1)/2; ...    
            0; L(1); sum(L(1:2)); sum(L(1:3))] ];

for i = 1:11
    nSet.nodes{i} = struct('x',coords(i,1),'y',coords(i,2));
end

% Adjacency list
lst = {[2,7,8],[3,6,7,9],[4,5,10,6],[5,11], ... 
       [],[],[],[7],[6,7,8],[5,6,9],[5,10]};

%Scan down the adjacency list of each node, draw the edges
hold on
for i = 1:11
    for j = 1:length(lst{i})
        other = nSet.nodes{lst{i}(j)};
        line([nSet.nodes{i}.x, other.x], [nSet.nodes{i}.y, other.y],'linewidth',3);
    end
end
plot(coords(:,1),coords(:,2),'ro','markerfacecolor','r','markersize',8);
axis equal
grid
%tcmgPlotFormat()