
function F = animateTruss(lst, coords, nSet, ts, theta_1, theta_2, u)
    % Truss geometry setup for animation
    % Truss treated as a directed graph: adjacency list, edges and shit
 
    F = struct('cdata',[],'colormap',[]);  
    %axis([-1 4 -1 11]);
     
    % Scan down the adjacency and dof lists of each node: update coordinates
    % and draw the edges
    c = 1;
    for k = 1:100:ts  
          
        % Update node coordinates with corresponding calculated
        % displacements
        for i = 1:11           
            nSet.nodes{i}.x = nSet.nodes{i}.x + u(nSet.nodes{i}.dofs(1),k);
            nSet.nodes{i}.y = nSet.nodes{i}.y + u(nSet.nodes{i}.dofs(2),k);
            coords(i,:) = [nSet.nodes{i}.x, nSet.nodes{i}.y];
        end
         
        % Draw the edges
        grid on; hold on
        for i = 1:11
            for j = 1:length(lst{i})
                other = nSet.nodes{lst{i}(j)};
                h(1) = line([nSet.nodes{i}.x, other.x], [nSet.nodes{i}.y, other.y],'linewidth',5,'color','c');
            end
        end
        h(2) = plot(coords(:,1),coords(:,2),'ro','markerfacecolor','r','markersize',10);
        
        % draw force directions
        f1x = 10*cos(theta_1) + nSet.nodes{4}.x;
        f1y = 10*sin(theta_1) + nSet.nodes{4}.y;
        f2x = 10*cos(theta_2) + nSet.nodes{11}.x;
        f2y = 10*sin(theta_2) + nSet.nodes{11}.y;
        plot([f1x,nSet.nodes{4}.x],[f1y, nSet.nodes{4}.y], 'm-.','linewidth',3);
        plot([f2x,nSet.nodes{11}.x],[f2y, nSet.nodes{11}.y], 'g-.','linewidth',3);
        
        %axis([-4 7 -1 14]);
        axis([-5 8 -1 12]);
        %axis equal
        drawnow 
        tcmgPlotFormat()
        F(c) = getframe(); c = c+1; clf
    end
    
    makeGif(F,'FEM_Truss.gif');
    movie2avi(F,'FEM_Truss.avi','Compression','None')
end