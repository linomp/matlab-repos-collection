
function [nodeLists,lst,nSet,coords] = get_the_lists(L, a)

    % Adjacency list and DOFs associated with each node
    lst = {[2,7,8], [3,6,7,9], [4,5,10,6], [5,11], ... 
           [], [], [], [7], [6,7,8], [5,6,9], [5,10]};
       
    dofsList = {[1,2], [3,4], [5,6], [7,8], [13,14], [11,12], [9,10], ...
            [15,16], [17,18], [19,20], [21,22]};

    % Initial node coords
    coords = [ [zeros(4,1); (a/2)*ones(3,1); a*ones(4,1)], ... 
               [0; L(1); sum(L(1:2)); sum(L(1:3)); ... 
                sum(L(1:2)) + L(3)/2 ; L(1) + L(2)/2; L(1)/2; ...    
                0; L(1); sum(L(1:2)); sum(L(1:3))] ];
    
    % Couple nodes with associated dofs        
    for i = 1:11
        nSet.nodes{i} = struct('x',coords(i,1),'y',coords(i,2));
        nSet.nodes{i}.dofs = dofsList{i}; 
    end       
    
    % Nodes associated with each element - for further use in main script
    nodeLists = {[1,8], [1,7], [2,7], [2,9], [2,6], [3,6], [3,10], [3,5], ...
               [4,5], [5,11], [10,5], [6,10], [9,6], [7,9], [8,7], [1,2], ...
               [2,3],[3,4], [11,4] , [10,11], [9,10], [8,9]};
            
end