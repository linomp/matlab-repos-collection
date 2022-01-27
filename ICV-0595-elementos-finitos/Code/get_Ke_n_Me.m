
function [k, m, Le, a] = get_Ke_n_Me(nodes)
    
    % Define a vector from element's node 1 to element's node 2
    elem_vector = [nodes{1}.x, nodes{1}.y] - [nodes{2}.x, nodes{2}.y];  
    
    % Calculate its norm and angle (ccw wrt to horizontal axis)
    Le = norm(elem_vector);
    a = atan(elem_vector(2)/elem_vector(1));
    
    k = (1/Le) * [  cos(a)*cos(a),   sin(a)*cos(a),  -cos(a)*cos(a),  -sin(a)*sin(a)*cos(a);
                    sin(a)*cos(a),   sin(a)*sin(a),  -sin(a)*cos(a),  -sin(a)*sin(a);
                   -cos(a)*cos(a),  -sin(a)*cos(a),   cos(a)*cos(a),   sin(a)*cos(a);
                   -sin(a)*cos(a),  -sin(a)*sin(a),   sin(a)*cos(a),   sin(a)*sin(a) ];
    
    m = (Le/3) * [  cos(a)*cos(a),        sin(a)*cos(a),        (1/2)*cos(a)*cos(a), (1/2)*sin(a)*sin(a)*cos(a);
                    sin(a)*cos(a),        sin(a)*sin(a),        (1/2)*sin(a)*cos(a), (1/2)*sin(a)*sin(a);
                    (1/2)*cos(a)*cos(a),  (1/2)*sin(a)*cos(a),  cos(a)*cos(a),       sin(a)*cos(a);
                    (1/2)*sin(a)*cos(a),  (1/2)*sin(a)*sin(a),  sin(a)*cos(a),       sin(a)*sin(a) ];
end

