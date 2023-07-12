function xi = twist_se3_coord_from_axis_point(w_R3,q_R3)
    % input: mat \in se(3) \subset R^{4x4} 
    % Example uses: 
    %   - rotation about a line
    [n1, m1] = size(w_R3);
    [n2, m2] = size(q_R3);
    assert(prod([n1,m1,n2,m2] == [3,1,3,1]), "[ERR] w_R3 and q_R3 \in R^3!");
    
    xi = [ - hat_so3(w_R3) * q_R3; w_R3 ];
    % return: xi \in se(3) coordinates vector R^{6}
end
