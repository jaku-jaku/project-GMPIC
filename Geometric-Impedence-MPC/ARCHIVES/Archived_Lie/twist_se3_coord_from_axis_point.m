function xi = twist_se3_coord_from_axis_point(w_R3,q_R3)
    % input: mat \in se(3) \subset R^{4x4} 
    % Example uses: 
    %   - rotation about a line
    validate.if_dimension("w_R3", w_R3, [3,1]);
    validate.if_dimension("q_R3", q_R3, [3,1]);
    
    xi = [ - hat_so3(w_R3) * q_R3; w_R3 ];
    % return: xi \in se(3) coordinates vector R^{6}
end
