function exp_mat = rodrigues_SO3_from_R3xR(w_R3, t_R)
    % > input: w_R3 \in so(3) coordinate \subset R^{3} unit vector, t_R \in R scalar 
    % - [rodrigues' formula]: axis-angle representation of rotation matrix
    % - so3 ---> SO3
    % < output: mat \in SO(3) rotation matrix R^{3x3}
    
    [ n, k ] = size(w_R3);
    assert( n == 3 , "[ERR] hat_w_so3 must be 3x1 matrix" );
    assert( norm(w_R3) == 1, "[ERR] w_R3 must be in UNIT so(3) coordinate");

    % pre-compute const, (minimize division operations)
    exp_mat = eye(3);
    hat_w_so3 = hat_so3(w_R3);    
    
    % [ rodrigues' formula unit vector x scalar ]:
    exp_mat = exp_mat + sin(t_R) * hat_w_so3 + (1 - cos(t_R)) * hat_w_so3^2;
    
    % [ validation ]:
    if isnumeric(exp_mat) 
        assert(det(exp_mat) == 1, ...
            "[ERR] exp_mat should have a determinant of |R|=1 to be in SO(3)")
    end
end