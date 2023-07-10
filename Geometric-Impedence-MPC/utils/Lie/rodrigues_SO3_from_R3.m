function exp_mat = rodrigues_SO3_from_R3(w_R3)
    % > input: w_R3 \in so(3) coordinate \subset R^{3} vector
    % - [rodrigues' formula]
    % - so3 ---> SO3
    % < output: mat \in SO(3) rotation matrix R^{3x3}
    
    [n,m] = size(w_R3);
    assert( n == 3 && m == 1 , "[ERR] hat_w_so3 must be 3x1 matrix" );

    % pre-compute const, (minimize division operations)
    w_abs = norm(w_R3,2);
    inv_w_abs = 1 / w_abs;
    inv_w_abs_sqr = inv_w_abs * inv_w_abs;
    hat_w_so3 = hat_so3(w_R3);
    
    % [ rodrigues' formula std ]:
    exp_mat = eye(3) + sin(w_abs) * hat_w_so3 * inv_w_abs + (1 - cos(w_abs)) * hat_w_so3^2 * inv_w_abs_sqr;
    
    % [ validation ]:
    if isnumeric(exp_mat) 
        assert(det(exp_mat) == 1, ...
            "[ERR] exp_mat should have a determinant of |R|=1 to be in SO(3)")
    end
end