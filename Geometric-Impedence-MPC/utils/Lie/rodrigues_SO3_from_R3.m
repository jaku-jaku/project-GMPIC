function exp_mat = rodrigues_SO3_from_R3(w_R3)
    % > input: w_R3 \in so(3) coordinate \subset R^{3} vector
    % - [rodrigues' formula]
    % - so3 ---> SO3
    % < output: mat \in SO(3) rotation matrix R^{3x3}
    
    [n,m] = size(w_R3);
    assert( n == 3 && m == 1 , "[ERR] hat_w_so3 must be 3x1 matrix" );

    w_abs = norm(w_R3,2);
    hat_w_so3 = hat_so3(w_R3);
    % normalized w_so3:
    hat_w_so3 = hat_w_so3 / w_abs;
    
    % [ rodrigues' formula std ]:
    exp_mat = eye(3) + sin(w_abs) * hat_w_so3 + (1 - cos(w_abs)) * hat_w_so3^2;

end