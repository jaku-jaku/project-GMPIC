function exp_mat = rodrigues_SO3_from_so3(hat_w_so3)
    % > input: hat_w_so3 \in so(3) \subset R^{3x3} skew-symmetric matrix
    % - [rodrigues' formula]
    % - so3 ---> SO3
    % < output: mat \in SO(3) rotation matrix R^{3x3}
    
    [n,m] = size(hat_w_so3)
    assert( n == 3 && m == 3 , "[ERR] hat_w_so3 must be 3x3 matrix" )

    w_R3 = vee_so3(hat_w_so3);
    w_abs = norm(w_R3,2);
    
    USE_FORMULATION_EFF = false;
    if USE_FORMULATION_EFF
        % [ efficient closed form expression ]: \cite{parkGeometricIntegrationEuclidean2005}
        w_abs_half = w_abs/2;
        s = sin(w_abs_half)/w_abs_half;
        c = cos(w_abs_half);
        
        alpha = s * c;
        beta = s * s;
        
        exp_mat = eye(3) + alpha * hat_w_so3 + 1/2 * beta * hat_w_so3^2;
    else
        % [ rodrigues' formula std ]:
        alpha = sin(w_abs) / w_abs;
        beta = (1 - cos(w_abs)) / (w_abs * w_abs);
        exp_mat = eye(3) + alpha * hat_w_so3 + beta * hat_w_so3^2;
    end
end
