function mat_SE3 = exp_map_SE3_from_R6(xi)
    [n, m] = size(xi);
    assert(n == 6, "[ERR] ad_se3_: xi \in se(3) coordinates vector R^{6xk} (v,w)!");
    v_R3 = xi(1:3,:);
    w_R3 = xi(4:6,:);
    
    % pre-compute const, (minimize division operations)
    w_abs = norm(w_R3,2);
    inv_w_abs = 1 / w_abs;
    inv_w_abs_sqr = inv_w_abs * inv_w_abs;
    hat_w_so3 = hat_so3(w_R3);
    
    % [ rodrigues' formula std ]:
    R = eye(3) + sin(w_abs) * hat_w_so3 * inv_w_abs + (1 - cos(w_abs)) * hat_w_so3^2 * inv_w_abs_sqr;
    
    % [ validation ]:
    if isnumeric(R) 
        assert(det(R) == 1, ...
            "[ERR] rotation matrix should have a determinant of |R|=1.")
    end

    % [ transplation ]:
    p = inv_w_abs_sqr * ((eye(3) - R) * hat_w_so3 + w_R3 * w_R3') * v_R3;
    
    % [ OUTPUT ]:
    mat_SE3 = [R, p; zeros(1,3), 1];
end