
function mat_SE3 = exp_map_SE3_from_R6xR(xi, t_R)
    [n,k] = size(xi);
    assert(n == 6, "[ERR] ad_se3_: xi \in se(3) coordinates vector R^{6} (v,w)!");
    v_R3 = xi(1:3);
    w_R3 = xi(4:6);
    assert(det(w_R3) == 1, "[ERR] w_R3 should be a unit vector!");
    assert(0 < t_R && t_R < 2*pi, "[ERR] t_R should be in [0, 2*pi]!"); % --> A(.) non-singular
    
    % pre-compute const, (minimize division operations):
    
    if norm(w_R3) == 0
        mat_SE3 = [eye(3), v_R3 * t_R; zeros(1,3), 1]; % (Murray 2.32)
    else
        % [ rodrigues' formula std ]:
        R = rodrigues_SO3_from_R3xR(w_R3, t_R);
    
        % [ transplation ]:
        p = ((eye(3) - R) * hat_so3(w_R3) + w_R3 * w_R3.' * t_R) * v_R3; % (Murray 2.36)
        
        % [ OUTPUT ]:
        mat_SE3 = [R, p; zeros(1,3), 1];
    end
end