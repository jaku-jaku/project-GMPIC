
function mat_SE3 = exp_map_SE3_from_R6xR(xi, t_R)
    validate.if_dimension(" xi \in se(3) coordinates", xi, [6,1]);
    v_R3 = xi(1:3);
    w_R3 = xi(4:6);
    
    validate.if_unitVector("w_R3", w_R3);
    validate.if_inRangeStrict("t_R", t_R, [0, 2*pi]); % --> A(.) non-singular
    
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